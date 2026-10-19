//! Ground physics on terrain: slopes, heights and curved ground.

use super::*;
use crate::environment::Environment;
use crate::trajectory::{ShotInput, ShotResult, simulate_flight, simulate_shot};
use crate::units;

const BALL: BallModel = BallModel::TOUR;
const R: f64 = 0.04267 / 2.0;

/// Terrain `y = h(x, z)` with normal from the gradient `(∂h/∂x, ∂h/∂z)`,
/// made of `surface(x, z)`.
fn heightfield(
    h: impl Fn(f64, f64) -> f64 + Copy,
    grad: impl Fn(f64, f64) -> (f64, f64) + Copy,
    surface: impl Fn(f64, f64) -> Surface + Copy,
) -> GroundPhysics<FnTerrain<impl Fn(f64, f64) -> GroundSample + Copy>> {
    GroundPhysics::new(FnTerrain(move |x: f64, z: f64| {
        let (hx, hz) = grad(x, z);
        GroundSample {
            height: h(x, z),
            normal: DVec3::new(-hx, 1.0, -hz).normalize(),
            surface: surface(x, z),
        }
    }))
}

/// The plane through the origin with unit normal `n`, made of `surface`.
fn plane(
    n: DVec3,
    surface: Surface,
) -> GroundPhysics<FnTerrain<impl Fn(f64, f64) -> GroundSample + Copy>> {
    let (a, b) = (-n.x / n.y, -n.z / n.y);
    heightfield(
        move |x, z| a * x + b * z,
        move |_, _| (a, b),
        move |_, _| surface,
    )
}

/// Unit normal of a plane inclined `deg` degrees, falling toward the
/// horizontal direction `azimuth_deg` (0 = downrange, +90 = right).
fn tilt(deg: f64, azimuth_deg: f64) -> DVec3 {
    let (theta, az) = (units::deg_to_rad(deg), units::deg_to_rad(azimuth_deg));
    DVec3::new(theta.sin() * az.cos(), theta.cos(), theta.sin() * az.sin())
}

/// Unit vector down the fall line of a plane with unit normal `n`.
fn fall_line(n: DVec3) -> DVec3 {
    tangential_gravity(n).normalize()
}

fn at_origin(vel: DVec3, spin: DVec3) -> BallState {
    BallState {
        time: 0.0,
        pos: DVec3::ZERO,
        vel,
        spin,
    }
}

/// Mechanical energy per unit mass: translation + rotation + height.
fn energy(vel: DVec3, spin: DVec3, y: f64) -> f64 {
    let k = BALL.inertia_factor;
    0.5 * vel.length_squared() + 0.5 * k * R * R * spin.length_squared() + G * y
}

fn shot(input: &ShotInput, ground: &dyn GroundModel) -> ShotResult {
    simulate_shot(input, &BALL, &Environment::SEA_LEVEL, ground)
}

/// A 7-iron with some curve: 2° right, 500 rpm fade sidespin.
fn curving_seven_iron() -> ShotInput {
    ShotInput {
        launch_azimuth_deg: 2.0,
        sidespin_rpm: 500.0,
        ..ShotInput::seven_iron()
    }
}

/// Rolling terrain: a gentle tilt plus bumps a few metres across, slopes
/// up to about 2°.
fn bumps(x: f64, z: f64) -> f64 {
    0.01 * x - 0.015 * z + 0.12 * (x / 6.0).sin() * (z / 4.0).cos()
}

fn bumps_grad(x: f64, z: f64) -> (f64, f64) {
    (
        0.01 + 0.02 * (x / 6.0).cos() * (z / 4.0).cos(),
        -0.015 - 0.03 * (x / 6.0).sin() * (z / 4.0).sin(),
    )
}

// ══════════════════════════════════════════════════════════════════
//  REST ON SLOPES — hold where tan θ ≤ mu_roll, roll downhill otherwise
// ══════════════════════════════════════════════════════════════════

#[test]
fn ball_at_rest_holds_on_gentle_slopes_and_rolls_down_steep_ones() {
    for surface in [Surface::FAIRWAY, Surface::GREEN] {
        let limit = units::rad_to_deg(surface.mu_roll.atan());
        for (factor, holds) in [(0.8, true), (1.2, false)] {
            let deg = units::rad_to_deg((factor * surface.mu_roll).atan());
            let n = tilt(deg, 30.0);
            let ground = GroundPhysics {
                max_time: 2.0,
                ..plane(n, surface)
            };
            let sample = ground.terrain.sample(0.0, 0.0);
            assert_eq!(sample.holds_at_rest(), holds, "{deg:.2}° vs {limit:.2}°");

            let g = ground.settle(&at_origin(DVec3::ZERO, DVec3::ZERO), &BALL);
            if holds {
                assert_eq!(g.termination, Termination::Rest);
                assert_eq!(g.rest.pos, DVec3::ZERO);
                assert_eq!(g.rest.time, 0.0);
            } else {
                assert_eq!(g.termination, Termination::MaxTime);
                let down = fall_line(n);
                let moved = g.rest.pos;
                assert!(moved.dot(down) > 0.05, "moved {moved:?} on {deg:.2}°");
                assert!((moved - moved.dot(down) * down).length() < 1e-9);
                assert!(g.rest.vel.dot(down) > 0.0);
                assert!(g.points.iter().all(|p| p.phase == Phase::Roll));
            }
        }
    }
}

/// Steps a ball with `GroundPhysics::step` until it holds at rest or `t`
/// passes `t_max`; returns the final state and how many steps ended at zero
/// speed.
fn step_until(
    ground: &GroundPhysics<impl Terrain>,
    start: BallState,
    t_max: f64,
) -> (BallState, usize) {
    let mut s = start;
    let mut stops = 0;
    while s.time < t_max {
        let (next, _) = ground.step(&s, &BALL);
        s = next;
        if s.vel == DVec3::ZERO {
            stops += 1;
            if ground.terrain.sample(s.pos.x, s.pos.z).holds_at_rest() {
                break;
            }
        }
    }
    (s, stops)
}

#[test]
fn ball_rolling_uphill_stops_once_then_rolls_back_or_holds() {
    // Rolling up a slope the ball slows to a stop exactly once. On a slope
    // steeper than atan(mu_roll) it then rolls back down, accelerating
    // away from zero speed without chattering; on a gentler one it stays.
    for (deg, holds) in [(6.0, false), (2.0, true)] {
        let n = tilt(deg, 0.0);
        let ground = plane(n, Surface::GREEN);
        let up = -fall_line(n);
        let v0 = 2.0;
        let start = at_origin(v0 * up, n.cross(v0 * up) / R);
        let (end, stops) = step_until(&ground, start, 6.0);
        assert_eq!(stops, 1, "{deg}°: stopped {stops} times");
        if holds {
            assert_eq!(end.vel, DVec3::ZERO);
            assert!(end.pos.dot(up) > 0.0);
            let g = ground.settle(&start, &BALL);
            assert_eq!(g.termination, Termination::Rest);
            assert!((g.rest.pos - end.pos).length() < 1e-9);
        } else {
            assert!(end.vel.dot(up) < -0.5, "rolling back down: {:?}", end.vel);
            assert!(end.pos.dot(up) < 0.0, "below the start: {:?}", end.pos);
            let g = GroundPhysics {
                max_time: 6.0,
                ..ground
            }
            .settle(&start, &BALL);
            assert_eq!(g.termination, Termination::MaxTime);
        }
    }
}

#[test]
fn long_downhill_reports_max_time() {
    // tan 5° > mu_roll on GREEN: the ball rolls on at its terminal speed
    // (where μ_eff = tan θ, ≈ 1.92 m/s) and never stops.
    let n = tilt(5.0, 0.0);
    let ground = GroundPhysics {
        max_time: 8.0,
        ..plane(n, Surface::GREEN)
    };
    let down = fall_line(n);
    let v0 = 1.9 * down;
    let g = ground.settle(&at_origin(v0, n.cross(v0) / R), &BALL);
    assert_eq!(g.termination, Termination::MaxTime);
    assert!((g.rest.time - 8.0).abs() < 1e-9);
    let mu = 5.0_f64.to_radians().tan() / Surface::GREEN.mu_roll;
    let terminal = Surface::GREEN.roll_v_ref * (mu - 1.0).sqrt();
    assert!((g.rest.vel.length() - terminal).abs() < 0.01);
}

// ══════════════════════════════════════════════════════════════════
//  SLIDE ON A SLOPE — the slide-to-roll event
// ══════════════════════════════════════════════════════════════════

#[test]
fn slide_to_roll_on_a_plane_matches_closed_form() {
    // On the fall line of a plane the slip stays along the fall line, so
    // d|u|/dt = −c ± g sin θ with c = ((1+k)/k)·μ·g·cos θ (+ moving
    // downhill, − uphill), and v(t) = v₀ + (∓μ g cos θ + g sin θ)·t along
    // the line until the slip vanishes.
    let k = BALL.inertia_factor;
    let theta = units::deg_to_rad(4.0);
    let n = tilt(4.0, 0.0);
    let down = fall_line(n);
    let mu = Surface::GREEN.mu_slide;
    let c = (1.0 + k) / k * mu * G * theta.cos();
    for (sign, v0, w) in [(1.0, 3.0, 300.0), (-1.0, 2.5, 250.0)] {
        // Moving along sign·down with backspin w about the across-slope axis.
        let dir = sign * down;
        let axis = dir.cross(n); // backspin axis for travel along `dir`
        let start = at_origin(v0 * dir, w * axis);
        let u0 = slip(start.vel, start.spin, n, R);
        assert!((u0.normalize() - dir).length() < 1e-12);
        let rate = c - sign * G * theta.sin();
        let t_roll = u0.length() / rate;
        let accel = -mu * G * theta.cos() + sign * G * theta.sin();

        let g = plane(n, Surface::GREEN).settle(&start, &BALL);
        let roll = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
        assert!(
            (roll.time - t_roll).abs() < 1e-12,
            "roll at {:.12} s, closed form {t_roll:.12} s",
            roll.time
        );
        let v_end = (v0 + accel * t_roll) * dir;
        assert!(
            (roll.vel - v_end).length() < 1e-12,
            "{:?} vs {v_end:?}",
            roll.vel
        );
        let s_end = v0 * t_roll + 0.5 * accel * t_roll * t_roll;
        assert!((roll.pos - s_end * dir).length() < 1e-12);
        assert!(slip(roll.vel, roll.spin, n, R).length() < 1e-12);
    }
}

#[test]
fn slide_across_a_slope_ends_rolling() {
    // Slip across the fall line: the slope turns the slip direction, so
    // the slide-to-roll time lies between |u₀|/(c + g sin θ) and
    // |u₀|/(c − g sin θ), and the ball leaves without slip.
    let k = BALL.inertia_factor;
    let theta = units::deg_to_rad(4.0);
    let n = tilt(4.0, 90.0);
    let across = n.cross(fall_line(n));
    let start = at_origin(3.0 * across, 250.0 * across.cross(n));
    let u0 = slip(start.vel, start.spin, n, R).length();
    let c = (1.0 + k) / k * Surface::GREEN.mu_slide * G * theta.cos();
    let g = plane(n, Surface::GREEN).settle(&start, &BALL);
    let roll = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
    let (fastest, slowest) = (u0 / (c + G * theta.sin()), u0 / (c - G * theta.sin()));
    assert!(fastest < roll.time && roll.time < slowest);
    assert!(slip(roll.vel, roll.spin, n, R).length() < 1e-12);
    assert!(roll.vel.dot(n).abs() < 1e-12);
    // It drifted downhill while sliding.
    assert!(roll.pos.dot(fall_line(n)) > 0.0);
}

// ══════════════════════════════════════════════════════════════════
//  SPIN ABOUT THE VERTICAL ON A SLOPE
// ══════════════════════════════════════════════════════════════════

#[test]
fn spin_about_vertical_kicks_on_a_slope_only() {
    // ω_y has a component in a sloped contact plane, so it slips the
    // contact point and kicks the ball. Rigid gripping ground: the kick is
    // Δv = (k/(1+k))·R·ω_y·(ŷ × n), zero on level ground. On a side slope
    // it acts along the line of play; on an up/downhill slope, across it.
    let k = BALL.inertia_factor;
    let rigid = Surface {
        restitution: Restitution::Constant(0.5),
        mu_slide: 10.0,
        compliance: 0.0,
        ..Surface::GREEN
    };
    let vel = DVec3::new(15.0, -10.0, 1.0);
    let w = 150.0;
    for (n, kicks) in [
        (DVec3::Y, false),
        (tilt(8.0, 90.0), true),
        (tilt(8.0, 0.0), true),
        (tilt(8.0, -135.0), true),
    ] {
        let at = GroundSample {
            height: 0.0,
            normal: n,
            surface: rigid,
        };
        let plain = impact(&at_origin(vel, DVec3::ZERO), &at, &BALL);
        let spun = impact(&at_origin(vel, w * DVec3::Y), &at, &BALL);
        let kick = spun.vel - plain.vel;
        let expected = k / (1.0 + k) * R * w * DVec3::Y.cross(n);
        assert!((kick - expected).length() < 1e-12, "{n:?}: {kick:?}");
        assert_eq!(kick.length() > 0.01, kicks);
    }
}

// ══════════════════════════════════════════════════════════════════
//  FULL SHOTS ON TERRAIN — landing, uphill/downhill, elevated green
// ══════════════════════════════════════════════════════════════════

#[test]
fn uphill_landing_checks_up_downhill_releases() {
    // Planes through the tee, rising or falling 2° downrange. Uphill the
    // ball lands sooner and runs less; downhill it lands later and runs
    // further.
    let input = ShotInput::seven_iron();
    let run = |deg: f64| {
        let r = shot(&input, &plane(tilt(deg, 180.0), Surface::FAIRWAY));
        assert_eq!(r.ground.termination, Termination::Rest);
        (r.flight.carry_m, r.final_pos.x - r.ground.contacts[0].pos.x)
    };
    let (up_carry, up_run) = run(2.0);
    let (flat_carry, flat_run) = run(0.0);
    let (down_carry, down_run) = run(-2.0);
    assert!(up_carry < flat_carry && flat_carry < down_carry);
    assert!(
        up_run < flat_run && flat_run < down_run,
        "runs: up {up_run:.2}, flat {flat_run:.2}, down {down_run:.2} m"
    );
}

#[test]
fn elevated_green_shortens_carry_and_lands_on_it() {
    // A green raised 3 m on a plateau from 150 m to 200 m downrange, with
    // smooth 10 m ramps, fairway elsewhere.
    const H: f64 = 3.0;
    let ramp = |t: f64| t * t * (3.0 - 2.0 * t);
    let ramp_slope = |t: f64| 6.0 * t * (1.0 - t);
    let h = move |x: f64, _z: f64| match x {
        x if (140.0..150.0).contains(&x) => H * ramp((x - 140.0) / 10.0),
        x if (150.0..=200.0).contains(&x) => H,
        x if (200.0..210.0).contains(&x) => H * ramp((210.0 - x) / 10.0),
        _ => 0.0,
    };
    let grad = move |x: f64, _z: f64| match x {
        x if (140.0..150.0).contains(&x) => (H * ramp_slope((x - 140.0) / 10.0) / 10.0, 0.0),
        x if (200.0..210.0).contains(&x) => (-H * ramp_slope((210.0 - x) / 10.0) / 10.0, 0.0),
        _ => (0.0, 0.0),
    };
    let surface = move |x: f64, _z: f64| {
        if (150.0..=200.0).contains(&x) {
            Surface::GREEN
        } else {
            Surface::FAIRWAY
        }
    };
    let ground = heightfield(h, grad, surface);

    let input = ShotInput::seven_iron();
    let flat = simulate_flight(&input, &BALL, &Environment::SEA_LEVEL);
    let r = shot(&input, &ground);
    let landing = r.flight.points.last().unwrap().pos;
    assert!(
        r.flight.carry_m < flat.carry_m - 1.0,
        "carry {:.2} vs flat {:.2}",
        r.flight.carry_m,
        flat.carry_m
    );
    assert!(
        (150.0..=200.0).contains(&landing.x),
        "landed at {landing:?}"
    );
    assert_eq!(landing.y, H);
    assert!((r.flight.carry_m - landing.x.hypot(landing.z)).abs() < 1e-12);
    let first = &r.ground.contacts[0];
    assert_eq!(first.surface, Surface::GREEN);
    assert_eq!(first.normal, DVec3::Y);
    assert_eq!(r.ground.termination, Termination::Rest);
    assert_eq!(r.final_pos.y, H);
}

#[test]
fn tee_sits_on_the_terrain() {
    let ground = heightfield(
        |x, _| 2.0 + 0.01 * x,
        |_, _| (0.01, 0.0),
        |_, _| Surface::FAIRWAY,
    );
    let r = shot(&ShotInput::driver(), &ground);
    let tee = r.flight.points[0].pos;
    assert_eq!(tee, DVec3::new(0.0, 2.0, 0.0));
    let top = r
        .flight
        .points
        .iter()
        .map(|p| p.pos.y)
        .fold(f64::MIN, f64::max);
    let above_tee = top - 2.0;
    assert!(r.flight.apex_m >= above_tee && r.flight.apex_m - above_tee < 0.01);
}

// ══════════════════════════════════════════════════════════════════
//  CROSSINGS ON CURVED TERRAIN — landing and hop arrivals
// ══════════════════════════════════════════════════════════════════

#[test]
fn landing_root_find_is_on_curved_terrain() {
    // Strongly curved terrain and segments far longer than a flight step:
    // the straight-line interpolation of the end clearances is far off,
    // and the root-find must still put the crossing on the surface.
    let wavy = |x: f64, z: f64| 0.5 * (3.0 * x).sin() + 0.2 * (2.0 * z).cos();
    let ground = heightfield(wavy, |_, _| (0.0, 0.0), |_, _| Surface::GREEN);
    let mut worst_linear = 0.0_f64;
    for i in 0..200 {
        let x0 = f64::from(i) * 0.137;
        let z0 = f64::from(i) * -0.071;
        let a = BallState {
            time: 1.0,
            pos: DVec3::new(x0, wavy(x0, z0) + 0.3, z0),
            vel: DVec3::new(8.0, -6.0, 1.0),
            spin: DVec3::new(0.0, -20.0, 300.0),
        };
        let b_pos = a.pos + DVec3::new(0.4, -0.6, 0.05);
        let b = BallState {
            time: 1.05,
            pos: DVec3::new(b_pos.x, b_pos.y.min(wavy(b_pos.x, b_pos.z) - 0.01), b_pos.z),
            vel: DVec3::new(8.0, -6.5, 1.0),
            spin: DVec3::new(0.0, -19.0, 299.0),
        };
        let l = ground.landing(&a, &b).expect("segment crosses the surface");
        let s = (l.time - a.time) / (b.time - a.time);
        let on_segment = a.pos.lerp(b.pos, s);
        let gap = on_segment.y - wavy(on_segment.x, on_segment.z);
        assert!(gap.abs() < 1e-9, "clearance {gap:.3e} at s = {s}");
        assert!((l.pos.x - on_segment.x).abs() < 1e-12 && (l.pos.z - on_segment.z).abs() < 1e-12);
        assert_eq!(l.pos.y, wavy(l.pos.x, l.pos.z));

        let fa = a.pos.y - wavy(a.pos.x, a.pos.z);
        let fb = b.pos.y - wavy(b.pos.x, b.pos.z);
        let p = a.pos.lerp(b.pos, fa / (fa - fb));
        worst_linear = worst_linear.max((p.y - wavy(p.x, p.z)).abs());
    }
    assert!(
        worst_linear > 1e-3,
        "linear guess already exact: {worst_linear:.3e}"
    );
}

#[test]
fn flight_and_hops_land_on_curved_terrain() {
    // The flight's landing and every hop arrival lie on the surface, and
    // each arrival is on the gravity parabola from the previous rebound
    // (the crossing is solved, not stepped). On the ground the ball stays
    // on the surface with its velocity in the tangent plane.
    let ground = heightfield(bumps, bumps_grad, |_, _| Surface::FAIRWAY);
    let r = shot(&ShotInput::driver(), &ground);
    let landing = r.flight.points.last().unwrap();
    assert_eq!(landing.pos.y, bumps(landing.pos.x, landing.pos.z));
    assert!(r.ground.contacts.len() > 5);
    for c in &r.ground.contacts {
        assert_eq!(c.pos.y, bumps(c.pos.x, c.pos.z));
    }
    for w in r.ground.contacts.windows(2) {
        let (prev, next) = (&w[0], &w[1]);
        let tau = next.time - prev.time;
        let on_parabola = prev.pos + prev.vel_out * tau + 0.5 * GRAVITY * tau * tau;
        assert!(
            (on_parabola - next.pos).length() < 1e-9,
            "arrival off the parabola by {:.3e} m",
            (on_parabola - next.pos).length()
        );
    }
    let mut on_ground = 0;
    for p in r.ground.points.iter().filter(|p| p.phase != Phase::Hop) {
        let (hx, hz) = bumps_grad(p.pos.x, p.pos.z);
        let n = DVec3::new(-hx, 1.0, -hz).normalize();
        assert_eq!(p.pos.y, bumps(p.pos.x, p.pos.z));
        assert!(p.vel.dot(n).abs() < 1e-12, "v·n = {:.3e}", p.vel.dot(n));
        on_ground += 1;
    }
    assert!(on_ground > 50);
    assert_eq!(r.ground.termination, Termination::Rest);
}

// ══════════════════════════════════════════════════════════════════
//  INVARIANTS ON TERRAIN — mirror, rotation, energy
// ══════════════════════════════════════════════════════════════════

#[test]
fn mirrored_terrain_and_shot_give_a_mirrored_result() {
    let ground = heightfield(bumps, bumps_grad, |_, _| Surface::FAIRWAY);
    let mirrored = heightfield(
        |x, z| bumps(x, -z),
        |x, z| {
            let (hx, hz) = bumps_grad(x, -z);
            (hx, -hz)
        },
        |_, _| Surface::FAIRWAY,
    );
    let input = curving_seven_iron();
    let flipped = ShotInput {
        launch_azimuth_deg: -input.launch_azimuth_deg,
        sidespin_rpm: -input.sidespin_rpm,
        ..input
    };
    let a = shot(&input, &ground);
    let b = shot(&flipped, &mirrored);
    let mirror = |v: DVec3| DVec3::new(v.x, v.y, -v.z);
    assert_eq!(a.ground.contacts.len(), b.ground.contacts.len());
    for (ca, cb) in a.ground.contacts.iter().zip(&b.ground.contacts) {
        assert!((mirror(ca.pos) - cb.pos).length() < 1e-9);
        assert!((mirror(ca.vel_out) - cb.vel_out).length() < 1e-9);
    }
    assert!((mirror(a.final_pos) - b.final_pos).length() < 1e-9);
    assert!(a.final_pos.z.abs() > 1.0);
    assert_eq!(a.ground.termination, b.ground.termination);
}

#[test]
fn rotating_shot_and_sloped_plane_rotates_the_result() {
    // Rotating the launch direction and a sloped plane together about Y
    // rotates the whole result. Positive azimuth ψ maps +X to
    // (cos ψ, 0, sin ψ).
    let rotate = |v: DVec3, psi: f64| {
        let (s, c) = psi.sin_cos();
        DVec3::new(c * v.x - s * v.z, v.y, s * v.x + c * v.z)
    };
    let n = tilt(2.0, 140.0);
    let input = curving_seven_iron();
    let base = shot(&input, &plane(n, Surface::FAIRWAY));
    for deg in [25.0, -70.0, 180.0] {
        let psi = units::deg_to_rad(deg);
        let turned = ShotInput {
            launch_azimuth_deg: input.launch_azimuth_deg + deg,
            ..input
        };
        let r = shot(&turned, &plane(rotate(n, psi), Surface::FAIRWAY));
        assert_eq!(r.ground.contacts.len(), base.ground.contacts.len());
        let expected = rotate(base.final_pos, psi);
        assert!(
            (r.final_pos - expected).length() < 1e-6,
            "{deg}°: {:?} vs {expected:?}",
            r.final_pos
        );
        assert_eq!(r.ground.termination, base.ground.termination);
    }
}

#[test]
fn energy_never_increases_on_slopes() {
    // Translation + rotation + height never grows at an impact or from one
    // ground point to the next, on planes tilted in several directions.
    let starts = [
        (DVec3::new(21.3, -18.5, 0.0), DVec3::new(0.0, 0.0, 212.0)),
        (DVec3::new(21.0, -18.0, 2.5), DVec3::new(0.0, -50.0, 210.0)),
        (DVec3::new(9.0, -20.0, -1.0), DVec3::new(10.0, 40.0, 700.0)),
        (DVec3::new(3.0, -0.5, 0.5), DVec3::new(-20.0, 30.0, 300.0)),
    ];
    let planes = [
        tilt(2.0, 180.0),
        tilt(2.0, 0.0),
        tilt(5.0, 90.0),
        tilt(3.0, -60.0),
        tilt(6.0, 30.0),
    ];
    for n in planes {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let ground = GroundPhysics {
                max_time: 10.0,
                ..plane(n, surface)
            };
            for (vel, spin) in starts {
                let g = ground.settle(&at_origin(vel, spin), &BALL);
                for c in &g.contacts {
                    let before = energy(c.vel_in, c.spin_in, c.pos.y);
                    let after = energy(c.vel_out, c.spin_out, c.pos.y);
                    assert!(after <= before + 1e-12 * before.abs().max(1.0));
                }
                let mut prev = energy(vel, spin, 0.0);
                for p in &g.points {
                    let e = energy(p.vel, p.spin, p.pos.y);
                    assert!(
                        e <= prev + 1e-9 * prev.abs().max(1.0),
                        "energy rose {prev:.9} → {e:.9} at t={:.3} ({:?}) on {n:?}",
                        p.time,
                        p.phase
                    );
                    prev = e;
                }
            }
        }
    }
}

#[test]
fn step_keeps_the_ball_on_a_slope() {
    let n = tilt(4.0, 45.0);
    let ground = plane(n, Surface::FAIRWAY);
    let across = n.cross(fall_line(n));
    let mut s = at_origin(2.0 * across, DVec3::new(10.0, 30.0, -50.0));
    let mut phases = Vec::new();
    for _ in 0..400 {
        let (next, phase) = ground.step(&s, &BALL);
        s = next;
        phases.push(phase);
        assert!(
            (s.pos.dot(n)).abs() < 1e-12,
            "off the plane by {:.3e}",
            s.pos.dot(n)
        );
        assert!(s.vel.dot(n).abs() < 1e-12);
    }
    assert!(phases.contains(&Phase::Slide) && phases.contains(&Phase::Roll));
}

// ══════════════════════════════════════════════════════════════════
//  HOPS ON SLOPES — the hop train dies out
// ══════════════════════════════════════════════════════════════════

/// Run past the landing point (yards) and termination of `input` on a
/// plane through the tee with unit normal `n`.
fn run_on(input: &ShotInput, n: DVec3, surface: Surface) -> (f64, f64, Termination) {
    let r = shot(input, &plane(n, surface));
    let landing = r.ground.contacts[0].pos;
    (
        units::meters_to_yards(r.final_pos.x - landing.x),
        units::meters_to_yards(r.final_pos.z - landing.z),
        r.ground.termination,
    )
}

#[test]
fn hops_die_out_only_where_the_slope_hands_back_less_than_the_impact_takes() {
    let state = |vel: DVec3| at_origin(vel, DVec3::ZERO);
    let h = 0.005;
    // Level ground: always, even for an elastic rebound.
    let level = DVec3::Y;
    assert!(hops_die_out(
        &state(DVec3::new(2.0, -0.5, 0.0)),
        &state(DVec3::new(2.0, 0.5, 0.0)),
        level,
        &BALL,
        h
    ));
    // A rebound up the slope lands higher: the slope returns nothing.
    let n = tilt(4.0, 0.0);
    let up = -fall_line(n);
    let rebound = 2.0 * up + 0.5 * n;
    assert!(hops_die_out(
        &state(2.0 * up - 0.5 * n),
        &state(rebound),
        n,
        &BALL,
        h
    ));
    // Down the slope: lossless or barely lossy impacts are sustained by
    // the slope; a lossy one dies out.
    let down = fall_line(n);
    let rebound = 2.0 * down + 0.5 * n;
    let w = rebound.dot(n);
    let t = 2.0 * w / (G * n.y);
    let released = -G * (rebound.y * t - 0.5 * G * t * t);
    assert!(released > 0.0);
    let arriving = |extra: f64| {
        let e = 0.5 * rebound.length_squared() + extra;
        state((2.0 * e).sqrt() * (2.0 * down - 0.5 * n).normalize())
    };
    assert!(!hops_die_out(&arriving(0.0), &state(rebound), n, &BALL, h));
    let threshold = released + G * h;
    assert!(!hops_die_out(
        &arriving(0.99 * threshold),
        &state(rebound),
        n,
        &BALL,
        h
    ));
    assert!(hops_die_out(
        &arriving(1.01 * threshold),
        &state(rebound),
        n,
        &BALL,
        h
    ));
}

#[test]
fn fairway_downslopes_end_at_rest() {
    // Up to 6° downhill (below FAIRWAY's 6.8° hold angle) every preset
    // comes to rest, and the driver's run grows with the slope.
    let mut last_driver = 0.0;
    for step in 0..=12 {
        let deg = 0.5 * f64::from(step);
        let n = tilt(deg, 0.0);
        for (name, input) in [
            ("driver", ShotInput::driver()),
            ("7-iron", ShotInput::seven_iron()),
            ("wedge", ShotInput::wedge()),
        ] {
            let (run, _, end) = run_on(&input, n, Surface::FAIRWAY);
            assert_eq!(end, Termination::Rest, "{name} on {deg}°: run {run:.1} yd");
            if name == "driver" {
                assert!(
                    run >= last_driver,
                    "driver {deg}°: {run:.1} < {last_driver:.1} yd"
                );
                last_driver = run;
            }
        }
    }
    assert!(last_driver < 120.0, "driver on 6°: {last_driver:.1} yd");
}

#[test]
fn green_downslopes_end_at_rest() {
    for deg in [1.0, 2.0, 2.5] {
        for input in [
            ShotInput::driver(),
            ShotInput::seven_iron(),
            ShotInput::wedge(),
        ] {
            let (run, _, end) = run_on(&input, tilt(deg, 0.0), Surface::GREEN);
            assert_eq!(end, Termination::Rest, "{deg}°: run {run:.1} yd");
        }
    }
}

#[test]
fn side_slope_rests_mirrored() {
    // 5° falling right or left across the line: the ball drifts downhill,
    // rests, and the two sides mirror exactly.
    for input in [
        ShotInput::driver(),
        ShotInput::seven_iron(),
        ShotInput::wedge(),
    ] {
        let (run_r, lat_r, end_r) = run_on(&input, tilt(5.0, 90.0), Surface::FAIRWAY);
        let (run_l, lat_l, end_l) = run_on(&input, tilt(5.0, -90.0), Surface::FAIRWAY);
        assert_eq!((end_r, end_l), (Termination::Rest, Termination::Rest));
        assert!(lat_r > 1.0 && lat_r < 30.0, "lateral run {lat_r:.1} yd");
        assert_eq!(lat_l, -lat_r);
        assert_eq!(run_l, run_r);
    }
}

// ══════════════════════════════════════════════════════════════════
//  FIRST CROSSING — thin walls and ridges inside one step
// ══════════════════════════════════════════════════════════════════

/// Fairway with a 20 cm wall over `x0 ≤ x < x0 + 0.1`, level elsewhere.
fn wall(x0: f64) -> GroundPhysics<FnTerrain<impl Fn(f64, f64) -> GroundSample + Copy>> {
    heightfield(
        move |x, _| {
            if (x0..x0 + 0.1).contains(&x) {
                0.2
            } else {
                0.0
            }
        },
        |_, _| (0.0, 0.0),
        |_, _| Surface::FAIRWAY,
    )
}

fn segment(from: DVec3, to: DVec3) -> (BallState, BallState) {
    let vel = to - from;
    (
        BallState {
            time: 0.0,
            pos: from,
            vel,
            spin: DVec3::ZERO,
        },
        BallState {
            time: 1.0,
            pos: to,
            vel,
            spin: DVec3::ZERO,
        },
    )
}

#[test]
fn landing_finds_the_first_crossing_of_a_wall_and_a_ridge() {
    // The segment ends below the ground at x = 0.5, but first meets a wall
    // or ridge around x = 0.3 that its ends never see.
    let (a, b) = segment(DVec3::new(0.0, 0.05, 0.0), DVec3::new(1.0, -0.05, 0.0));

    let l = wall(0.3).landing(&a, &b).expect("hits the wall");
    assert!((l.pos.x - 0.3).abs() < 1e-9, "landed at {:?}", l.pos);
    assert_eq!(l.pos.y, 0.2);

    let ridge = |x: f64, _z: f64| 0.1 * (-((x - 0.35) / 0.03).powi(2)).exp();
    let ridge_grad = |x: f64, z: f64| (-2.0 * (x - 0.35) / (0.03 * 0.03) * ridge(x, z), 0.0);
    let ground = heightfield(ridge, ridge_grad, |_, _| Surface::FAIRWAY);
    let l = ground.landing(&a, &b).expect("hits the ridge");
    assert!(l.pos.x < 0.35, "landed at {:?}", l.pos);
    let s = l.time;
    let on_segment = a.pos.lerp(b.pos, s);
    assert!((on_segment.y - ridge(on_segment.x, 0.0)).abs() < 1e-9);
    // Nothing earlier on the segment is below the ridge.
    for i in 0..1000 {
        let p = a.pos.lerp(b.pos, s * f64::from(i) / 1000.0);
        assert!(p.y > ridge(p.x, 0.0));
    }
}

#[test]
fn landing_needs_a_segment_starting_above_the_ground() {
    let ground = plane(tilt(3.0, 0.0), Surface::FAIRWAY);
    let (a, b) = segment(DVec3::new(0.0, -0.01, 0.0), DVec3::new(0.1, -0.2, 0.0));
    assert!(ground.landing(&a, &b).is_none());
    let (a, b) = segment(DVec3::new(0.0, -0.01, 0.0), DVec3::new(0.1, 0.2, 0.0));
    assert!(ground.landing(&a, &b).is_none());
    // Above the ground at every sub-step: no landing.
    let (a, b) = segment(DVec3::new(0.0, 0.05, 0.0), DVec3::new(1.0, 0.01, 0.0));
    assert!(Surface::FAIRWAY.landing(&a, &b).is_none());
}

#[test]
fn hop_finds_a_wall_inside_a_step() {
    // 50 ms steps carry the hop 0.5 m per step, far past a 10 cm wall
    // standing 0.2 m downrange; the ball must still strike it.
    let rigid = Surface {
        restitution: Restitution::Constant(0.8),
        compliance: 0.0,
        ..Surface::FAIRWAY
    };
    let ground = GroundPhysics {
        dt: 0.05,
        ..heightfield(
            |x, _| if (0.2..0.3).contains(&x) { 0.2 } else { 0.0 },
            |_, _| (0.0, 0.0),
            move |_, _| rigid,
        )
    };
    let g = ground.settle(&at_origin(DVec3::new(10.0, -1.0, 0.0), DVec3::ZERO), &BALL);
    assert!(g.contacts.len() >= 2);
    let hit = &g.contacts[1];
    assert!(
        (hit.pos.x - 0.2).abs() < 1e-9,
        "second contact at {:?}",
        hit.pos
    );
    assert_eq!(hit.pos.y, 0.2);
}

// ══════════════════════════════════════════════════════════════════
//  SURFACE BOUNDARIES AND CURVED GROUND ON THE GROUND
// ══════════════════════════════════════════════════════════════════

/// Positions at the start and end of the first step that brings the ball
/// to a stop, stepping from `start` with `GroundPhysics::step`.
fn first_stop(ground: &GroundPhysics<impl Terrain>, start: BallState) -> (DVec3, DVec3) {
    let mut s = start;
    for _ in 0..100_000 {
        let (next, _) = ground.step(&s, &BALL);
        if next.vel == DVec3::ZERO {
            return (s.pos, next.pos);
        }
        s = next;
    }
    panic!("never stopped");
}

#[test]
fn rest_is_judged_where_the_ball_stops() {
    // A 4° upslope: FAIRWAY holds a ball there (tan 4° < 0.12), GREEN does
    // not (tan 4° > 0.055). The ball rolls up and stops one step across a
    // surface boundary; whether it stays is the surface it stopped on.
    let n = tilt(4.0, 180.0);
    let up = -fall_line(n);
    let v0 = 1.5;
    let start = at_origin(v0 * up, n.cross(v0 * up) / R);
    let split = |below: Surface, above: Surface, at: f64| {
        let (a, b) = (-n.x / n.y, -n.z / n.y);
        heightfield(
            move |x, z| a * x + b * z,
            move |_, _| (a, b),
            move |x, _| if x < at { below } else { above },
        )
    };

    // FAIRWAY below, GREEN above: the ball stops on GREEN and rolls back
    // down onto FAIRWAY, where it rests.
    let (before, stop) = first_stop(&plane(n, Surface::FAIRWAY), start);
    let boundary = 0.5 * (before.x + stop.x);
    assert!(before.x < boundary && boundary < stop.x);
    let g = split(Surface::FAIRWAY, Surface::GREEN, boundary).settle(&start, &BALL);
    assert_eq!(g.termination, Termination::Rest);
    assert!(g.rest.pos.x < boundary, "rested at {:?}", g.rest.pos);
    assert_eq!(g.rest.vel, DVec3::ZERO);

    // GREEN below, FAIRWAY above: the ball stops on FAIRWAY and stays.
    let (before, stop) = first_stop(&plane(n, Surface::GREEN), start);
    let boundary = 0.5 * (before.x + stop.x);
    let g = split(Surface::GREEN, Surface::FAIRWAY, boundary).settle(&start, &BALL);
    assert_eq!(g.termination, Termination::Rest);
    assert!(
        (g.rest.pos - stop).length() < 1e-12,
        "{:?} vs {stop:?}",
        g.rest.pos
    );
}

#[test]
fn slide_ends_rolling_on_curved_ground_without_flicker() {
    // A ball sliding over curved ground goes Slide → Roll once: the step
    // that ends the slide leaves it locked to rolling on the surface there,
    // and every slide step but the last is a full step. Its spin about the
    // normal makes this sensitive: where the normal turns from one step to
    // the next, part of that spin lies in the new contact plane and slips
    // the contact point unless the spin is re-locked.
    let ground = heightfield(bumps, bumps_grad, |_, _| Surface::FAIRWAY);
    let (x, z) = (3.0, 1.0);
    let (hx, hz) = bumps_grad(x, z);
    let n = DVec3::new(-hx, 1.0, -hz).normalize();
    let dir = (DVec3::X - DVec3::X.dot(n) * n).normalize();
    let start = BallState {
        time: 0.0,
        pos: DVec3::new(x, bumps(x, z), z),
        vel: 3.0 * dir,
        spin: 40.0 * n,
    };
    let mut s = start;
    let mut phases = Vec::new();
    while s.time < 1.0 {
        let (next, phase) = ground.step(&s, &BALL);
        phases.push((phase, next.time - s.time));
        s = next;
    }
    let slides = phases
        .iter()
        .take_while(|(p, _)| *p == Phase::Slide)
        .count();
    assert!(slides > 10);
    assert!(phases[slides..].iter().all(|(p, _)| *p == Phase::Roll));
    assert!(
        phases[..slides - 1]
            .iter()
            .all(|(_, h)| (h - ground.dt).abs() < 1e-12)
    );

    // The same in `settle`: one switch from Slide to Roll.
    let g = ground.settle(&start, &BALL);
    let switches = g
        .points
        .windows(2)
        .filter(|w| w[0].phase != w[1].phase)
        .count();
    assert_eq!(g.points[0].phase, Phase::Slide);
    assert_eq!(
        switches,
        1,
        "{:?}",
        g.points.iter().map(|p| p.phase).collect::<Vec<_>>()
    );
    assert_eq!(g.termination, Termination::Rest);
}
