use super::regression_data::REFERENCES;
use super::*;
use crate::aero::Tour;
use crate::air::Wind;
use crate::units;

const BALL: Tour = Tour::SEA_LEVEL;
const R: f64 = 0.04267 / 2.0;

/// Gravity vector (m/s²).
const GRAVITY: Vec3 = Vec3::new(0.0, -G, 0.0);

static NO_AIR: Tour = Tour::new(0.0, Wind::None);
static STILL_AIR: Tour = Tour::SEA_LEVEL;

/// No air: ballistic hops, gravity only.
const VACUUM: &Tour = &NO_AIR;

/// Still air at sea level.
fn air() -> &'static Tour {
    &STILL_AIR
}

fn rpm_to_rads(rpm: f64) -> f64 {
    units::rpm_to_rads(rpm)
}

fn state(vel: Vec3, spin: Vec3) -> BallState {
    BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel,
        spin,
    }
}

/// Spin about the across-track axis for a ball moving along +X:
/// positive = backspin.
fn backspin(w: f64) -> Vec3 {
    Vec3::new(0.0, 0.0, w)
}

/// Impact on flat `surface` with incoming `vel` and spin.
fn hit(vel: Vec3, spin: Vec3, surface: &Surface) -> BallState {
    impact(&state(vel, spin), &surface.sample(0.0, 0.0), &BALL)
}

/// Mechanical energy per unit mass: translation + rotation + height.
fn energy(vel: Vec3, spin: Vec3, y: f64) -> f64 {
    let k = BALL.inertia_factor();
    0.5 * vel.length_squared() + 0.5 * k * R * R * spin.length_squared() + G * y
}

/// Ground phase only: a ball already on flat ground moving horizontally.
fn roll(vel: Vec3, spin: Vec3, surface: &Surface) -> GroundResult {
    surface.settle(&state(vel, spin), air())
}

/// Rolling spin for a ball moving along +X at `v` (no slip).
fn rolling(v: f64) -> Vec3 {
    Vec3::new(0.0, 0.0, -v / R)
}

/// Asserts that a ground result stopped bouncing by the 5 mm hop rule: every
/// rebound but the last rises at least 5 mm, the last leaves along the
/// ground, and the ball came to rest.
fn assert_ends_by_hop_height(g: &GroundResult) {
    assert_eq!(g.termination, Termination::Rest);
    let (last, rest) = g.contacts.split_last().unwrap();
    for c in rest {
        let hop = c.vel_out.y * c.vel_out.y / (2.0 * G);
        assert!(hop >= 0.005, "hop of {:.2} mm", hop * 1000.0);
    }
    assert!(last.vel_out.y.abs() < 1e-12);
}

// ══════════════════════════════════════════════════════════════════
//  REGRESSION — the 3D impact reduces to the planar Penner solution
// ══════════════════════════════════════════════════════════════════

#[test]
fn planar_impacts_match_reference() {
    // Zero sidespin, zero azimuth, flat ground, ballistic hops: the slip is
    // along the track and the spin has no component about the tilted
    // normal, so the vector impact must reproduce the planar (scalar)
    // solution to float tolerance, impact by impact.
    for r in REFERENCES {
        let landing = BallState {
            time: 0.0,
            pos: Vec3::from(r.landing_pos),
            vel: Vec3::from(r.landing_vel),
            spin: backspin(r.backspin),
        };
        let g = r.surface.settle(&landing, VACUUM);
        let label = format!(
            "{} on {}",
            r.name,
            if r.surface == Surface::GREEN {
                "GREEN"
            } else {
                "FAIRWAY"
            }
        );
        assert_eq!(g.contacts.len(), r.impacts.len(), "{label}: impact count");
        let last = r.impacts.len() - 1;
        for (i, (c, &(vel, omega, pos))) in g.contacts.iter().zip(r.impacts).enumerate() {
            let mut expected = Vec3::from(vel);
            let mut expected_spin = backspin(-omega);
            if i == last {
                // The last rebound is below the hop threshold, so the ball
                // stays on the ground, and `vel_out`/`spin_out` are what it
                // slides or rolls away with: the reference rebound projected
                // onto the ground, with a slip at or below the roll
                // threshold removed. The reference stops at the projection.
                expected.y = 0.0;
                if slip(expected, expected_spin, Vec3::Y, R).length() <= DEFAULT_ROLL_THRESHOLD {
                    (expected, expected_spin) =
                        lock_rolling(expected, expected_spin, Vec3::Y, R, BALL.inertia_factor());
                }
            }
            let scale = expected.length().max(1.0);
            assert!(
                (c.vel_out - expected).length() <= 1e-9 * scale,
                "{label} impact {i}: vel_out {:?}, reference {expected:?}",
                c.vel_out
            );
            assert!(
                c.spin_out.x == 0.0 && c.spin_out.y == 0.0,
                "{label} impact {i}: spin stays about Z: {:?}",
                c.spin_out
            );
            assert!(
                (c.spin_out.z - expected_spin.z).abs() <= 1e-9 * expected_spin.z.abs().max(1.0),
                "{label} impact {i}: ω_z {}, reference {}",
                c.spin_out.z,
                expected_spin.z
            );
            assert!(
                (c.pos - Vec3::from(pos)).length() <= 1e-8,
                "{label} impact {i}: pos {:?}, reference {pos:?}",
                c.pos
            );
        }
        // The reference integrates the slide and roll with 10 µs Euler
        // steps, this crate with 1 ms RK4 steps and an exact slide-to-roll
        // switch; the rest positions agree to within 0.01 mm.
        let rest = Vec3::from(r.rest);
        assert!(
            (g.rest.pos - rest).length() < 1e-4,
            "{label}: rest {:?}, reference {rest:?}",
            g.rest.pos
        );
    }
}

// ══════════════════════════════════════════════════════════════════
//  IMPACT — structural invariants
// ══════════════════════════════════════════════════════════════════

#[test]
fn rebound_goes_up() {
    let r = hit(Vec3::new(0.0, -10.0, 0.0), Vec3::ZERO, &Surface::FAIRWAY);
    assert!(r.vel.y > 0.0, "rebound should go up: vy={:.3}", r.vel.y);
}

#[test]
fn impact_keeps_time_and_position() {
    let s = BallState {
        time: 4.2,
        pos: Vec3::new(100.0, 0.0, -3.0),
        vel: Vec3::new(20.0, -15.0, 1.0),
        spin: Vec3::new(10.0, -60.0, 300.0),
    };
    let out = impact(&s, &Surface::GREEN.sample(0.0, 0.0), &BALL);
    assert!(out.time == s.time && out.pos == s.pos);
}

#[test]
fn no_lateral_without_sidespin() {
    let r = hit(
        Vec3::new(20.0, -15.0, 0.0),
        backspin(rpm_to_rads(3000.0)),
        &Surface::GREEN,
    );
    assert!(
        r.vel.z == 0.0 && r.spin.x == 0.0 && r.spin.y == 0.0,
        "no lateral without sidespin: vz={:.6}",
        r.vel.z
    );
}

#[test]
fn vertical_drop_no_horizontal() {
    // Ball dropping straight down with no spin should bounce straight up.
    let r = hit(Vec3::new(0.0, -5.0, 0.0), Vec3::ZERO, &Surface::GREEN);
    assert!(r.vel.y > 0.0);
    assert!(r.vel.x.abs() < 1e-12 && r.vel.z.abs() < 1e-12);
}

#[test]
fn grazing_impact_applies_no_impulse() {
    for vel in [Vec3::new(10.0, 0.0, 0.0), Vec3::new(10.0, 0.5, -2.0)] {
        let s = state(vel, Vec3::new(30.0, -80.0, 250.0));
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let out = impact(&s, &surface.sample(0.0, 0.0), &BALL);
            assert_eq!(out, s, "grazing {vel:?} should be untouched");
        }
    }
}

#[test]
fn energy_never_increases_at_impact() {
    let surfaces = [
        Surface::FAIRWAY,
        Surface::GREEN,
        Surface {
            restitution: Restitution::Constant(1.0),
            mu_slide: 0.0,
            ..Surface::GREEN
        },
        Surface {
            restitution: Restitution::Constant(1.0),
            compliance: 3.0,
            mu_slide: 1.0,
            pivot: 0.01,
            ..Surface::GREEN
        },
    ];
    let vels = [
        Vec3::new(20.0, -15.0, 0.0),
        Vec3::new(25.0, -10.0, 4.0),
        Vec3::new(-5.0, -3.0, -2.0),
        Vec3::new(3.0, -20.0, 0.5),
        Vec3::new(0.0, -8.0, 0.0),
    ];
    let spins = [
        Vec3::ZERO,
        backspin(800.0),
        Vec3::new(0.0, -100.0, 300.0),
        Vec3::new(200.0, 50.0, -100.0),
        Vec3::new(-150.0, 0.0, 0.0),
    ];
    for surface in &surfaces {
        for &v in &vels {
            for &w in &spins {
                let out = hit(v, w, surface);
                let e_in = energy(v, w, 0.0);
                let e_out = energy(out.vel, out.spin, 0.0);
                assert!(
                    e_out <= e_in * (1.0 + 1e-12),
                    "energy rose {e_in:.6} → {e_out:.6} for v={v:?} ω={w:?} on {surface:?}"
                );
            }
        }
    }
}

#[test]
fn impact_mirror_symmetry() {
    // Mirroring the impact in z (a pseudovector spin flips x and y) mirrors
    // the outcome; mirroring in x (spin flips y and z) likewise.
    let sideways = |v: Vec3| Vec3::new(v.x, v.y, -v.z);
    let sideways_spin = |w: Vec3| Vec3::new(-w.x, -w.y, w.z);
    let backward = |v: Vec3| Vec3::new(-v.x, v.y, v.z);
    let backward_spin = |w: Vec3| Vec3::new(w.x, -w.y, -w.z);
    for (vel, spin) in [
        (Vec3::new(5.0, -3.0, 0.0), backspin(-5.0 / R)),
        (Vec3::new(15.0, -10.0, 2.0), Vec3::new(20.0, -60.0, 300.0)),
        (Vec3::new(23.8, -16.0, -1.0), Vec3::new(0.0, 80.0, 900.0)),
    ] {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let fwd = hit(vel, spin, &surface);
            let z = hit(sideways(vel), sideways_spin(spin), &surface);
            assert!((z.vel - sideways(fwd.vel)).length() < 1e-12);
            assert!((z.spin - sideways_spin(fwd.spin)).length() < 1e-9);
            let x = hit(backward(vel), backward_spin(spin), &surface);
            assert!((x.vel - backward(fwd.vel)).length() < 1e-12);
            assert!((x.spin - backward_spin(fwd.spin)).length() < 1e-9);
        }
    }
}

#[test]
fn impact_frame_invariance_about_y() {
    let a = 0.7_f64;
    let rot = |v: Vec3| {
        Vec3::new(
            v.x * a.cos() - v.z * a.sin(),
            v.y,
            v.x * a.sin() + v.z * a.cos(),
        )
    };
    for (vel, spin) in [
        (Vec3::new(20.0, -15.0, 1.5), Vec3::new(10.0, -60.0, 250.0)),
        (Vec3::new(6.0, -2.0, -0.5), Vec3::new(-30.0, 40.0, -200.0)),
    ] {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let base = hit(vel, spin, &surface);
            let turned = hit(rot(vel), rot(spin), &surface);
            assert!((turned.vel - rot(base.vel)).length() < 1e-12);
            assert!((turned.spin - rot(base.spin)).length() < 1e-9);
        }
    }
}

// ══════════════════════════════════════════════════════════════════
//  RESTITUTION — Penner (2002) eq. 5
// ══════════════════════════════════════════════════════════════════

#[test]
fn penner2002_restitution_matches_eq5() {
    // Penner (2002) eq. 5: e = 0.510 − 0.0375·v + 0.000903·v² for
    // v ≤ 20 m/s, 0.120 above.
    let e = Restitution::Penner2002;
    for (v, expected) in [
        (0.0, 0.510),
        (5.0, 0.345_075),
        (10.0, 0.225_3),
        (15.0, 0.150_675),
        (20.0, 0.121_2),
    ] {
        assert!(
            (e.coefficient_at(v) - expected).abs() < 1e-12,
            "e({v}) = {:.6}, expected {expected}",
            e.coefficient_at(v)
        );
    }
    assert!((e.coefficient_at(-10.0) - e.coefficient_at(10.0)).abs() < f64::EPSILON);
}

#[test]
fn penner2002_restitution_constant_above_20_mps() {
    let e = Restitution::Penner2002;
    for v in [20.001, 25.0, 40.0] {
        assert!(
            (e.coefficient_at(v) - 0.120).abs() < f64::EPSILON,
            "e({v}) = {}",
            e.coefficient_at(v)
        );
    }
    // The two branches meet within 0.0012 at 20 m/s.
    assert!((e.coefficient_at(20.0) - e.coefficient_at(20.001)).abs() < 0.002);
}

#[test]
fn constant_restitution_ignores_speed() {
    let e = Restitution::Constant(0.3);
    assert!((e.coefficient_at(1.0) - 0.3).abs() < f64::EPSILON);
    assert!((e.coefficient_at(30.0) - 0.3).abs() < f64::EPSILON);
}

// ══════════════════════════════════════════════════════════════════
//  BACKSPIN AT IMPACT
// ══════════════════════════════════════════════════════════════════

#[test]
fn fairway_bounces_lower_than_green_and_resists_rolling_more() {
    // The same landing leaves FAIRWAY lower than GREEN (32.5° against
    // 35.0° at 15 m/s, 34° and 3000 rpm): the measured greens are softer
    // and grip a little harder. On the ground the fairway resists rolling
    // more, stimp 5 ft 8 in against 10.5 ft.
    let vel = Vec3::new(15.0, -10.0, 0.0);
    let spin = backspin(rpm_to_rads(3000.0));
    let angle = |v: Vec3| v.y.atan2(v.x);
    let fairway = hit(vel, spin, &Surface::FAIRWAY);
    let green = hit(vel, spin, &Surface::GREEN);
    assert!(
        angle(fairway.vel) < angle(green.vel),
        "{:?} vs {:?}",
        fairway.vel,
        green.vel
    );
    assert!(Surface::FAIRWAY.stimp_feet() < 0.6 * Surface::GREEN.stimp_feet());
}

#[test]
fn backspin_reduces_forward_speed() {
    // Up to the friction bound, more backspin leaves less forward speed: on
    // GREEN a 20 m/s, 53° landing leaves at 3.35 m/s with 2000 rpm and at
    // 2.22 m/s with 8000 rpm, where the impact slides throughout.
    let vel = Vec3::new(12.0, -16.0, 0.0);
    let low = hit(vel, backspin(rpm_to_rads(2000.0)), &Surface::GREEN);
    let high = hit(vel, backspin(rpm_to_rads(8000.0)), &Surface::GREEN);
    assert!(
        high.vel.x < low.vel.x,
        "high spin should reduce forward speed: high={:.3}, low={:.3}",
        high.vel.x,
        low.vel.x
    );
}

#[test]
fn very_high_backspin_can_reverse() {
    // Wedge at 20 m/s, 55° with 10000 rpm on GREEN: the first impact slides
    // throughout and the ball hops forward still carrying backspin, then
    // spins back and comes to rest behind its pitch mark.
    let angle = units::deg_to_rad(55.0);
    let vel = 20.0 * Vec3::new(angle.cos(), -angle.sin(), 0.0);
    let spin = backspin(rpm_to_rads(10000.0));
    let first = hit(vel, spin, &Surface::GREEN);
    assert!(first.vel.x > 0.0 && first.spin.z > 0.0, "{first:?}");
    let g = Surface::GREEN.settle(&state(vel, spin), air());
    assert!(g.rest.pos.x < 0.0, "rests at {:.2} m", g.rest.pos.x);
}

#[test]
fn driver_bounces_forward() {
    let angle = units::deg_to_rad(30.0);
    let vel = 30.0 * Vec3::new(angle.cos(), -angle.sin(), 0.0);
    let r = hit(vel, backspin(rpm_to_rads(2500.0)), &Surface::FAIRWAY);
    assert!(
        r.vel.x > 0.0,
        "driver should bounce forward: vx={:.3}",
        r.vel.x
    );
}

#[test]
fn seven_iron_checks_on_green() {
    let angle = units::deg_to_rad(45.0);
    let vel = 25.0 * Vec3::new(angle.cos(), -angle.sin(), 0.0);
    let spin = hit(vel, backspin(rpm_to_rads(7000.0)), &Surface::GREEN);
    let no_spin = hit(vel, Vec3::ZERO, &Surface::GREEN);
    assert!(spin.vel.x < no_spin.vel.x);
}

// ══════════════════════════════════════════════════════════════════
//  PENNER COMPLIANCE — the tilted effective contact plane
// ══════════════════════════════════════════════════════════════════

fn plane(restitution: f64, mu: f64, compliance: f64) -> Surface {
    Surface {
        restitution: Restitution::Constant(restitution),
        mu_slide: mu,
        compliance,
        ..Surface::GREEN
    }
}

#[test]
fn compliance_steepens_rebound() {
    let vel = Vec3::new(20.0, -15.0, 0.0);
    let spin = backspin(rpm_to_rads(5000.0));
    let firm = hit(vel, spin, &plane(0.35, 0.4, 0.2));
    let soft = hit(vel, spin, &plane(0.35, 0.4, 1.0));
    let angle = |v: Vec3| v.y.atan2(v.x.abs().max(0.001));
    assert!(angle(soft.vel) > angle(firm.vel));
}

#[test]
fn compliance_rotates_mirror_rebound_by_two_beta() {
    // Elastic, frictionless impact at 45°: the effective plane is tilted by
    // the fitted contact angle toward the ball, so the specular rebound is the flat-ground mirror direction rotated up by
    // 2β, at full speed. The same holds for any heading.
    for heading in [0.0_f64, 0.4, -2.0] {
        let h = Vec3::new(heading.cos(), 0.0, heading.sin());
        let vel = 10.0 * h - Vec3::new(0.0, 10.0, 0.0);
        let r = hit(vel, Vec3::ZERO, &plane(1.0, 0.0, 1.0));
        let incidence = std::f64::consts::FRAC_PI_4;
        let beta = fitted_tilt(10.0, 10.0);
        let rebound_angle = r.vel.y.atan2(r.vel.dot(h));
        assert!((rebound_angle - (incidence + 2.0 * beta)).abs() < 1e-12);
        assert!((r.vel.length() - vel.length()).abs() < 1e-12);
        assert!(
            r.vel.dot(Vec3::Y.cross(h)).abs() < 1e-12,
            "stays in the plane of travel"
        );
    }
}

#[test]
fn no_compliance_is_exact_mirror() {
    let r = hit(
        Vec3::new(10.0, -10.0, 0.0),
        Vec3::ZERO,
        &plane(1.0, 0.0, 0.0),
    );
    assert!((r.vel - Vec3::new(10.0, 10.0, 0.0)).length() < 1e-12);
}

#[test]
fn compliance_tilt_capped_at_incidence_angle() {
    // 10 × 11.3° = 113° would overshoot; β is capped at the 45° incidence
    // angle, so the ball meets the effective plane head-on and an elastic
    // frictionless impact sends it straight back.
    let vel = Vec3::new(10.0, -10.0, 0.0);
    let r = hit(vel, Vec3::ZERO, &plane(1.0, 0.0, 10.0));
    assert!((r.vel + vel).length() < 1e-9, "got {:?}", r.vel);
}

/// The contact-angle law at compliance 1.0, written out independently of
/// `Surface::contact_angle`: tan β = (`v_t` / 50 m/s)·(1 − exp(−(θ / 20°)³)).
fn fitted_tilt(v_t: f64, v_n: f64) -> f64 {
    let theta = v_n.atan2(v_t);
    let factor = 1.0 - (-(theta / units::deg_to_rad(20.0)).powi(3)).exp();
    (v_t / 50.0 * factor).atan()
}

#[test]
fn contact_angle_follows_the_fitted_law() {
    // The four landings Penner (2002) publishes, a green impact, two drives
    // and a 9-iron (speed, degrees from the vertical), at compliance 1.0:
    // 14.6°, 20.1°, 20.7° and 14.9°. The presets multiply the angle by their compliance.
    let firm = Surface {
        compliance: 1.0,
        ..Surface::FAIRWAY
    };
    for (speed, deg, tilt) in [
        (18.6, 44.4, 14.589),
        (26.8, 68.0, 20.086),
        (26.2, 46.2, 20.716),
        (23.8, 34.0, 14.905),
    ] {
        let phi = units::deg_to_rad(deg);
        let (v_t, v_n) = (speed * phi.sin(), speed * phi.cos());
        let beta = firm.contact_angle(v_t, v_n);
        assert!((beta - fitted_tilt(v_t, v_n)).abs() < 1e-12);
        assert!(
            (units::rad_to_deg(beta) - tilt).abs() < 1e-3,
            "{speed} m/s at {deg}°: β {}",
            units::rad_to_deg(beta)
        );
        for surface in [Surface::FAIRWAY, Surface::GREEN, Surface::FIRM_GREEN] {
            let scaled = surface.contact_angle(v_t, v_n);
            assert!((scaled - surface.compliance * beta).abs() < 1e-12);
        }
    }
    // A steep landing tilts by atan(v_t / 50 m/s) whatever its normal speed:
    // at 45° and steeper the factor in θ is within 2e-5 of 1.
    for (v_t, v_n) in [(10.0, 10.0), (10.0, 30.0), (30.0, 40.0)] {
        let beta = firm.contact_angle(v_t, v_n);
        assert!(
            ((beta.tan() * 50.0 / v_t) - 1.0).abs() < 2e-5,
            "{v_t}, {v_n}"
        );
    }
    // A thin shot landing at 35 m/s and 10°: the factor in θ is
    // 1 − exp(−1/8) = 0.1175, so β = 4.63°, where the steep law alone would
    // give 34.6°.
    let theta = units::deg_to_rad(10.0);
    let (v_t, v_n) = (35.0 * theta.cos(), 35.0 * theta.sin());
    let beta = firm.contact_angle(v_t, v_n);
    assert!((units::rad_to_deg(beta) - 4.631).abs() < 1e-3, "β {beta}");
    // `compliance` scales the angle, not its tangent.
    let half = Surface {
        compliance: 0.5,
        ..firm
    };
    assert!((half.contact_angle(v_t, v_n) - 0.5 * beta).abs() < 1e-12);
    // The tilt vanishes as the normal speed goes to zero.
    assert!(firm.contact_angle(30.0, 0.01) < 1e-6);
    assert_eq!(firm.contact_angle(0.0, 10.0), 0.0);
}

#[test]
fn fast_low_landing_hops_forward() {
    // A landing at 46 m/s and 15° (normal speed 11.9 m/s) never turns
    // backward and its first hop stays under 10 m (about 5 m without spin).
    for surface in [Surface::FAIRWAY, Surface::GREEN, Surface::FIRM_GREEN] {
        for rpm in [0.0, 2000.0, 4000.0] {
            let theta = units::deg_to_rad(15.0);
            let vel = Vec3::new(46.0 * theta.cos(), -46.0 * theta.sin(), 0.0);
            let g = surface.settle(&state(vel, backspin(rpm_to_rads(rpm))), air());
            let label = format!("{rpm} rpm");
            for p in &g.points {
                assert!(p.vel.x >= 0.0, "{label}: moves backward {:?}", p.vel);
            }
            let end = g.contacts.get(1).map_or(f64::INFINITY, |c| c.time);
            let hop = g
                .points
                .iter()
                .filter(|p| p.phase == Phase::Hop && p.time <= end)
                .map(|p| p.pos.y)
                .fold(0.0, f64::max);
            assert!(hop < 10.0, "{label}: first hop {hop:.2} m");
        }
    }
}

#[test]
fn shallow_landings_skip_forward_and_lower() {
    // Landings at 2–10° from 20, 35 and 46 m/s with 0, 2000 and 4000 rpm of
    // backspin on every preset: the first impact keeps at least 80% of the
    // forward speed, no rebound leaves faster along the normal than the
    // ball arrived, and each hop is no higher than the one before.
    for surface in [Surface::FAIRWAY, Surface::GREEN, Surface::FIRM_GREEN] {
        for deg in [2.0, 4.0, 6.0, 8.0, 10.0] {
            for speed in [20.0, 35.0, 46.0] {
                for rpm in [0.0, 2000.0, 4000.0] {
                    let theta = units::deg_to_rad(deg);
                    let vel = Vec3::new(speed * theta.cos(), -speed * theta.sin(), 0.0);
                    let g = surface.settle(&state(vel, backspin(rpm_to_rads(rpm))), air());
                    let label = format!("{speed} m/s at {deg}°, {rpm} rpm");
                    let first = &g.contacts[0];
                    assert!(
                        first.vel_out.x >= 0.8 * vel.x,
                        "{label}: keeps {:.2} of {:.2} m/s",
                        first.vel_out.x,
                        vel.x
                    );
                    let mut last_up = f64::INFINITY;
                    for c in &g.contacts {
                        assert!(
                            c.vel_out.y <= -c.vel_in.y,
                            "{label}: rebound {:.3} m/s from {:.3} m/s",
                            c.vel_out.y,
                            -c.vel_in.y
                        );
                        assert!(
                            c.vel_out.y <= last_up,
                            "{label}: hop at {:.3} m/s after {last_up:.3} m/s",
                            c.vel_out.y
                        );
                        last_up = c.vel_out.y;
                    }
                    assert_eq!(g.termination, Termination::Rest, "{label}");
                    assert!(g.rest.pos.x > 0.0, "{label}: rests behind the mark");
                }
            }
        }
    }
}

#[test]
fn low_skimming_hops_lose_height() {
    // A ball skimming along at 2–4 m/s in low hops (normal speed 0.3–1.5
    // m/s, rising 5–115 mm) gets little tilt from its flat contacts, where
    // the landing-angle factor of the contact angle is small, so each
    // rebound leaves lower than the ball arrived and the hops die out
    // within a few impacts.
    for surface in [Surface::FAIRWAY, Surface::GREEN] {
        for v_t in [2.0, 3.0, 4.0] {
            for w in [0.3, 0.5, 1.0, 1.5] {
                let r = hit(Vec3::new(v_t, -w, 0.0), rolling(v_t), &surface);
                assert!(
                    r.vel.y < 0.85 * w,
                    "{v_t} m/s skimming at {w} m/s rebounds at {:.3} m/s",
                    r.vel.y
                );
            }
        }
    }
}

#[test]
fn spin_back_into_ground_rolls_instead_of_rebounding_down() {
    // Dead turf (e=0.05), strong tilt and extreme backspin: the ball grips
    // and leaves the effective plane moving backward fast enough that the
    // rebound points into the ground. It must stay on the ground and roll
    // back, never rebound downward.
    let surface = Surface {
        restitution: Restitution::Constant(0.05),
        mu_slide: 1.0,
        compliance: 2.0,
        ..Surface::GREEN
    };
    let vel = Vec3::new(10.0, -10.0, 0.0);
    let raw = hit(vel, backspin(1500.0), &surface);
    assert!(raw.vel.y < 0.0, "impulse alone points into the ground");
    let g = surface.settle(&state(vel, backspin(1500.0)), air());
    assert_eq!(g.contacts.len(), 1);
    assert!(g.contacts[0].vel_out.y == 0.0);
    assert!(g.contacts[0].vel_out.x < 0.0);
    assert!(g.rest.pos.x < 0.0, "rolls back: x={:.3}", g.rest.pos.x);
}

#[test]
fn ball_rolling_backward_keeps_rolling_backward() {
    // A ball moving back toward the tee with matching rolling spin grips
    // and leaves still moving backward, spinning in the backward-rolling
    // sense (+Z for travel along −X).
    let r = hit(
        Vec3::new(-5.0, -3.0, 0.0),
        backspin(5.0 / R),
        &Surface::GREEN,
    );
    assert!(
        r.vel.x < -4.0 && r.vel.y > 0.0 && r.spin.z > 150.0,
        "leaves rolling backward: v={:?}, ω={:?}",
        r.vel,
        r.spin
    );
}

// ══════════════════════════════════════════════════════════════════
//  BOUNCE SEQUENCE
// ══════════════════════════════════════════════════════════════════

#[test]
fn bounce_sequence_terminates() {
    let g = Surface::FAIRWAY.settle(
        &state(Vec3::new(20.0, -15.0, 0.0), backspin(rpm_to_rads(3000.0))),
        air(),
    );
    assert_ends_by_hop_height(&g);
}

#[test]
fn bounce_sequence_stops_below_threshold() {
    // Straight drop at 15 m/s: no tilt, so each rebound is e(v)·v with
    // Penner's e rising toward 0.51 as the ball slows. Rebounds of 2.26,
    // 0.97, 0.46 and 0.23 m/s rise 260, 48, 10.8 and 2.6 mm: the fourth is
    // below the 5 mm transition height.
    let g = Surface::FAIRWAY.settle(&state(Vec3::new(0.0, -15.0, 0.0), Vec3::ZERO), air());
    assert_eq!(g.contacts.len(), 4);
    assert_ends_by_hop_height(&g);
    assert!(g.rest.pos.length() < 1e-9);
}

#[test]
fn bounce_sequence_safety_bound() {
    // Elastic, frictionless, rigid, no air: every rebound returns at full
    // speed, so only the impact limit ends the sequence (each ~1 s hop
    // would hit the default 120 s time limit first, so it is raised here).
    let physics = GroundPhysics {
        max_time: 1000.0,
        ..GroundPhysics::new(plane(1.0, 0.0, 0.0))
    };
    let start = state(Vec3::new(5.0, -5.0, 0.0), Vec3::ZERO);
    let g = physics.settle(&start, VACUUM);
    assert_eq!(g.termination, Termination::MaxImpacts);
    assert_eq!(g.contacts.len(), 200);
    assert!(g.contacts.iter().all(|c| c.vel_out.y > 4.99));
    // The path ends with a point at the last contact, like every other exit.
    let end = g.points.last().unwrap();
    assert_eq!(*end, GroundPoint::new(&g.rest, Phase::Hop));
    assert_eq!(end.time, g.contacts[199].time);
    // With the default limits the same ball runs out of time instead.
    let timed = plane(1.0, 0.0, 0.0).settle(&start, VACUUM);
    assert_eq!(timed.termination, Termination::MaxTime);
    assert!((timed.rest.time - 120.0).abs() < 1e-9);
}

#[test]
fn bounce_sequence_height_decreases() {
    let g = Surface::FAIRWAY.settle(
        &state(Vec3::new(15.0, -10.0, 0.0), backspin(rpm_to_rads(2000.0))),
        air(),
    );
    for w in g.contacts.windows(2) {
        assert!(w[1].vel_out.y <= w[0].vel_out.y + 0.01);
    }
}

#[test]
fn contacts_record_true_arrival_velocity() {
    // Each contact's vel_in is the velocity arriving at the ground: on flat
    // ground with ballistic hops, the previous rebound with its vertical
    // component reversed.
    let g = Surface::GREEN.settle(
        &state(Vec3::new(12.0, -9.0, 1.0), Vec3::new(0.0, -40.0, 400.0)),
        VACUUM,
    );
    assert!(g.contacts.len() > 3);
    for w in g.contacts.windows(2) {
        let expected = Vec3::new(w[0].vel_out.x, -w[0].vel_out.y, w[0].vel_out.z);
        assert!((w[1].vel_in - expected).length() < 1e-9);
        assert!(w[1].vel_in.y < 0.0);
        assert_eq!(w[1].spin_in, w[0].spin_out);
        assert!(w[1].time > w[0].time);
    }
}

#[test]
fn soft_drop_ends_by_hop_height() {
    let g = Surface::GREEN.settle(
        &state(Vec3::new(5.0, -1.0, 0.0), backspin(rpm_to_rads(500.0))),
        air(),
    );
    assert_ends_by_hop_height(&g);
}

// ══════════════════════════════════════════════════════════════════
//  HOPS THROUGH AIR — drag and Magnus between impacts
// ══════════════════════════════════════════════════════════════════

/// A single hop on `FAIRWAY`: a ball leaving the ground at `vel` (moving
/// up, so the first contact applies no impulse) with `spin`, through the
/// air of `env`. Returns the arrival contact.
fn hop_with(vel: Vec3, spin: Vec3, env: &Tour) -> Contact {
    let g = Surface::FAIRWAY.settle(&state(vel, spin), env);
    assert_eq!(g.contacts[0].vel_out, vel, "no impulse leaving the ground");
    g.contacts[1]
}

#[test]
fn vacuum_hops_are_bit_identical() {
    // At zero density `Surface` and `GroundPhysics` give the same result to
    // the bit, and every arrival is on the parabola from the rebound before
    // it, with the spin carried unchanged.
    let start = state(Vec3::new(21.3, -18.5, 0.0), backspin(212.0));
    let physics = GroundPhysics::new(Surface::FAIRWAY).settle(&start, VACUUM);
    let surface = Surface::FAIRWAY.settle(&start, VACUUM);
    assert_eq!(physics, surface);
    assert!(surface.contacts.len() > 5);
    for w in surface.contacts.windows(2) {
        let tau = w[1].time - w[0].time;
        let on_parabola = w[0].pos + w[0].vel_out * tau + 0.5 * GRAVITY * tau * tau;
        assert!((on_parabola - w[1].pos).length() < 1e-9);
        assert!((w[1].vel_in - (w[0].vel_out + GRAVITY * tau)).length() < 1e-9);
        assert_eq!(w[1].spin_in, w[0].spin_out);
    }
    assert_ne!(surface, Surface::FAIRWAY.settle(&start, air()));
}

#[test]
fn hops_lose_speed_to_drag() {
    // The same landing with and without air: the first impact is the same,
    // and the hop through air arrives earlier, nearer and slower, on the
    // ground, still carrying the rebound's spin exactly.
    let start = state(Vec3::new(21.3, -18.5, 0.0), backspin(212.0));
    let vacuum = Surface::FAIRWAY.settle(&start, VACUUM);
    let air = Surface::FAIRWAY.settle(&start, air());
    assert_eq!(vacuum.contacts[0], air.contacts[0]);
    let (v, a) = (&vacuum.contacts[1], &air.contacts[1]);
    assert!(a.time < v.time, "{} vs {} s", a.time, v.time);
    assert!(a.pos.x < v.pos.x, "{} vs {} m", a.pos.x, v.pos.x);
    assert!(a.vel_in.length() < v.vel_in.length());
    for g in [&vacuum, &air] {
        for c in &g.contacts {
            assert_eq!(c.pos.y, 0.0);
        }
        for w in g.contacts.windows(2) {
            assert_eq!(w[1].spin_in, w[0].spin_out);
        }
    }
    assert!(air.rest.pos.x < vacuum.rest.pos.x);
}

#[test]
fn hop_drag_scales_with_density() {
    // Thinner air takes less of the run: the loss against ballistic hops
    // shrinks from sea level to Denver to Mexico City.
    let start = state(Vec3::new(21.3, -18.5, 0.0), backspin(212.0));
    let run = |env: &Tour| Surface::FAIRWAY.settle(&start, env).rest.pos.x;
    let ballistic = run(VACUUM);
    let loss = |env: Tour| ballistic - run(&env);
    let at = |altitude_m, temperature_c| {
        Tour::new(
            crate::air::standard_air_density(altitude_m, temperature_c),
            Wind::None,
        )
    };
    let (sea, denver, mexico) = (
        loss(Tour::SEA_LEVEL),
        loss(at(1609.0, 17.0)),
        loss(at(2240.0, 16.0)),
    );
    assert!(
        sea > denver && denver > mexico && mexico > 0.0,
        "{sea:.3}, {denver:.3}, {mexico:.3} m"
    );
}

#[test]
fn topspin_hop_lands_sooner_than_spinless() {
    // Magnus lift is spin × velocity: topspin (forward rolling spin) pushes
    // a hop down, backspin holds it up. The same spin rate either way, so
    // the drag is the same.
    let vel = Vec3::new(6.0, 2.0, 0.0);
    let w = 6.0 / R;
    let top = hop_with(vel, backspin(-w), air());
    let none = hop_with(vel, Vec3::ZERO, air());
    let back = hop_with(vel, backspin(w), air());
    assert!(
        top.time < none.time && none.time < back.time,
        "{:.6} / {:.6} / {:.6} s",
        top.time,
        none.time,
        back.time
    );
    assert!(top.pos.x < back.pos.x);
    // The spin rides through the hop unchanged.
    assert_eq!(top.spin_in, backspin(-w));
    // Without air the spin has no effect on the hop.
    let vacuum = |spin| hop_with(vel, spin, VACUUM);
    assert_eq!(vacuum(backspin(-w)).time, vacuum(backspin(w)).time);
}

// ══════════════════════════════════════════════════════════════════
//  SPIN VECTOR ON THE GROUND — kicks and drift
// ══════════════════════════════════════════════════════════════════

#[test]
fn kick_matches_grip_formula() {
    // Rigid ground (no tilt), grip regime: the ball leaves rolling, and its
    // heading becomes atan((5/7)v_z / ((5/7)v_x − (2/7)Rω_b)) for a ball
    // arriving with velocity (v_x, ·, v_z) and backspin ω_b about +Z.
    let surface = plane(0.5, 1.0, 0.0);
    for (vel, w) in [
        (Vec3::new(20.0, -12.0, 3.0), 250.0),
        (Vec3::new(15.0, -14.0, -2.0), 300.0),
        (Vec3::new(30.0, -10.0, 4.0), 150.0),
    ] {
        let out = hit(vel, backspin(w), &surface);
        let expected = ((5.0 / 7.0) * vel.z).atan2((5.0 / 7.0) * vel.x - (2.0 / 7.0) * R * w);
        let heading = out.vel.z.atan2(out.vel.x);
        assert!(
            (heading - expected).abs() < 1e-12,
            "heading {:.6}°, formula {:.6}°",
            units::rad_to_deg(heading),
            units::rad_to_deg(expected)
        );
        // The kick enlarges the heading angle (fade right, draw left).
        let heading_in = vel.z.atan2(vel.x);
        assert!(heading.abs() > heading_in.abs() && heading.signum() == heading_in.signum());
        // Grip: the contact point leaves without slip.
        assert!(slip(out.vel, out.spin, Vec3::Y, R).length() < 1e-12);
    }
}

#[test]
fn spin_about_normal_gives_no_kick_on_rigid_ground() {
    // Pure ω_y cannot move the contact point of a ball on flat rigid
    // ground: no lateral velocity at impact, and pivot friction only
    // decays it.
    let surface = Surface {
        compliance: 0.0,
        ..Surface::FAIRWAY
    };
    for w_y in [-150.0, 150.0] {
        let out = hit(
            Vec3::new(15.0, -10.0, 0.0),
            Vec3::new(0.0, w_y, 0.0),
            &surface,
        );
        assert!(out.vel.z == 0.0, "vz={}", out.vel.z);
        assert!(out.spin.y.abs() < w_y.abs() && out.spin.y.signum() == w_y.signum());
    }
}

#[test]
fn spin_about_normal_gives_no_drift_on_the_ground() {
    // Sliding (with backspin) and rolling balls carrying spin about the
    // vertical go straight; pivot friction removes the vertical spin.
    for spin in [
        Vec3::new(0.0, 120.0, -4.0 / R),
        Vec3::new(0.0, -120.0, 400.0),
    ] {
        let g = roll(Vec3::new(4.0, 0.0, 0.0), spin, &Surface::FAIRWAY);
        assert!(g.points.iter().all(|p| p.pos.z == 0.0));
        assert!(g.rest.pos.z == 0.0);
        assert!(g.rest.spin.y.abs() < spin.y.abs());
    }
}

#[test]
fn pivot_friction_at_impact() {
    // Δω_n = pivot · j_n / (k R²) about the surface normal. Rigid ground
    // keeps n' = ŷ, so j_n and the expected change are closed-form.
    let surface = Surface {
        restitution: Restitution::Constant(0.4),
        compliance: 0.0,
        ..Surface::FAIRWAY
    };
    let vel = Vec3::new(15.0, -10.0, 0.0);
    let out = hit(vel, Vec3::new(0.0, -200.0, 0.0), &surface);
    let j_n = 1.4 * 10.0;
    let expected = -200.0 + surface.pivot * j_n / (0.4 * R * R);
    assert!(
        (out.spin.y - expected).abs() < 1e-9,
        "ω_y {} vs {expected}",
        out.spin.y
    );
    // Small spins stop rather than reverse.
    let out = hit(vel, Vec3::new(0.0, 5.0, 0.0), &surface);
    assert!(out.spin.y == 0.0);

    // On tilted (compliant) impacts it still acts about the true normal ŷ:
    // the velocity is untouched and only ω_y changes, so the spin along the
    // heading that drives a curving ball's kick is kept.
    for spin in [
        Vec3::new(30.0, -60.0, 300.0),
        Vec3::new(-80.0, 400.0, 900.0),
    ] {
        let vel = Vec3::new(20.0, -14.0, 3.0);
        let free = hit(
            vel,
            spin,
            &Surface {
                pivot: 0.0,
                ..Surface::GREEN
            },
        );
        let with = hit(vel, spin, &Surface::GREEN);
        assert_eq!(with.vel, free.vel);
        assert!(with.spin.x == free.spin.x && with.spin.z == free.spin.z);
        assert!(with.spin.y.abs() < free.spin.y.abs());
    }
}

#[test]
fn tilted_impact_throws_against_sidespin() {
    // On compliant turf the impact is against Penner's tilted plane, so
    // part of a spin about the vertical lies in the contact plane and slips
    // the contact point sideways. A ball spinning clockwise from above
    // (fade sidespin, ω_y < 0) is thrown left, like side spin off a cushion.
    let vel = Vec3::new(18.0, -15.0, 0.0);
    let fade = hit(vel, Vec3::new(0.0, -60.0, 0.0), &Surface::GREEN);
    let draw = hit(vel, Vec3::new(0.0, 60.0, 0.0), &Surface::GREEN);
    assert!(fade.vel.z < 0.0 && draw.vel.z > 0.0);
    assert!((fade.vel.z + draw.vel.z).abs() < 1e-12);
}

#[test]
fn tilt_throw_magnitude() {
    // A ball moving along +X with pure spin ω_y = w meets the effective
    // plane n' = (−sin β, cos β, 0). Its contact point slips sideways at
    // u_z = −R·w·sin β, alongside the in-plane slip v − (v·n')n'. Friction
    // along −û gives the lateral impulse
    //   grip:  Δv_z = (k/(1+k))·R·w·sin β
    //   slip:  Δv_z = μ·j_n·R·w·sin β / |u|,  j_n = (1+e)·|v·n'|
    let k = BALL.inertia_factor();
    let e = 0.4;
    let vel = Vec3::new(18.0, -15.0, 0.0);
    let beta = fitted_tilt(vel.x, -vel.y);
    let n_eff = Vec3::new(-beta.sin(), beta.cos(), 0.0);
    let v_n = vel.dot(n_eff);
    for w in [-60.0, 60.0, -150.0] {
        let lateral_slip = -R * w * beta.sin();
        let grip = hit(vel, Vec3::new(0.0, w, 0.0), &plane(e, 10.0, 1.0));
        let expected = -k / (1.0 + k) * lateral_slip;
        assert!(
            (grip.vel.z - expected).abs() < 1e-12,
            "grip, ω_y {w}: v_z {:.9}, expected {expected:.9}",
            grip.vel.z
        );

        let mu = 0.05;
        let sliding = hit(vel, Vec3::new(0.0, w, 0.0), &plane(e, mu, 1.0));
        let in_plane = (vel - v_n * n_eff).length();
        let slip_speed = in_plane.hypot(lateral_slip);
        assert!(mu * (1.0 + e) * -v_n < slip_speed * k / (1.0 + k), "slides");
        let expected = -mu * (1.0 + e) * -v_n * lateral_slip / slip_speed;
        assert!(
            (sliding.vel.z - expected).abs() < 1e-12,
            "slip, ω_y {w}: v_z {:.9}, expected {expected:.9}",
            sliding.vel.z
        );
        // Fade sidespin (ω_y < 0) is thrown left.
        assert!(grip.vel.z * w > 0.0 && sliding.vel.z * w > 0.0);
    }
}

#[test]
fn rifle_spin_drifts_laterally_while_sliding() {
    // Spin about the direction of travel slips the contact point sideways:
    // ω_x > 0 moves the bottom of the ball toward −Z, so friction pushes the
    // ball toward +Z while it slides, then it rolls on its new heading.
    for w_x in [200.0, -200.0] {
        let g = roll(
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::new(w_x, 0.0, -4.0 / R),
            &Surface::GREEN,
        );
        let slide_end = g
            .points
            .iter()
            .find(|p| p.phase == Phase::Roll)
            .expect("ball ends up rolling");
        assert!(slide_end.vel.z * w_x > 0.0, "drifting toward ω_x sign");
        assert!(g.rest.pos.z * w_x > 0.0, "rest z {:.3}", g.rest.pos.z);
    }
}

// ══════════════════════════════════════════════════════════════════
//  GROUND PHASE — slide, roll, rest
// ══════════════════════════════════════════════════════════════════

#[test]
fn slide_ends_exactly_at_roll_time() {
    // On flat ground the slip falls linearly at ((1+k)/k)·μ·g; the slide is
    // stepped exactly to its end, where the ball leaves without slip.
    let vel = Vec3::new(3.0, 0.0, 0.5);
    let spin = Vec3::new(-20.0, 30.0, 300.0);
    let u0 = slip(vel, spin, Vec3::Y, R).length();
    let expected = u0 / (3.5 * 0.40 * G);
    let g = roll(vel, spin, &TEST_TURF);
    let roll_start = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
    assert!(
        (roll_start.time - expected).abs() < 1e-12,
        "roll starts at {:.9} s, expected {expected:.9} s",
        roll_start.time
    );
    assert!(slip(roll_start.vel, roll_start.spin, Vec3::Y, R).length() < 1e-12);
    // Sliding friction acts along the initial slip direction throughout:
    // velocity change is exactly −μ g t û.
    let u_hat = slip(vel, spin, Vec3::Y, R).normalize();
    let dv = roll_start.vel - vel;
    let expected_dv = -0.40 * G * expected * u_hat;
    assert!(
        (dv - expected_dv).length() < 1e-9,
        "Δv {dv:?} vs {expected_dv:?}"
    );
}

#[test]
fn stimpmeter_roll_matches_closed_form() {
    // A stimpmeter roll leaves the ramp at 1.83 m/s, already rolling. With
    // a constant rolling resistance a = μ_roll·g/(1+k) the distance is
    // (1+k)·v₀²/(2·μ_roll·g).
    let v0 = STIMPMETER_SPEED;
    let g = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &TEST_TURF);
    let expected = 1.4 * v0 * v0 / (2.0 * TEST_TURF.mu_roll * G);
    assert!(
        (g.rest.pos.x - expected).abs() < 1e-5,
        "roll {:.6} m, closed form {expected:.6} m",
        g.rest.pos.x
    );
    assert!((units::meters_to_feet(g.rest.pos.x) - TEST_TURF.stimp_feet()).abs() < 1e-4);
    // With plowing c as well, a = (μ_roll·g + c·v)/(1+k) and the distance
    // is (1+k)·[v₀/c − (μ_roll·g/c²)·ln(1 + c·v₀/(μ_roll·g))].
    for surface in [Surface::FAIRWAY, Surface::GREEN] {
        assert_eq!(surface.roll_law, RollLaw::Constant);
        let (mu_g, c) = (surface.mu_roll * G, surface.plow);
        let expected = 1.4 * (v0 / c - mu_g / (c * c) * (c * v0 / mu_g).ln_1p());
        let g = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &surface);
        assert!(
            (g.rest.pos.x - expected).abs() < 1e-5,
            "roll {:.6} m, closed form {expected:.6} m",
            g.rest.pos.x
        );
        assert!((units::meters_to_feet(g.rest.pos.x) - surface.stimp_feet()).abs() < 1e-4);
    }
}

#[test]
fn stimpmeter_fairway_shorter() {
    let v0 = STIMPMETER_SPEED;
    let green = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::GREEN);
    let fairway = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::FAIRWAY);
    assert!(fairway.rest.pos.x < green.rest.pos.x);
}

#[test]
fn preset_stimp_readings() {
    for (surface, feet) in [
        (Surface::GREEN, 10.5),
        (Surface::FIRM_GREEN, 10.5),
        (Surface::FAIRWAY, 5.0 + 8.0 / 12.0),
    ] {
        let stimp = surface.stimp_feet();
        assert!((stimp - feet).abs() < 0.01, "stimp {stimp:.3} ft");
    }
}

/// Level-ground roll (m) from `v0` under a constant rolling resistance
/// `mu_roll` and plowing `c`: (1+k)·[v₀/c − (μg/c²)·ln(1 + c·v₀/(μg))].
fn plowed_roll(mu_roll: f64, c: f64, v0: f64) -> f64 {
    let mu_g = mu_roll * G;
    1.4 * (v0 / c - mu_g / (c * c) * (c * v0 / mu_g).ln_1p())
}

#[test]
fn stimp_feet_with_plow_matches_closed_form() {
    let v0 = STIMPMETER_SPEED;
    for (mu_roll, plow) in [(0.12, 0.36), (0.055, 0.36), (0.12, 0.05), (0.2, 2.0)] {
        let surface = Surface {
            roll_law: RollLaw::Constant,
            mu_roll,
            plow,
            ..Surface::FAIRWAY
        };
        let expected = units::meters_to_feet(plowed_roll(mu_roll, plow, v0));
        assert!(
            (surface.stimp_feet() - expected).abs() < 1e-12,
            "μ {mu_roll}, plow {plow}: {} vs {expected}",
            surface.stimp_feet()
        );
    }
    // Plowing reads shorter than none, and a tiny rate is the unplowed
    // reading.
    let unplowed = Surface {
        plow: 0.0,
        ..Surface::FAIRWAY
    };
    let tiny = Surface {
        plow: 1e-9,
        ..Surface::FAIRWAY
    };
    assert!(Surface::FAIRWAY.stimp_feet() < unplowed.stimp_feet());
    assert!((tiny.stimp_feet() / unplowed.stimp_feet() - 1.0).abs() < 1e-8);
    // With no rolling resistance plowing alone stops the ball within
    // (1+k)·v₀/plow; a reading at or past that has no mu_roll.
    let free = Surface {
        mu_roll: 0.0,
        ..Surface::FAIRWAY
    };
    let longest = units::meters_to_feet(1.4 * v0 / 0.36);
    assert!((free.stimp_feet() - longest).abs() < 1e-12);
    assert_eq!(
        Surface::FAIRWAY.with_stimp(1.001 * longest),
        Surface::FAIRWAY
    );
    assert_eq!(Surface::FAIRWAY.with_stimp(30.0), Surface::FAIRWAY);
}

/// Ground speeds along a result never rise and never turn backward.
fn assert_slows_forward(g: &GroundResult, label: &str) {
    let mut last = f64::INFINITY;
    for p in &g.points {
        let speed = p.vel.length();
        assert!(speed <= last + 1e-12, "{label}: speeds up to {speed} m/s");
        assert!(p.vel.x >= 0.0, "{label}: moves backward {:?}", p.vel);
        last = speed;
    }
}

/// A straight slide from `v0` with no spin on level ground under constant
/// rolling resistance: the slide-to-roll time and the rest distance.
/// v' = −a − c·v with a = `mu_slide` · g; the friction spins the ball up at
/// a/(kR), so the slip v − Rω reaches zero at the t* where v(t*) = a·t*/k;
/// it then rolls from v(t*).
fn slide_then_roll(surface: &Surface, v0: f64) -> (f64, f64) {
    let (a, c) = (surface.mu_slide * G, surface.plow);
    let v = |t: f64| (v0 + a / c) * (-c * t).exp() - a / c;
    let (mut lo, mut hi) = (0.0, v0 / a);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if v(mid) - a * mid / 0.4 > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let t_roll = 0.5 * (lo + hi);
    let slid = (v0 + a / c) * (1.0 - (-c * t_roll).exp()) / c - a * t_roll / c;
    (t_roll, slid + plowed_roll(surface.mu_roll, c, v(t_roll)))
}

#[test]
fn fast_roll_and_slide_on_fairway_slow_to_the_closed_form() {
    let fairway = Surface::FAIRWAY;
    let c = fairway.plow;
    let v0 = 30.0;

    // Rolling from 30 m/s.
    let rolled = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &fairway);
    assert_slows_forward(&rolled, "roll");
    let expected = plowed_roll(fairway.mu_roll, c, v0);
    assert!(
        (rolled.rest.pos.x / expected - 1.0).abs() < 0.01,
        "roll {:.3} m vs {expected:.3} m",
        rolled.rest.pos.x
    );

    // Sliding from 30 m/s with no spin.
    let (t_roll, expected) = slide_then_roll(&fairway, v0);
    let g = roll(Vec3::new(v0, 0.0, 0.0), Vec3::ZERO, &fairway);
    assert_slows_forward(&g, "slide");
    let roll_start = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
    assert!(
        (roll_start.time - t_roll).abs() < 1e-9,
        "rolls from {} s vs {t_roll} s",
        roll_start.time
    );
    assert!(
        (g.rest.pos.x / expected - 1.0).abs() < 0.01,
        "slide and roll {:.3} m vs {expected:.3} m",
        g.rest.pos.x
    );
}

#[test]
fn heavy_plowing_and_sliding_friction_end_the_slide_at_zero_slip() {
    // Off the preset, the slip falls fast enough to pass zero within one
    // step; the slide still ends at the analytic zero-slip time.
    let surfaces = [
        Surface {
            plow: 6.0,
            ..Surface::FAIRWAY
        },
        Surface {
            mu_slide: 0.8,
            ..Surface::FAIRWAY
        },
    ];
    for surface in &surfaces {
        for i in 0..=76 {
            let v0 = 2.0 + 0.5 * f64::from(i);
            let (t_roll, expected) = slide_then_roll(surface, v0);
            let g = roll(Vec3::new(v0, 0.0, 0.0), Vec3::ZERO, surface);
            let roll_start = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
            assert!(
                (roll_start.time - t_roll).abs() < 1e-9,
                "plow {}, μ {}, {v0} m/s: rolls from {} s vs {t_roll} s",
                surface.plow,
                surface.mu_slide,
                roll_start.time
            );
            assert!(
                (g.rest.pos.x - expected).abs() < 1e-5,
                "plow {}, μ {}, {v0} m/s: rests at {} m vs {expected} m",
                surface.plow,
                surface.mu_slide,
                g.rest.pos.x
            );
        }
    }
}

/// The laws `with_stimp` is checked under: constant, and exponential rising
/// and falling with speed.
const LAWS: [RollLaw; 4] = [
    RollLaw::Constant,
    RollLaw::Exponential { v_ref: 1.6 },
    RollLaw::Exponential { v_ref: 7.0 },
    RollLaw::Exponential { v_ref: -4.7 },
];

#[test]
fn with_stimp_round_trips_under_every_law() {
    for law in LAWS {
        for base in [Surface::GREEN, Surface::FAIRWAY] {
            let surface = Surface {
                roll_law: law,
                ..base
            };
            for feet in [3.0, 6.0, 8.5, 10.0, 11.4, 12.0, 13.5, 14.0, 20.0] {
                let tuned = surface.with_stimp(feet);
                assert!(
                    (tuned.stimp_feet() - feet).abs() < 1e-9,
                    "{law:?}: {feet} ft reads {} ft",
                    tuned.stimp_feet()
                );
                assert_eq!(
                    tuned,
                    Surface {
                        mu_roll: tuned.mu_roll,
                        ..surface
                    },
                    "only mu_roll changes"
                );
            }
        }
    }
    // Faster greens have less rolling resistance.
    assert!(Surface::GREEN.with_stimp(13.0).mu_roll < Surface::GREEN.with_stimp(9.0).mu_roll);
    // Constant rolling resistance and no plowing: d = (1+k)·v₀²/(2·μ·g).
    let with_law = |roll_law| Surface {
        roll_law,
        ..TEST_TURF
    };
    let constant = with_law(RollLaw::Constant).with_stimp(10.0);
    let expected =
        1.4 * STIMPMETER_SPEED * STIMPMETER_SPEED / (2.0 * G * units::feet_to_meters(10.0));
    assert!((constant.mu_roll - expected).abs() < 1e-15);
    // The exponential law with an infinite scale is the constant law.
    for law in [
        RollLaw::Exponential {
            v_ref: f64::INFINITY,
        },
        RollLaw::Exponential {
            v_ref: f64::NEG_INFINITY,
        },
    ] {
        assert_eq!(
            with_law(law).with_stimp(10.0).mu_roll,
            constant.mu_roll,
            "{law:?}"
        );
    }
    // A very long scale approaches it smoothly (the small-s₀ series).
    for v_ref in [1e3, -1e3, 1e8] {
        let mu = with_law(RollLaw::Exponential { v_ref })
            .with_stimp(10.0)
            .mu_roll;
        let expected = constant.mu_roll * (1.0 - 2.0 * STIMPMETER_SPEED / (3.0 * v_ref));
        assert!(
            (mu - expected).abs() < 1e-6 * constant.mu_roll,
            "v_ref {v_ref}: {mu}"
        );
    }
}

#[test]
fn with_stimp_rolls_the_requested_distance() {
    // A ball leaving the stimpmeter, stepped with `GroundPhysics::step`
    // on level ground until it stops, and the same roll through `settle`.
    // Without plowing the constant law is integrated exactly; with
    // `GREEN`'s plowing the RK4 steps carry a small error.
    let v0 = STIMPMETER_SPEED;
    for (base, tolerance) in [(TEST_TURF, 1e-9), (Surface::GREEN, 1e-7)] {
        for feet in [6.0, 8.0, 10.0, 12.0, 14.0] {
            let surface = base.with_stimp(feet);
            let physics = GroundPhysics::new(surface);
            let mut s = state(Vec3::new(v0, 0.0, 0.0), rolling(v0));
            while s.vel != Vec3::ZERO {
                assert!(s.time < 30.0, "{feet} ft: still rolling");
                s = physics.step(&s, &BALL).0;
            }
            let stepped = units::meters_to_feet(s.pos.x);
            assert!(
                (stepped - feet).abs() < tolerance * feet,
                "{feet} ft: stepped roll {stepped:.9} ft"
            );
            let settled = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &surface);
            assert_eq!(settled.termination, Termination::Rest);
            let settled = units::meters_to_feet(settled.rest.pos.x);
            assert!(
                (settled - feet).abs() < tolerance * feet,
                "{feet} ft: settled roll {settled:.9} ft"
            );
        }
    }
}

#[test]
fn with_stimp_ignores_invalid_readings() {
    for feet in [0.0, -0.0, -8.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(Surface::GREEN.with_stimp(feet), Surface::GREEN, "{feet}");
    }
    let frictionless = Surface {
        mu_roll: 0.0,
        ..TEST_TURF
    };
    assert_eq!(frictionless.stimp_feet(), f64::INFINITY);
}

#[test]
fn invalid_v_ref_leaves_with_stimp_unchanged() {
    for v_ref in [0.0, -0.0, f64::NAN] {
        let roll_law = RollLaw::Exponential { v_ref };
        let surface = Surface {
            roll_law,
            ..Surface::GREEN
        };
        let tuned = surface.with_stimp(10.0);
        assert_eq!(tuned.mu_roll, surface.mu_roll, "{roll_law:?}");
        assert!(surface.stimp_feet().is_nan(), "{roll_law:?}");
    }
}

#[test]
fn invalid_v_ref_rolls_as_the_constant_law() {
    // Outside the documented range a zero or NaN scale gives f = 1 rather
    // than a NaN that would stall the roll until the time limit.
    let v0 = 3.0;
    let constant = Surface {
        roll_law: RollLaw::Constant,
        ..Surface::GREEN
    };
    let expected = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &constant);
    for v_ref in [0.0, f64::NAN] {
        let roll_law = RollLaw::Exponential { v_ref };
        assert_eq!(roll_law.speed_factor(0.0), 1.0, "{roll_law:?}");
        assert_eq!(roll_law.speed_factor(v0), 1.0, "{roll_law:?}");
        let surface = Surface {
            roll_law,
            ..Surface::GREEN
        };
        let g = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &surface);
        assert_eq!(g.termination, Termination::Rest, "{roll_law:?}");
        assert_eq!(g.rest.pos, expected.rest.pos, "{roll_law:?}");
    }
}

#[test]
fn roll_laws_agree_at_zero_speed() {
    for law in LAWS {
        assert_eq!(law.speed_factor(0.0), 1.0, "{law:?}");
        let surface = Surface {
            roll_law: law,
            ..Surface::GREEN
        };
        assert_eq!(surface.rolling_friction(0.0), surface.mu_roll);
    }
    // Away from zero each law has its own shape.
    let v = 2.0;
    assert_eq!(RollLaw::Constant.speed_factor(v), 1.0);
    let e = RollLaw::Exponential { v_ref: 4.0 }.speed_factor(v);
    assert!((e - 0.5_f64.exp()).abs() < 1e-12);
    assert!(RollLaw::Exponential { v_ref: -4.0 }.speed_factor(v) < 1.0);
}

#[test]
fn exponential_rollout_matches_closed_form() {
    // a = μ₀·e^(v/v_ref)·g/(1+k), so a ball rolling from v₀ stops after
    // (1+k)·v_ref²·(1 − e^(−s₀)(1 + s₀))/(μ₀·g), s₀ = v₀/v_ref.
    let v0 = STIMPMETER_SPEED;
    for v_ref in [1.6, 3.8, 7.0, -4.7] {
        let surface = Surface {
            roll_law: RollLaw::Exponential { v_ref },
            ..TEST_TURF
        };
        let g = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &surface);
        assert_eq!(g.termination, Termination::Rest);
        let s0 = v0 / v_ref;
        let expected =
            1.4 * v_ref * v_ref * (1.0 - (-s0).exp() * (1.0 + s0)) / (surface.mu_roll * G);
        assert!(
            (g.rest.pos.x - expected).abs() < 1e-5,
            "v_ref {v_ref}: roll {:.6} m, closed form {expected:.6} m",
            g.rest.pos.x
        );
        let feet = units::meters_to_feet(g.rest.pos.x);
        assert!((feet - surface.stimp_feet()).abs() < 1e-4, "v_ref {v_ref}");
    }
}

#[test]
fn constant_law_reads_penner_stimp() {
    // Penner's (2002) constant rolling resistance ρ_g = 0.131 is a ~6.0 ft
    // green.
    let penner = Surface {
        mu_roll: 0.131,
        ..TEST_TURF
    };
    let stimp = penner.stimp_feet();
    assert!((stimp - 6.0).abs() < 0.05, "stimp {stimp:.3} ft");
}

#[test]
fn exponential_with_positive_rate_shortens_fast_rollouts() {
    // Same stimp, so the same roll from the stimpmeter; from 5 m/s the
    // exponential law resists more and stops sooner, nearer.
    let feet = Surface::GREEN.stimp_feet();
    let with_law = |roll_law| {
        Surface {
            roll_law,
            ..Surface::GREEN
        }
        .with_stimp(feet)
    };
    let v0 = 5.0;
    let constant = roll(
        Vec3::new(v0, 0.0, 0.0),
        rolling(v0),
        &with_law(RollLaw::Constant),
    );
    let exponential = roll(
        Vec3::new(v0, 0.0, 0.0),
        rolling(v0),
        &with_law(RollLaw::Exponential { v_ref: 3.8 }),
    );
    assert!(exponential.rest.time < constant.rest.time);
    assert!(exponential.rest.pos.x < constant.rest.pos.x);
}

#[test]
fn rollout_stops_and_moves_forward() {
    let g = roll(Vec3::new(10.0, 0.0, 0.0), rolling(10.0), &Surface::FAIRWAY);
    assert_eq!(g.termination, Termination::Rest);
    assert_eq!(g.rest.vel, Vec3::ZERO);
    assert!(g.rest.pos.x > 0.0);
    assert_eq!(
        g.contacts.len(),
        1,
        "a ball already on the ground has one contact"
    );
    assert_eq!(g.contacts[0].vel_out, g.contacts[0].vel_in);
}

#[test]
fn rollout_on_ground_plane() {
    let g = roll(Vec3::new(10.0, 0.0, 2.0), Vec3::ZERO, &Surface::GREEN);
    assert!(g.points.iter().all(|p| p.pos.y == 0.0 && p.vel.y == 0.0));
}

#[test]
fn speed_monotonically_decreases_pure_roll() {
    let g = roll(Vec3::new(8.0, 0.0, 0.0), rolling(8.0), &Surface::GREEN);
    assert!(g.points.iter().all(|p| p.phase == Phase::Roll));
    for w in g.points.windows(2) {
        assert!(w[1].vel.length() <= w[0].vel.length());
    }
}

#[test]
fn energy_never_increases_over_the_ground_phase() {
    // Translation + rotation + height never grows from one ground point to
    // the next, through hops, impacts, slides and rolls.
    let starts = [
        (Vec3::new(21.3, -18.5, 0.0), backspin(212.0)),
        (Vec3::new(21.0, -18.0, 2.5), Vec3::new(0.0, -50.0, 210.0)),
        (Vec3::new(9.0, -20.0, -1.0), Vec3::new(10.0, 40.0, 700.0)),
        (Vec3::new(3.0, 0.0, 0.5), Vec3::new(-20.0, 30.0, 300.0)),
    ];
    for (vel, spin) in starts {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let g = surface.settle(&state(vel, spin), air());
            let mut prev = energy(vel, spin, 0.0);
            for p in &g.points {
                let e = energy(p.vel, p.spin, p.pos.y);
                assert!(
                    e <= prev + 1e-9 * prev.max(1.0),
                    "energy rose {prev:.9} → {e:.9} at t={:.3} ({:?})",
                    p.time,
                    p.phase
                );
                prev = e;
            }
        }
    }
}

#[test]
fn green_rolls_farther_than_fairway() {
    let fairway = roll(Vec3::new(5.0, 0.0, 0.0), rolling(5.0), &Surface::FAIRWAY);
    let green = roll(Vec3::new(5.0, 0.0, 0.0), rolling(5.0), &Surface::GREEN);
    assert!(green.rest.pos.x > fairway.rest.pos.x);
}

#[test]
fn fairway_rolls_less_than_green_from_the_same_landing() {
    // The presets bounce a little differently, so compare the roll from each
    // one's last impact: the fairway's higher rolling resistance stops the
    // ball in 0.51–0.61 of the green's roll. Landings from a lob wedge, a
    // pitch, a low chip, a 7-iron and a drive.
    for (speed, descent, rpm) in [
        (19.8, 44.7, 5651.0),
        (17.6, 39.1, 3599.0),
        (12.4, 19.1, 2413.0),
        (24.2, 49.5, 5611.0),
        (28.8, 39.7, 1914.0),
    ] {
        let a = units::deg_to_rad(descent);
        let landing = state(
            Vec3::new(speed * a.cos(), -speed * a.sin(), 0.0),
            backspin(rpm_to_rads(rpm)),
        );
        let roll_out = |g: &GroundResult| g.rest.pos.x - g.contacts.last().unwrap().pos.x;
        let f = roll_out(&Surface::FAIRWAY.settle(&landing, air()));
        let g = roll_out(&Surface::GREEN.settle(&landing, air()));
        assert!(
            f > 0.0 && f < 0.7 * g,
            "{speed} m/s at {descent}°: fairway {f:.2} m, green {g:.2} m"
        );
    }
}

#[test]
fn backspin_reduces_rollout() {
    let no_spin = roll(Vec3::new(5.0, 0.0, 0.0), rolling(5.0), &Surface::GREEN);
    let with_spin = roll(
        Vec3::new(5.0, 0.0, 0.0),
        backspin(rpm_to_rads(5000.0)),
        &Surface::GREEN,
    );
    assert!(with_spin.rest.pos.x < no_spin.rest.pos.x);
}

#[test]
fn extreme_backspin_can_reverse() {
    let g = roll(
        Vec3::new(2.0, 0.0, 0.0),
        backspin(rpm_to_rads(8000.0)),
        &Surface::GREEN,
    );
    assert!(g.rest.pos.x < 0.0, "should reverse: x={:.3}", g.rest.pos.x);
}

#[test]
fn driver_rollout_reasonable() {
    // ~10 m/s forward, nearly rolling: the speed-dependent rolling
    // resistance keeps a firm-fairway roll short.
    let g = roll(Vec3::new(10.0, 0.0, 0.0), rolling(8.0), &Surface::FAIRWAY);
    let yards = units::meters_to_yards(g.rest.pos.x);
    assert!(yards > 5.0 && yards < 60.0, "rollout {yards:.1} yds");
}

#[test]
fn ground_physics_step() {
    let physics = GroundPhysics::new(Surface::GREEN);
    let sliding = state(Vec3::new(3.0, 0.0, 0.0), backspin(300.0));
    let (next, phase) = physics.step(&sliding, &BALL);
    assert_eq!(phase, Phase::Slide);
    assert!((next.time - 0.001).abs() < 1e-15 && next.vel.x < 3.0);

    let rolling_ball = state(Vec3::new(3.0, 0.0, 0.0), rolling(3.0));
    let (next, phase) = physics.step(&rolling_ball, &BALL);
    assert_eq!(phase, Phase::Roll);
    assert!((next.spin.z + next.vel.x / R).abs() < 1e-9);

    // A ball slower than one step's deceleration stops inside the step.
    let creeping = state(Vec3::new(1e-4, 0.0, 0.0), rolling(1e-4));
    let (next, phase) = physics.step(&creeping, &BALL);
    assert_eq!(phase, Phase::Roll);
    assert_eq!(next.vel, Vec3::ZERO);
    assert!(next.time < 0.001 && next.pos.x > 0.0);

    // The step uses the instance's tunables: its `dt`, and its
    // `stop_speed`, which here stops a ball one default step would not.
    let tuned = GroundPhysics {
        dt: 0.002,
        stop_speed: 0.5,
        ..physics
    };
    let (next, _) = tuned.step(&rolling_ball, &BALL);
    assert!((next.time - 0.002).abs() < 1e-15);
    let slow = state(Vec3::new(0.3, 0.0, 0.0), rolling(0.3));
    assert_ne!(physics.step(&slow, &BALL).0.vel, Vec3::ZERO);
    assert_eq!(tuned.step(&slow, &BALL).0.vel, Vec3::ZERO);

    // It samples its own terrain under the ball.
    let split = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| {
        if x < 1.0 {
            Surface::GREEN.sample(x, 0.0)
        } else {
            Surface::FAIRWAY.sample(x, 0.0)
        }
    }));
    let far = BallState {
        pos: Vec3::new(5.0, 0.0, 0.0),
        ..rolling_ball
    };
    let on_fairway = GroundPhysics::new(Surface::FAIRWAY).step(&far, &BALL);
    assert_eq!(split.step(&far, &BALL), on_fairway);
    assert_eq!(
        split.step(&rolling_ball, &BALL),
        physics.step(&rolling_ball, &BALL)
    );
}

/// The tour ball in still sea-level air with inertia factor `.0`.
struct Hollow(f64);

impl AeroModel for Hollow {
    fn acceleration(&self, ball: &BallState, height: f64) -> Vec3 {
        STILL_AIR.acceleration(ball, height)
    }
    fn spin_after(&self, ball: &BallState, dt: f64) -> Vec3 {
        STILL_AIR.spin_after(ball, dt)
    }
    fn radius_m(&self) -> f64 {
        STILL_AIR.radius_m()
    }
    fn inertia_factor(&self) -> f64 {
        self.0
    }
}

#[test]
fn hollow_ball_follows_generalized_formulas() {
    // A thin-shelled ball, k = 2/3, follows the general-k forms of the
    // grip kick, the slide-to-roll time and the stimpmeter distance.
    let k = 2.0 / 3.0;
    let hollow = Hollow(k);

    // Grip on rigid ground: heading atan(v_z / (v_x − k·R·ω_b)).
    let rigid = plane(0.5, 1.0, 0.0);
    for (vel, w) in [
        (Vec3::new(20.0, -12.0, 3.0), 250.0),
        (Vec3::new(15.0, -14.0, -2.0), 300.0),
    ] {
        let out = impact(&state(vel, backspin(w)), &rigid.sample(0.0, 0.0), &hollow);
        let expected = vel.z.atan2(vel.x - k * R * w);
        assert!((out.vel.z.atan2(out.vel.x) - expected).abs() < 1e-12);
        assert!(slip(out.vel, out.spin, Vec3::Y, R).length() < 1e-12);
    }

    // Slide to roll after |u₀| / (((1+k)/k)·μ·g).
    let vel = Vec3::new(3.0, 0.0, 0.5);
    let spin = Vec3::new(-20.0, 30.0, 300.0);
    let expected = slip(vel, spin, Vec3::Y, R).length() / ((1.0 + k) / k * 0.40 * G);
    let g = TEST_TURF.settle(&state(vel, spin), &hollow);
    let roll_start = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
    assert!((roll_start.time - expected).abs() < 1e-12);

    // Stimpmeter, constant law without plowing: (1+k)·v₀²/(2·μ_roll·g).
    let v0 = STIMPMETER_SPEED;
    let start = state(Vec3::new(v0, 0.0, 0.0), rolling(v0));
    let expected = (1.0 + k) * v0 * v0 / (2.0 * TEST_TURF.mu_roll * G);
    let g = TEST_TURF.settle(&start, &hollow);
    assert!(
        (g.rest.pos.x - expected).abs() < 1e-5,
        "roll {:.6} m, closed form {expected:.6} m",
        g.rest.pos.x
    );
    // Constant law with plowing c (FAIRWAY):
    // (1+k)·[v₀/c − (μ_roll·g/c²)·ln(1 + c·v₀/(μ_roll·g))].
    let fairway = Surface::FAIRWAY;
    assert_eq!(fairway.roll_law, RollLaw::Constant);
    let (mu_g, c) = (fairway.mu_roll * G, fairway.plow);
    let expected = (1.0 + k) * (v0 / c - mu_g / (c * c) * (c * v0 / mu_g).ln_1p());
    let g = fairway.settle(&start, &hollow);
    assert!(
        (g.rest.pos.x - expected).abs() < 1e-5,
        "roll {:.6} m, closed form {expected:.6} m",
        g.rest.pos.x
    );
}

#[test]
fn ground_point_times_are_monotonic() {
    // Shortened slide steps restart the 10 ms sample grid at the switch to
    // rolling; the points stay in time order, never more than one sample
    // interval apart, and end at the final state, for every exit. The last
    // gap can be longer: a ball below `stop_speed` finishes at constant
    // deceleration in one step (up to ~26 ms on `GREEN`).
    let mut results = vec![
        roll(
            Vec3::new(3.0, 0.0, 0.5),
            Vec3::new(-20.0, 30.0, 300.0),
            &Surface::GREEN,
        ),
        roll(
            Vec3::new(2.0, 0.0, 0.0),
            backspin(rpm_to_rads(8000.0)),
            &Surface::GREEN,
        ),
        GroundPhysics {
            max_time: 0.0237,
            ..GroundPhysics::new(Surface::GREEN)
        }
        .settle(&state(Vec3::new(3.0, 0.0, 0.0), backspin(300.0)), air()),
        GroundPhysics {
            max_impacts: 3,
            ..GroundPhysics::new(Surface::FAIRWAY)
        }
        .settle(&state(Vec3::new(21.3, -18.5, 0.0), backspin(212.0)), air()),
    ];
    for (vel, spin) in [
        (Vec3::new(21.3, -18.5, 0.0), backspin(212.0)),
        (Vec3::new(21.0, -18.0, 2.5), Vec3::new(0.0, -50.0, 210.0)),
        (Vec3::new(9.0, -20.0, -1.0), Vec3::new(10.0, 40.0, 700.0)),
    ] {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            results.push(surface.settle(&state(vel, spin), air()));
        }
    }
    for g in &results {
        let n = g.points.len();
        for (i, w) in g.points.windows(2).enumerate() {
            let gap = w[1].time - w[0].time;
            let max_gap = if i + 2 == n {
                f64::INFINITY
            } else {
                SAMPLE_INTERVAL + TIME_EPS
            };
            assert!(
                (0.0..=max_gap).contains(&gap),
                "gap {gap:.6} s at t={:.4} ({:?})",
                w[0].time,
                g.termination
            );
        }
        let end = g.points.last().unwrap();
        assert!(end.time == g.rest.time && end.pos == g.rest.pos);
    }
}

// ══════════════════════════════════════════════════════════════════
//  TERRAIN — per-point surfaces and limits
// ══════════════════════════════════════════════════════════════════

#[test]
fn roll_from_fairway_onto_green_ends_between() {
    // Rolling from fairway onto a green that starts 2 m ahead: the ball
    // feels the green mid-roll and stops between the all-fairway and
    // all-green distances. At 3 m/s the presets' green resists less than
    // their fairway at every speed the ball passes through.
    let v0 = 3.0;
    let mixed = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| GroundSample {
        height: 0.0,
        normal: Vec3::Y,
        surface: if x < 2.0 {
            Surface::FAIRWAY
        } else {
            Surface::GREEN
        },
    }));
    let start = state(Vec3::new(v0, 0.0, 0.0), rolling(v0));
    let g = mixed.settle(&start, air());
    let fairway = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::FAIRWAY);
    let green = roll(Vec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::GREEN);
    assert_eq!(g.termination, Termination::Rest);
    assert!(
        fairway.rest.pos.x < g.rest.pos.x && g.rest.pos.x < green.rest.pos.x,
        "fairway {:.2} < mixed {:.2} < green {:.2}",
        fairway.rest.pos.x,
        g.rest.pos.x,
        green.rest.pos.x
    );
    assert!(g.rest.pos.x > 2.0);
    assert_eq!(g.contacts[0].surface, Surface::FAIRWAY);
}

#[test]
fn landing_on_green_from_fn_terrain_uses_green() {
    // The surface at each contact comes from the terrain at that point.
    let terrain = FnTerrain(|x: f64, _z: f64| GroundSample {
        height: 0.0,
        normal: Vec3::Y,
        surface: if x > 0.5 {
            Surface::GREEN
        } else {
            Surface::FAIRWAY
        },
    });
    let g = GroundPhysics::new(terrain)
        .settle(&state(Vec3::new(10.0, -8.0, 0.0), backspin(300.0)), air());
    assert_eq!(g.contacts[0].surface, Surface::FAIRWAY);
    assert!(g.contacts[1..].iter().all(|c| c.surface == Surface::GREEN));
}

#[test]
fn long_roll_reports_max_time() {
    // No rolling resistance: the ball would roll forever.
    let physics = GroundPhysics {
        max_time: 2.0,
        ..GroundPhysics::new(Surface {
            mu_roll: 0.0,
            ..TEST_TURF
        })
    };
    let g = physics.settle(&state(Vec3::new(3.0, 0.0, 0.0), rolling(3.0)), air());
    assert_eq!(g.termination, Termination::MaxTime);
    assert!((g.rest.time - 2.0).abs() < 1e-9);
    assert!((g.rest.pos.x - 6.0).abs() < 1e-6);
    assert!(g.rest.vel.x > 2.99);
}

#[test]
fn ground_physics_tunables() {
    let start = state(Vec3::new(15.0, -10.0, 0.0), backspin(250.0));
    let default = GroundPhysics::new(Surface::FAIRWAY).settle(&start, air());
    assert_eq!(default, Surface::FAIRWAY.settle(&start, air()));
    // A higher hop transition ends the bouncing sooner.
    let fewer = GroundPhysics {
        min_hop_height: 0.05,
        ..GroundPhysics::new(Surface::FAIRWAY)
    }
    .settle(&start, air());
    assert!(fewer.contacts.len() < default.contacts.len());
    // A coarser step barely moves the rest position.
    let coarse = GroundPhysics {
        dt: 0.004,
        ..GroundPhysics::new(Surface::FAIRWAY)
    }
    .settle(&start, air());
    assert!((coarse.rest.pos - default.rest.pos).length() < 1e-3);
}

#[test]
fn default_landing_interpolates_crossing() {
    let a = BallState {
        time: 1.0,
        pos: Vec3::new(10.0, 0.3, 1.0),
        vel: Vec3::new(10.0, -3.0, 1.0),
        spin: backspin(100.0),
    };
    let b = BallState {
        time: 1.1,
        pos: Vec3::new(11.0, -0.1, 1.1),
        vel: Vec3::new(10.0, -4.0, 1.0),
        spin: backspin(90.0),
    };
    let l = Surface::GREEN.landing(&a, &b).unwrap();
    assert!((l.time - 1.075).abs() < 1e-12);
    assert!((l.pos - Vec3::new(10.75, 0.0, 1.075)).length() < 1e-12 && l.pos.y == 0.0);
    assert!((l.vel.y + 3.75).abs() < 1e-12 && (l.spin.z - 92.5).abs() < 1e-12);
    // Leaving the ground is not a landing.
    assert!(Surface::GREEN.landing(&b, &a).is_none());
    let tee = BallState {
        pos: Vec3::ZERO,
        ..a
    };
    assert!(Surface::GREEN.landing(&tee, &a).is_none());
}
