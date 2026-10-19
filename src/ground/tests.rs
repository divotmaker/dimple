use super::regression_data::REFERENCES;
use super::*;
use crate::units;

const BALL: BallModel = BallModel::TOUR;
const R: f64 = 0.04267 / 2.0;

fn rpm_to_rads(rpm: f64) -> f64 {
    units::rpm_to_rads(rpm)
}

fn state(vel: DVec3, spin: DVec3) -> BallState {
    BallState {
        time: 0.0,
        pos: DVec3::ZERO,
        vel,
        spin,
    }
}

/// Spin about the across-track axis for a ball moving along +X:
/// positive = backspin.
fn backspin(w: f64) -> DVec3 {
    DVec3::new(0.0, 0.0, w)
}

/// Impact on flat `surface` with incoming `vel` and spin.
fn hit(vel: DVec3, spin: DVec3, surface: &Surface) -> BallState {
    impact(&state(vel, spin), &surface.sample(0.0, 0.0), &BALL)
}

/// Impact at `speed` m/s, `from_vertical_deg` from the vertical, moving +X.
fn penner_impact(speed: f64, from_vertical_deg: f64) -> DVec3 {
    let phi = units::deg_to_rad(from_vertical_deg);
    DVec3::new(speed * phi.sin(), -speed * phi.cos(), 0.0)
}

/// Mechanical energy per unit mass: translation + rotation + height.
fn energy(vel: DVec3, spin: DVec3, y: f64) -> f64 {
    let k = BALL.inertia_factor;
    0.5 * vel.length_squared() + 0.5 * k * R * R * spin.length_squared() + G * y
}

/// Ground phase only: a ball already on flat ground moving horizontally.
fn roll(vel: DVec3, spin: DVec3, surface: &Surface) -> GroundResult {
    surface.settle(&state(vel, spin), &BALL)
}

/// Rolling spin for a ball moving along +X at `v` (no slip).
fn rolling(v: f64) -> DVec3 {
    DVec3::new(0.0, 0.0, -v / R)
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
    // Zero sidespin, zero azimuth, flat ground: the slip is along the
    // track and the spin has no component about the tilted normal, so the
    // vector impact must reproduce the planar (scalar) solution to float
    // tolerance, impact by impact.
    for r in REFERENCES {
        let landing = BallState {
            time: 0.0,
            pos: DVec3::from(r.landing_pos),
            vel: DVec3::from(r.landing_vel),
            spin: backspin(r.backspin),
        };
        let g = r.surface.settle(&landing, &BALL);
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
            let mut expected = DVec3::from(vel);
            let mut expected_spin = backspin(-omega);
            if i == last {
                // The last rebound is below the hop threshold, so the ball
                // stays on the ground, and `vel_out`/`spin_out` are what it
                // slides or rolls away with: the reference rebound projected
                // onto the ground, with a slip at or below the roll
                // threshold removed. The reference stops at the projection.
                expected.y = 0.0;
                if slip(expected, expected_spin, DVec3::Y, R).length() <= DEFAULT_ROLL_THRESHOLD {
                    (expected, expected_spin) =
                        lock_rolling(expected, expected_spin, DVec3::Y, R, BALL.inertia_factor);
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
                (c.pos - DVec3::from(pos)).length() <= 1e-8,
                "{label} impact {i}: pos {:?}, reference {pos:?}",
                c.pos
            );
        }
        // The rest positions differ from the reference by 1–16 mm (largest
        // for the driver on GREEN), mainly because the reference integrates
        // the slide and roll with a first-order Euler step and this crate
        // uses RK4 with an exact slide-to-roll switch; the different stop
        // thresholds account for at most 2 mm.
        let rest = DVec3::from(r.rest);
        assert!(
            (g.rest.pos - rest).length() < 0.025,
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
    let r = hit(DVec3::new(0.0, -10.0, 0.0), DVec3::ZERO, &Surface::FAIRWAY);
    assert!(r.vel.y > 0.0, "rebound should go up: vy={:.3}", r.vel.y);
}

#[test]
fn impact_keeps_time_and_position() {
    let s = BallState {
        time: 4.2,
        pos: DVec3::new(100.0, 0.0, -3.0),
        vel: DVec3::new(20.0, -15.0, 1.0),
        spin: DVec3::new(10.0, -60.0, 300.0),
    };
    let out = impact(&s, &Surface::GREEN.sample(0.0, 0.0), &BALL);
    assert!(out.time == s.time && out.pos == s.pos);
}

#[test]
fn no_lateral_without_sidespin() {
    let r = hit(
        DVec3::new(20.0, -15.0, 0.0),
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
    let r = hit(DVec3::new(0.0, -5.0, 0.0), DVec3::ZERO, &Surface::GREEN);
    assert!(r.vel.y > 0.0);
    assert!(r.vel.x.abs() < 1e-12 && r.vel.z.abs() < 1e-12);
}

#[test]
fn grazing_impact_applies_no_impulse() {
    for vel in [DVec3::new(10.0, 0.0, 0.0), DVec3::new(10.0, 0.5, -2.0)] {
        let s = state(vel, DVec3::new(30.0, -80.0, 250.0));
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
        DVec3::new(20.0, -15.0, 0.0),
        DVec3::new(25.0, -10.0, 4.0),
        DVec3::new(-5.0, -3.0, -2.0),
        DVec3::new(3.0, -20.0, 0.5),
        DVec3::new(0.0, -8.0, 0.0),
    ];
    let spins = [
        DVec3::ZERO,
        backspin(800.0),
        DVec3::new(0.0, -100.0, 300.0),
        DVec3::new(200.0, 50.0, -100.0),
        DVec3::new(-150.0, 0.0, 0.0),
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
    let mz = |v: DVec3| DVec3::new(v.x, v.y, -v.z);
    let mz_spin = |w: DVec3| DVec3::new(-w.x, -w.y, w.z);
    let mx = |v: DVec3| DVec3::new(-v.x, v.y, v.z);
    let mx_spin = |w: DVec3| DVec3::new(w.x, -w.y, -w.z);
    for (vel, spin) in [
        (DVec3::new(5.0, -3.0, 0.0), backspin(-5.0 / R)),
        (DVec3::new(15.0, -10.0, 2.0), DVec3::new(20.0, -60.0, 300.0)),
        (DVec3::new(23.8, -16.0, -1.0), DVec3::new(0.0, 80.0, 900.0)),
    ] {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let fwd = hit(vel, spin, &surface);
            let z = hit(mz(vel), mz_spin(spin), &surface);
            assert!((z.vel - mz(fwd.vel)).length() < 1e-12);
            assert!((z.spin - mz_spin(fwd.spin)).length() < 1e-9);
            let x = hit(mx(vel), mx_spin(spin), &surface);
            assert!((x.vel - mx(fwd.vel)).length() < 1e-12);
            assert!((x.spin - mx_spin(fwd.spin)).length() < 1e-9);
        }
    }
}

#[test]
fn impact_frame_invariance_about_y() {
    let a = 0.7_f64;
    let rot = |v: DVec3| {
        DVec3::new(
            v.x * a.cos() - v.z * a.sin(),
            v.y,
            v.x * a.sin() + v.z * a.cos(),
        )
    };
    for (vel, spin) in [
        (DVec3::new(20.0, -15.0, 1.5), DVec3::new(10.0, -60.0, 250.0)),
        (DVec3::new(6.0, -2.0, -0.5), DVec3::new(-30.0, 40.0, -200.0)),
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
fn green_rebound_steeper_and_slower_than_fairway() {
    // Both presets share Penner's restitution and friction; the green's
    // larger compliance tilts the effective contact plane further, which
    // steepens the rebound and takes more forward speed out of it.
    let vel = DVec3::new(15.0, -10.0, 0.0);
    let spin = backspin(rpm_to_rads(3000.0));
    let rf = hit(vel, spin, &Surface::FAIRWAY);
    let rg = hit(vel, spin, &Surface::GREEN);
    assert!(rg.vel.x < rf.vel.x);
    assert!(rg.vel.y.atan2(rg.vel.x) > rf.vel.y.atan2(rf.vel.x));
}

#[test]
fn backspin_reduces_forward_speed() {
    let vel = DVec3::new(15.0, -10.0, 0.0);
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
    // Wedge shot: steep angle, very high backspin on green → spin back.
    let angle = units::deg_to_rad(50.0);
    let vel = 25.0 * DVec3::new(angle.cos(), -angle.sin(), 0.0);
    let r = hit(vel, backspin(rpm_to_rads(10000.0)), &Surface::GREEN);
    assert!(r.vel.x < 0.0, "should spin back: vx={:.3}", r.vel.x);
}

#[test]
fn driver_bounces_forward() {
    let angle = units::deg_to_rad(30.0);
    let vel = 30.0 * DVec3::new(angle.cos(), -angle.sin(), 0.0);
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
    let vel = 25.0 * DVec3::new(angle.cos(), -angle.sin(), 0.0);
    let spin = hit(vel, backspin(rpm_to_rads(7000.0)), &Surface::GREEN);
    let no_spin = hit(vel, DVec3::ZERO, &Surface::GREEN);
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
    let vel = DVec3::new(20.0, -15.0, 0.0);
    let spin = backspin(rpm_to_rads(5000.0));
    let firm = hit(vel, spin, &plane(0.35, 0.4, 0.2));
    let soft = hit(vel, spin, &plane(0.35, 0.4, 1.0));
    let angle = |v: DVec3| v.y.atan2(v.x.abs().max(0.001));
    assert!(angle(soft.vel) > angle(firm.vel));
}

#[test]
fn compliance_rotates_mirror_rebound_by_two_beta() {
    // Elastic, frictionless impact at 45°: the effective plane is tilted by
    // Penner's θc = 15.4° · (v / 18.6 m/s) · (φ / 44.4°) toward the ball, so
    // the specular rebound is the flat-ground mirror direction rotated up by
    // 2β, at full speed. The same holds for any heading.
    for heading in [0.0_f64, 0.4, -2.0] {
        let h = DVec3::new(heading.cos(), 0.0, heading.sin());
        let vel = 10.0 * h - DVec3::new(0.0, 10.0, 0.0);
        let r = hit(vel, DVec3::ZERO, &plane(1.0, 0.0, 1.0));
        let incidence = std::f64::consts::FRAC_PI_4;
        let beta =
            units::deg_to_rad(15.4) * (vel.length() / 18.6) * (incidence / units::deg_to_rad(44.4));
        let rebound_angle = r.vel.y.atan2(r.vel.dot(h));
        assert!((rebound_angle - (incidence + 2.0 * beta)).abs() < 1e-12);
        assert!((r.vel.length() - vel.length()).abs() < 1e-12);
        assert!(
            r.vel.dot(DVec3::Y.cross(h)).abs() < 1e-12,
            "stays in the plane of travel"
        );
    }
}

#[test]
fn no_compliance_is_exact_mirror() {
    let r = hit(
        DVec3::new(10.0, -10.0, 0.0),
        DVec3::ZERO,
        &plane(1.0, 0.0, 0.0),
    );
    assert!((r.vel - DVec3::new(10.0, 10.0, 0.0)).length() < 1e-12);
}

#[test]
fn compliance_tilt_capped_at_incidence_angle() {
    // 10 × θc ≈ 119° would overshoot; β is capped at the 45° incidence
    // angle, so the ball meets the effective plane head-on and an elastic
    // frictionless impact sends it straight back.
    let vel = DVec3::new(10.0, -10.0, 0.0);
    let r = hit(vel, DVec3::ZERO, &plane(1.0, 0.0, 10.0));
    assert!((r.vel + vel).length() < 1e-9, "got {:?}", r.vel);
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
    let vel = DVec3::new(10.0, -10.0, 0.0);
    let raw = hit(vel, backspin(1500.0), &surface);
    assert!(raw.vel.y < 0.0, "impulse alone points into the ground");
    let g = surface.settle(&state(vel, backspin(1500.0)), &BALL);
    assert_eq!(g.contacts.len(), 1);
    assert!(g.contacts[0].vel_out.y == 0.0);
    assert!(g.contacts[0].vel_out.x < 0.0);
    assert!(g.rest.pos.x < 0.0, "rolls back: x={:.3}", g.rest.pos.x);
}

// ══════════════════════════════════════════════════════════════════
//  PENNER (2002) VALIDATION — GREEN reproduces the published model
// ══════════════════════════════════════════════════════════════════

#[test]
fn penner_drive_first_bounce_heights_on_green() {
    // Penner (2002) §3.1: drives land at 26.8 m/s, 68.0° from vertical with
    // 224 rad/s backspin and at 26.2 m/s, 46.2° with 315 rad/s; the first
    // bounce rises 2.82 m and 1.17 m respectively, and the ball leaves with
    // topspin.
    for (speed, angle, w, height) in [(26.8, 68.0, 224.0, 2.82), (26.2, 46.2, 315.0, 1.17)] {
        let r = hit(penner_impact(speed, angle), backspin(w), &Surface::GREEN);
        let peak = r.vel.y * r.vel.y / (2.0 * G);
        assert!(
            (peak - height).abs() < 0.03,
            "{speed} m/s at {angle}°: first bounce peak {peak:.3} m, Penner {height} m"
        );
        assert!(
            r.vel.x > 0.0 && r.spin.z < 0.0,
            "leaves forward with topspin"
        );
    }
}

#[test]
fn penner_nine_iron_bounces_backward_on_green() {
    // Penner (2002) §3.2: a 9-iron landing at 23.8 m/s, 34.0° from vertical
    // with 809 rad/s backspin bounces backward once the rebound is rotated
    // back to the ground frame.
    let r = hit(penner_impact(23.8, 34.0), backspin(809.0), &Surface::GREEN);
    assert!(r.vel.x < 0.0, "bounces backward: vx={:.3}", r.vel.x);
    assert!(r.vel.y * r.vel.y / (2.0 * G) >= 0.005);
}

/// Net run (m) along x from the landing point: hops, slide and roll.
fn run_out(vel: DVec3, w: f64, surface: &Surface) -> f64 {
    let g = surface.settle(&state(vel, backspin(w)), &BALL);
    assert_ends_by_hop_height(&g);
    g.rest.pos.x
}

#[test]
fn penner_drive_runs_match_fig4() {
    // Penner (2002) Fig. 4: the two drives stop about 37 m and 19 m past
    // the landing point on his green.
    let long = run_out(penner_impact(26.8, 68.0), 224.0, &Surface::GREEN);
    let short = run_out(penner_impact(26.2, 46.2), 315.0, &Surface::GREEN);
    assert!((long - 37.0).abs() < 4.0, "68° drive runs {long:.1} m");
    assert!((short - 19.0).abs() < 3.0, "46° drive runs {short:.1} m");
}

#[test]
fn penner_nine_iron_spin_back_grows_with_backspin() {
    // Penner (2002) Fig. 7: the 9-iron landing on his green with 809, 1000
    // and 1200 rad/s backspin. More backspin carries the ball further back.
    let runs: Vec<f64> = [809.0, 1000.0, 1200.0]
        .iter()
        .map(|&w| run_out(penner_impact(23.8, 34.0), w, &Surface::GREEN))
        .collect();
    assert!(
        runs[0] <= 0.0 && runs[1] < runs[0] - 1.0 && runs[2] < runs[1] - 1.0,
        "backward run should grow with backspin: {runs:.2?} m"
    );
}

#[test]
fn penner_nine_iron_firm_green_bounces_forward_then_spins_back() {
    // Penner (2002) Fig. 8: with μ = 0.25 the friction cannot check the
    // backspin on the first impact, so the ball bounces forward still
    // carrying backspin, then runs back; with 1000 and 1200 rad/s it
    // finishes behind its landing point.
    let firm = Surface {
        mu_slide: 0.25,
        ..Surface::GREEN
    };
    for w in [809.0, 1000.0, 1200.0] {
        let first = hit(penner_impact(23.8, 34.0), backspin(w), &firm);
        assert!(
            first.vel.x > 0.0 && first.spin.z > 0.0,
            "{w} rad/s: first bounce forward with backspin"
        );
    }
    let low = run_out(penner_impact(23.8, 34.0), 1000.0, &firm);
    let high = run_out(penner_impact(23.8, 34.0), 1200.0, &firm);
    assert!(
        low < 0.0 && high < low,
        "spin-back: {low:.2} m, {high:.2} m"
    );
}

#[test]
fn ball_rolling_backward_keeps_rolling_backward() {
    // A ball moving back toward the tee with matching rolling spin grips
    // and leaves still moving backward, spinning in the backward-rolling
    // sense (+Z for travel along −X).
    let r = hit(
        DVec3::new(-5.0, -3.0, 0.0),
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
        &state(DVec3::new(20.0, -15.0, 0.0), backspin(rpm_to_rads(3000.0))),
        &BALL,
    );
    assert_ends_by_hop_height(&g);
}

#[test]
fn bounce_sequence_stops_below_threshold() {
    // Straight drop at 15 m/s: no tilt, so each rebound is e(v)·v with
    // Penner's e rising toward 0.51 as the ball slows. Rebounds of 2.26,
    // 0.97, 0.46 and 0.23 m/s rise 260, 48, 10.8 and 2.6 mm: the fourth is
    // below the 5 mm transition height.
    let g = Surface::FAIRWAY.settle(&state(DVec3::new(0.0, -15.0, 0.0), DVec3::ZERO), &BALL);
    assert_eq!(g.contacts.len(), 4);
    assert_ends_by_hop_height(&g);
    assert!(g.rest.pos.length() < 1e-9);
}

#[test]
fn bounce_sequence_safety_bound() {
    // Elastic, frictionless, rigid: every rebound returns at full speed,
    // so only the impact limit ends the sequence (each ~1 s hop would hit
    // the default 30 s time limit first, so it is raised here).
    let physics = GroundPhysics {
        max_time: 1000.0,
        ..GroundPhysics::new(plane(1.0, 0.0, 0.0))
    };
    let start = state(DVec3::new(5.0, -5.0, 0.0), DVec3::ZERO);
    let g = physics.settle(&start, &BALL);
    assert_eq!(g.termination, Termination::MaxImpacts);
    assert_eq!(g.contacts.len(), 200);
    assert!(g.contacts.iter().all(|c| c.vel_out.y > 4.99));
    // The path ends with a point at the last contact, like every other exit.
    let end = g.points.last().unwrap();
    assert_eq!(*end, GroundPoint::new(&g.rest, Phase::Hop));
    assert_eq!(end.time, g.contacts[199].time);
    // With the default limits the same ball runs out of time instead.
    let timed = plane(1.0, 0.0, 0.0).settle(&start, &BALL);
    assert_eq!(timed.termination, Termination::MaxTime);
    assert!((timed.rest.time - 30.0).abs() < 1e-9);
}

#[test]
fn bounce_sequence_height_decreases() {
    let g = Surface::FAIRWAY.settle(
        &state(DVec3::new(15.0, -10.0, 0.0), backspin(rpm_to_rads(2000.0))),
        &BALL,
    );
    for w in g.contacts.windows(2) {
        assert!(w[1].vel_out.y <= w[0].vel_out.y + 0.01);
    }
}

#[test]
fn contacts_record_true_arrival_velocity() {
    // Each contact's vel_in is the velocity arriving at the ground: on flat
    // ground, the previous rebound with its vertical component reversed.
    let g = Surface::GREEN.settle(
        &state(DVec3::new(12.0, -9.0, 1.0), DVec3::new(0.0, -40.0, 400.0)),
        &BALL,
    );
    assert!(g.contacts.len() > 3);
    for w in g.contacts.windows(2) {
        let expected = DVec3::new(w[0].vel_out.x, -w[0].vel_out.y, w[0].vel_out.z);
        assert!((w[1].vel_in - expected).length() < 1e-9);
        assert!(w[1].vel_in.y < 0.0);
        assert_eq!(w[1].spin_in, w[0].spin_out);
        assert!(w[1].time > w[0].time);
    }
}

#[test]
fn soft_drop_ends_by_hop_height() {
    let g = Surface::GREEN.settle(
        &state(DVec3::new(5.0, -1.0, 0.0), backspin(rpm_to_rads(500.0))),
        &BALL,
    );
    assert_ends_by_hop_height(&g);
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
        (DVec3::new(20.0, -12.0, 3.0), 250.0),
        (DVec3::new(15.0, -14.0, -2.0), 300.0),
        (DVec3::new(30.0, -10.0, 4.0), 150.0),
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
        assert!(slip(out.vel, out.spin, DVec3::Y, R).length() < 1e-12);
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
            DVec3::new(15.0, -10.0, 0.0),
            DVec3::new(0.0, w_y, 0.0),
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
        DVec3::new(0.0, 120.0, -4.0 / R),
        DVec3::new(0.0, -120.0, 400.0),
    ] {
        let g = roll(DVec3::new(4.0, 0.0, 0.0), spin, &Surface::FAIRWAY);
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
    let vel = DVec3::new(15.0, -10.0, 0.0);
    let out = hit(vel, DVec3::new(0.0, -200.0, 0.0), &surface);
    let j_n = 1.4 * 10.0;
    let expected = -200.0 + surface.pivot * j_n / (0.4 * R * R);
    assert!(
        (out.spin.y - expected).abs() < 1e-9,
        "ω_y {} vs {expected}",
        out.spin.y
    );
    // Small spins stop rather than reverse.
    let out = hit(vel, DVec3::new(0.0, 5.0, 0.0), &surface);
    assert!(out.spin.y == 0.0);

    // On tilted (compliant) impacts it still acts about the true normal ŷ:
    // the velocity is untouched and only ω_y changes, so the spin along the
    // heading that drives a curving ball's kick is kept.
    for spin in [
        DVec3::new(30.0, -60.0, 300.0),
        DVec3::new(-80.0, 400.0, 900.0),
    ] {
        let vel = DVec3::new(20.0, -14.0, 3.0);
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
    let vel = DVec3::new(18.0, -15.0, 0.0);
    let fade = hit(vel, DVec3::new(0.0, -60.0, 0.0), &Surface::GREEN);
    let draw = hit(vel, DVec3::new(0.0, 60.0, 0.0), &Surface::GREEN);
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
    let k = BALL.inertia_factor;
    let e = 0.4;
    let vel = DVec3::new(18.0, -15.0, 0.0);
    let incidence = vel.x.atan2(-vel.y);
    let beta =
        units::deg_to_rad(15.4) * (vel.length() / 18.6) * (incidence / units::deg_to_rad(44.4));
    let n_eff = DVec3::new(-beta.sin(), beta.cos(), 0.0);
    let v_n = vel.dot(n_eff);
    for w in [-60.0, 60.0, -150.0] {
        let lateral_slip = -R * w * beta.sin();
        let grip = hit(vel, DVec3::new(0.0, w, 0.0), &plane(e, 10.0, 1.0));
        let expected = -k / (1.0 + k) * lateral_slip;
        assert!(
            (grip.vel.z - expected).abs() < 1e-12,
            "grip, ω_y {w}: v_z {:.9}, expected {expected:.9}",
            grip.vel.z
        );

        let mu = 0.05;
        let sliding = hit(vel, DVec3::new(0.0, w, 0.0), &plane(e, mu, 1.0));
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
            DVec3::new(4.0, 0.0, 0.0),
            DVec3::new(w_x, 0.0, -4.0 / R),
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
    let vel = DVec3::new(3.0, 0.0, 0.5);
    let spin = DVec3::new(-20.0, 30.0, 300.0);
    let u0 = slip(vel, spin, DVec3::Y, R).length();
    let expected = u0 / (3.5 * 0.40 * G);
    let g = roll(vel, spin, &Surface::GREEN);
    let roll_start = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
    assert!(
        (roll_start.time - expected).abs() < 1e-12,
        "roll starts at {:.9} s, expected {expected:.9} s",
        roll_start.time
    );
    assert!(slip(roll_start.vel, roll_start.spin, DVec3::Y, R).length() < 1e-12);
    // Sliding friction acts along the initial slip direction throughout:
    // velocity change is exactly −μ g t û.
    let u_hat = slip(vel, spin, DVec3::Y, R).normalize();
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
    // a = μ_roll(1 + (v/v_ref)²)·g/(1+k) the distance is
    // (1+k)·v_ref²/(2·μ_roll·g) · ln(1 + v₀²/v_ref²).
    let v0 = 1.83_f64;
    for surface in [Surface::GREEN, Surface::FAIRWAY] {
        let g = roll(DVec3::new(v0, 0.0, 0.0), rolling(v0), &surface);
        let vr = surface.roll_v_ref;
        let expected =
            1.4 * vr * vr / (2.0 * surface.mu_roll * G) * (1.0 + v0 * v0 / (vr * vr)).ln();
        assert!(
            (g.rest.pos.x - expected).abs() < 1e-5,
            "roll {:.6} m, closed form {expected:.6} m",
            g.rest.pos.x
        );
    }
    // The green preset is a stimp ≈ 11.4 ft green.
    let green = roll(DVec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::GREEN);
    let feet = units::meters_to_feet(green.rest.pos.x);
    assert!((feet - 11.4).abs() < 0.1, "stimp {feet:.2} ft");
}

#[test]
fn stimpmeter_fairway_shorter() {
    let v0 = 1.83;
    let green = roll(DVec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::GREEN);
    let fairway = roll(DVec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::FAIRWAY);
    assert!(fairway.rest.pos.x < green.rest.pos.x);
}

#[test]
fn rollout_stops_and_moves_forward() {
    let g = roll(DVec3::new(10.0, 0.0, 0.0), rolling(10.0), &Surface::FAIRWAY);
    assert_eq!(g.termination, Termination::Rest);
    assert_eq!(g.rest.vel, DVec3::ZERO);
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
    let g = roll(DVec3::new(10.0, 0.0, 2.0), DVec3::ZERO, &Surface::GREEN);
    assert!(g.points.iter().all(|p| p.pos.y == 0.0 && p.vel.y == 0.0));
}

#[test]
fn speed_monotonically_decreases_pure_roll() {
    let g = roll(DVec3::new(8.0, 0.0, 0.0), rolling(8.0), &Surface::GREEN);
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
        (DVec3::new(21.3, -18.5, 0.0), backspin(212.0)),
        (DVec3::new(21.0, -18.0, 2.5), DVec3::new(0.0, -50.0, 210.0)),
        (DVec3::new(9.0, -20.0, -1.0), DVec3::new(10.0, 40.0, 700.0)),
        (DVec3::new(3.0, 0.0, 0.5), DVec3::new(-20.0, 30.0, 300.0)),
    ];
    for (vel, spin) in starts {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let g = surface.settle(&state(vel, spin), &BALL);
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
    let fairway = roll(DVec3::new(5.0, 0.0, 0.0), rolling(5.0), &Surface::FAIRWAY);
    let green = roll(DVec3::new(5.0, 0.0, 0.0), rolling(5.0), &Surface::GREEN);
    assert!(green.rest.pos.x > fairway.rest.pos.x);
}

#[test]
fn backspin_reduces_rollout() {
    let no_spin = roll(DVec3::new(5.0, 0.0, 0.0), rolling(5.0), &Surface::GREEN);
    let with_spin = roll(
        DVec3::new(5.0, 0.0, 0.0),
        backspin(rpm_to_rads(5000.0)),
        &Surface::GREEN,
    );
    assert!(with_spin.rest.pos.x < no_spin.rest.pos.x);
}

#[test]
fn extreme_backspin_can_reverse() {
    let g = roll(
        DVec3::new(2.0, 0.0, 0.0),
        backspin(rpm_to_rads(8000.0)),
        &Surface::GREEN,
    );
    assert!(g.rest.pos.x < 0.0, "should reverse: x={:.3}", g.rest.pos.x);
}

#[test]
fn driver_rollout_reasonable() {
    // ~10 m/s forward, nearly rolling: the speed-dependent rolling
    // resistance (Kensrud et al. 2024) keeps a firm-fairway roll short.
    let g = roll(DVec3::new(10.0, 0.0, 0.0), rolling(8.0), &Surface::FAIRWAY);
    let yards = units::meters_to_yards(g.rest.pos.x);
    assert!(yards > 5.0 && yards < 60.0, "rollout {yards:.1} yds");
}

#[test]
fn ground_physics_step() {
    let physics = GroundPhysics::new(Surface::GREEN);
    let sliding = state(DVec3::new(3.0, 0.0, 0.0), backspin(300.0));
    let (next, phase) = physics.step(&sliding, &BALL);
    assert_eq!(phase, Phase::Slide);
    assert!((next.time - 0.001).abs() < 1e-15 && next.vel.x < 3.0);

    let rolling_ball = state(DVec3::new(3.0, 0.0, 0.0), rolling(3.0));
    let (next, phase) = physics.step(&rolling_ball, &BALL);
    assert_eq!(phase, Phase::Roll);
    assert!((next.spin.z + next.vel.x / R).abs() < 1e-9);

    // A ball slower than one step's deceleration stops inside the step.
    let creeping = state(DVec3::new(1e-4, 0.0, 0.0), rolling(1e-4));
    let (next, phase) = physics.step(&creeping, &BALL);
    assert_eq!(phase, Phase::Roll);
    assert_eq!(next.vel, DVec3::ZERO);
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
    let slow = state(DVec3::new(0.3, 0.0, 0.0), rolling(0.3));
    assert_ne!(physics.step(&slow, &BALL).0.vel, DVec3::ZERO);
    assert_eq!(tuned.step(&slow, &BALL).0.vel, DVec3::ZERO);

    // It samples its own terrain under the ball.
    let split = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| {
        if x < 1.0 {
            Surface::GREEN.sample(x, 0.0)
        } else {
            Surface::FAIRWAY.sample(x, 0.0)
        }
    }));
    let far = BallState {
        pos: DVec3::new(5.0, 0.0, 0.0),
        ..rolling_ball
    };
    let on_fairway = GroundPhysics::new(Surface::FAIRWAY).step(&far, &BALL);
    assert_eq!(split.step(&far, &BALL), on_fairway);
    assert_eq!(
        split.step(&rolling_ball, &BALL),
        physics.step(&rolling_ball, &BALL)
    );
}

#[test]
fn hollow_ball_follows_generalised_formulas() {
    // A thin-shelled ball, k = 2/3, follows the general-k forms of the
    // grip kick, the slide-to-roll time and the stimpmeter distance.
    let k = 2.0 / 3.0;
    let hollow = BallModel {
        inertia_factor: k,
        ..BallModel::TOUR
    };

    // Grip on rigid ground: heading atan(v_z / (v_x − k·R·ω_b)).
    let rigid = plane(0.5, 1.0, 0.0);
    for (vel, w) in [
        (DVec3::new(20.0, -12.0, 3.0), 250.0),
        (DVec3::new(15.0, -14.0, -2.0), 300.0),
    ] {
        let out = impact(&state(vel, backspin(w)), &rigid.sample(0.0, 0.0), &hollow);
        let expected = vel.z.atan2(vel.x - k * R * w);
        assert!((out.vel.z.atan2(out.vel.x) - expected).abs() < 1e-12);
        assert!(slip(out.vel, out.spin, DVec3::Y, R).length() < 1e-12);
    }

    // Slide to roll after |u₀| / (((1+k)/k)·μ·g).
    let vel = DVec3::new(3.0, 0.0, 0.5);
    let spin = DVec3::new(-20.0, 30.0, 300.0);
    let expected = slip(vel, spin, DVec3::Y, R).length() / ((1.0 + k) / k * 0.40 * G);
    let g = Surface::GREEN.settle(&state(vel, spin), &hollow);
    let roll_start = g.points.iter().find(|p| p.phase == Phase::Roll).unwrap();
    assert!((roll_start.time - expected).abs() < 1e-12);

    // Stimpmeter: (1+k)·v_ref²/(2·μ_roll·g)·ln(1 + v₀²/v_ref²).
    let v0 = 1.83_f64;
    for surface in [Surface::GREEN, Surface::FAIRWAY] {
        let g = surface.settle(&state(DVec3::new(v0, 0.0, 0.0), rolling(v0)), &hollow);
        let vr = surface.roll_v_ref;
        let expected =
            (1.0 + k) * vr * vr / (2.0 * surface.mu_roll * G) * (1.0 + v0 * v0 / (vr * vr)).ln();
        assert!(
            (g.rest.pos.x - expected).abs() < 1e-5,
            "roll {:.6} m, closed form {expected:.6} m",
            g.rest.pos.x
        );
    }
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
            DVec3::new(3.0, 0.0, 0.5),
            DVec3::new(-20.0, 30.0, 300.0),
            &Surface::GREEN,
        ),
        roll(
            DVec3::new(2.0, 0.0, 0.0),
            backspin(rpm_to_rads(8000.0)),
            &Surface::GREEN,
        ),
        GroundPhysics {
            max_time: 0.0237,
            ..GroundPhysics::new(Surface::GREEN)
        }
        .settle(&state(DVec3::new(3.0, 0.0, 0.0), backspin(300.0)), &BALL),
        GroundPhysics {
            max_impacts: 3,
            ..GroundPhysics::new(Surface::FAIRWAY)
        }
        .settle(&state(DVec3::new(21.3, -18.5, 0.0), backspin(212.0)), &BALL),
    ];
    for (vel, spin) in [
        (DVec3::new(21.3, -18.5, 0.0), backspin(212.0)),
        (DVec3::new(21.0, -18.0, 2.5), DVec3::new(0.0, -50.0, 210.0)),
        (DVec3::new(9.0, -20.0, -1.0), DVec3::new(10.0, 40.0, 700.0)),
    ] {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            results.push(surface.settle(&state(vel, spin), &BALL));
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
    // all-green distances.
    let v0 = 6.0;
    let mixed = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| GroundSample {
        height: 0.0,
        normal: DVec3::Y,
        surface: if x < 2.0 {
            Surface::FAIRWAY
        } else {
            Surface::GREEN
        },
    }));
    let start = state(DVec3::new(v0, 0.0, 0.0), rolling(v0));
    let g = mixed.settle(&start, &BALL);
    let fairway = roll(DVec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::FAIRWAY);
    let green = roll(DVec3::new(v0, 0.0, 0.0), rolling(v0), &Surface::GREEN);
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
        normal: DVec3::Y,
        surface: if x > 0.5 {
            Surface::GREEN
        } else {
            Surface::FAIRWAY
        },
    });
    let g = GroundPhysics::new(terrain)
        .settle(&state(DVec3::new(10.0, -8.0, 0.0), backspin(300.0)), &BALL);
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
            ..Surface::GREEN
        })
    };
    let g = physics.settle(&state(DVec3::new(3.0, 0.0, 0.0), rolling(3.0)), &BALL);
    assert_eq!(g.termination, Termination::MaxTime);
    assert!((g.rest.time - 2.0).abs() < 1e-9);
    assert!((g.rest.pos.x - 6.0).abs() < 1e-6);
    assert!(g.rest.vel.x > 2.99);
}

#[test]
fn ground_physics_tunables() {
    let start = state(DVec3::new(15.0, -10.0, 0.0), backspin(250.0));
    let default = GroundPhysics::new(Surface::FAIRWAY).settle(&start, &BALL);
    assert_eq!(default, Surface::FAIRWAY.settle(&start, &BALL));
    // A higher hop transition ends the bouncing sooner.
    let fewer = GroundPhysics {
        min_hop_height: 0.05,
        ..GroundPhysics::new(Surface::FAIRWAY)
    }
    .settle(&start, &BALL);
    assert!(fewer.contacts.len() < default.contacts.len());
    // A coarser step barely moves the rest position.
    let coarse = GroundPhysics {
        dt: 0.004,
        ..GroundPhysics::new(Surface::FAIRWAY)
    }
    .settle(&start, &BALL);
    assert!((coarse.rest.pos - default.rest.pos).length() < 1e-3);
}

#[test]
fn default_landing_interpolates_crossing() {
    let a = BallState {
        time: 1.0,
        pos: DVec3::new(10.0, 0.3, 1.0),
        vel: DVec3::new(10.0, -3.0, 1.0),
        spin: backspin(100.0),
    };
    let b = BallState {
        time: 1.1,
        pos: DVec3::new(11.0, -0.1, 1.1),
        vel: DVec3::new(10.0, -4.0, 1.0),
        spin: backspin(90.0),
    };
    let l = Surface::GREEN.landing(&a, &b).unwrap();
    assert!((l.time - 1.075).abs() < 1e-12);
    assert!((l.pos - DVec3::new(10.75, 0.0, 1.075)).length() < 1e-12 && l.pos.y == 0.0);
    assert!((l.vel.y + 3.75).abs() < 1e-12 && (l.spin.z - 92.5).abs() < 1e-12);
    // Leaving the ground is not a landing.
    assert!(Surface::GREEN.landing(&b, &a).is_none());
    let tee = BallState {
        pos: DVec3::ZERO,
        ..a
    };
    assert!(Surface::GREEN.landing(&tee, &a).is_none());
}
