//! [`Surface::absorption`]: the presets, the definition, its limits and the
//! sand example of `docs/AERO.md`.

use super::*;
use crate::aero::Tour;
use crate::units;

const AIR: Tour = Tour::SEA_LEVEL;

/// The sand example of `docs/AERO.md`.
const SAND: Surface = Surface {
    restitution: Restitution::Constant(0.0),
    compliance: 0.0,
    absorption: 0.9,
    mu_roll: 0.3,
    ..Surface::FAIRWAY
};

/// A landing at `speed` m/s, `descent_deg` below the horizontal, heading
/// +X, with `rpm` of backspin (negative: topspin) and `side_rpm` of
/// sidespin (positive curves right).
fn landing(speed: f64, descent_deg: f64, rpm: f64, side_rpm: f64) -> BallState {
    let a = units::deg_to_rad(descent_deg);
    BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel: Vec3::new(speed * a.cos(), -speed * a.sin(), 0.0),
        spin: Vec3::new(0.0, -units::rpm_to_rads(side_rpm), units::rpm_to_rads(rpm)),
    }
}

/// Whitlark & Pringle's (2012) fairway sand wedge landing: 49 mph, 51°,
/// 9100 rpm.
fn wedge() -> BallState {
    landing(units::mph_to_mps(49.0), 51.0, 9100.0, 0.0)
}

/// Their fairway 7-iron landing: 55 mph, 46°, 5100 rpm.
fn seven_iron() -> BallState {
    landing(units::mph_to_mps(55.0), 46.0, 5100.0, 0.0)
}

/// The preset driver's landing: 28.5 m/s, 40°, 2018 rpm.
fn driver() -> BallState {
    landing(28.5, 40.0, 2018.0, 0.0)
}

/// A thin shot: 35 m/s, 10°, 2000 rpm.
fn thin() -> BallState {
    landing(35.0, 10.0, 2000.0, 0.0)
}

/// Kinetic energy per unit mass, translation plus spin.
fn energy(vel: Vec3, spin: Vec3) -> f64 {
    let r = AIR.radius_m();
    0.5 * vel.length_squared() + 0.5 * AIR.inertia_factor() * r * r * spin.length_squared()
}

fn hops(g: &GroundResult) -> bool {
    g.points.iter().any(|p| p.phase == Phase::Hop)
}

#[test]
fn presets_absorb_nothing() {
    for s in [Surface::FAIRWAY, Surface::GREEN, Surface::FIRM_GREEN] {
        assert_eq!(s.absorption, 0.0);
    }
}

#[test]
fn absorption_scales_the_rebound_by_the_landing_angle() {
    let at = |s: Surface| s.sample(0.0, 0.0);
    for l in [wedge(), driver(), landing(35.0, 10.0, 2000.0, 1500.0)] {
        let turf = impact(&l, &at(Surface::FAIRWAY), &AIR);
        let sin = -l.vel.dot(Vec3::Y) / l.vel.length();
        for a in [0.25, 0.5, 0.9, 1.0] {
            let bed = impact(
                &l,
                &at(Surface {
                    absorption: a,
                    ..Surface::FAIRWAY
                }),
                &AIR,
            );
            let kept = (1.0 - a).powf(sin);
            assert_eq!(bed.vel, turf.vel * kept, "absorption {a}");
            assert_eq!(bed.spin, turf.spin * kept, "absorption {a}");
            assert_eq!((bed.time, bed.pos), (l.time, l.pos));
        }
    }
}

#[test]
fn full_absorption_stops_the_ball_on_its_mark() {
    let stop = Surface {
        absorption: 1.0,
        ..Surface::FAIRWAY
    };
    for l in [wedge(), seven_iron(), driver(), thin()] {
        let g = stop.settle(&l, &AIR);
        assert_eq!(g.termination, Termination::Rest);
        assert_eq!(g.contacts.len(), 1);
        assert_eq!(g.rest.pos, l.pos);
        assert_eq!(g.rest.vel, Vec3::ZERO);
        assert_eq!(g.rest.spin, Vec3::ZERO);
    }
}

#[test]
fn sand_does_not_bounce() {
    for speed in [8.0, 15.0, 25.0, 35.0, 45.0, 55.0] {
        for descent in [3.0, 10.0, 20.0, 30.0, 45.0, 60.0, 75.0] {
            for (rpm, side) in [(9100.0, 0.0), (3000.0, 800.0), (-1000.0, -500.0)] {
                let g = SAND.settle(&landing(speed, descent, rpm, side), &AIR);
                let case = format!("{speed} m/s, {descent}°, {rpm} rpm, side {side}");
                assert_eq!(g.termination, Termination::Rest, "{case}");
                assert_eq!(g.contacts.len(), 1, "{case}");
                assert!(!hops(&g), "{case}");
                let r = g.rest;
                assert!(
                    [r.pos, r.vel, r.spin]
                        .iter()
                        .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite()),
                    "{case}"
                );
                for c in &g.contacts {
                    assert!(
                        energy(c.vel_out, c.spin_out) <= energy(c.vel_in, c.spin_in),
                        "{case}"
                    );
                }
            }
        }
    }
}

#[test]
fn steeper_landings_keep_less() {
    let runs: Vec<(f64, f64)> = [10.0, 25.0, 40.0, 55.0, 70.0]
        .iter()
        .map(|&descent| {
            let g = SAND.settle(&landing(25.0, descent, 3000.0, 0.0), &AIR);
            (g.contacts[0].vel_out.length(), g.rest.pos.x)
        })
        .collect();
    for w in runs.windows(2) {
        assert!(w[1].0 < w[0].0 && w[1].1 < w[0].1, "{runs:?}");
    }
    assert!(runs[4].1.abs() < 0.1, "{runs:?}");
}

#[test]
fn sand_skips() {
    // The wedge hops 0.297 m and spins back to −0.23 m on FAIRWAY. On sand
    // it does not hop, keeps 3% of its spin and stops 0.22 m past its mark.
    let turf = Surface::FAIRWAY.settle(&wedge(), &AIR);
    assert!(hops(&turf));
    let g = SAND.settle(&wedge(), &AIR);
    assert!(!hops(&g));
    let c = g.contacts[0];
    assert!(c.spin_out.length() < 0.05 * c.spin_in.length());
    assert!(g.rest.pos.x > 0.0 && g.rest.pos.x < 0.5, "{}", g.rest.pos.x);

    // The 7-iron (0.66 m) and the driver landing (1.46 m) leave a track.
    for l in [seven_iron(), driver()] {
        let x = SAND.settle(&l, &AIR).rest.pos.x;
        assert!((0.3..3.0).contains(&x), "{x}");
    }

    // A thin shot keeps 21.5 of its 35 m/s.
    let out = SAND.settle(&thin(), &AIR).contacts[0].vel_out.length();
    assert!(out > 15.0, "{out}");
}
