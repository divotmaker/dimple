use super::*;
use crate::aero::Tour;
use crate::ground::{BallState, GroundModel, Surface};

const BALL: Tour = Tour::SEA_LEVEL;
const CUP: Cup = Cup::new(0.0, 0.0);
const R: f64 = REGULATION_RADIUS_M;

/// Hogan & Antali's (2025) values: ball 0.021 m, hole 0.054 m, solid
/// sphere, g 9.81.
const HOGAN: Rim = Rim {
    cup: 0.054,
    ball: 0.021,
    inertia: 0.4,
    gravity: 9.81,
};

/// Mahoney & Connaughton's (2017) values: ball 1.68 in, hole 4.25 in, g 9.8.
const MAHONEY: Rim = Rim {
    cup: 4.25 * 0.0254 / 2.0,
    ball: 0.84 * 0.0254,
    inertia: 0.4,
    gravity: 9.8,
};

fn outcome(speed: f64, offset_frac: f64) -> Outcome {
    CUP.outcome(speed, offset_frac * R, &BALL)
}

/// The highest speed at which a ballistic ball `d` off center drops in.
fn far_rim_limit(rim: &Rim, d: f64) -> f64 {
    let (mut lo, mut hi) = (rim.c1(d), 6.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if rim.far_rim_holds(mid, d) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

#[test]
fn the_regulation_hole_is_four_and_a_quarter_inches_across() {
    assert!((2.0 * REGULATION_RADIUS_M - 4.25 * 0.0254).abs() < 1e-12);
    assert!((Cup::new(1.0, 2.0).radius_m - REGULATION_RADIUS_M).abs() < 1e-15);
    assert!((Cup::new(1.0, 2.0).with_radius(0.1).radius_m - 0.1).abs() < 1e-15);
}

#[test]
fn rim_curves_match_hogan_and_antali() {
    let r = HOGAN.cup;
    // P23: c2 and c3 meet at δ/R 0.940, v 0.316.
    assert!((HOGAN.c2(0.94 * r) - 0.316).abs() < 0.006);
    assert!((HOGAN.c3(0.94 * r) - 0.316).abs() < 0.008);
    // P13: c3 ends on c1 at δ/R 0.748, v 0.689.
    assert!((HOGAN.c1(0.748 * r) - 0.689).abs() < 0.006);
    assert!((HOGAN.c3(0.748 * r) - 0.689).abs() < 0.02);
    // Fig. 5: a rim lip-out and a ball that pitches in.
    assert_eq!(HOGAN.outcome(0.46, 0.93 * r), Outcome::LipOut);
    assert_eq!(HOGAN.outcome(0.575, 0.74 * r), Outcome::Holed);
    // Unbounded where region 4 has no upper edge, and inside c2's asymptote.
    assert!(HOGAN.c3(0.6 * r).is_infinite());
    assert!(HOGAN.c2(0.4 * r).is_infinite());
}

#[test]
fn centered_far_rim_limit_matches_mahoney_and_connaughton() {
    // Their Case A with no bounce and no rim friction: θ 18.41°, V 1.607 m/s.
    // Holmes (1991) gives 1.626 m/s, consistent with his experiments.
    let v = far_rim_limit(&MAHONEY, 0.0);
    assert!((v - 1.607).abs() < 0.001, "{v}");
    // The ball dropping its own radius across the hole is 1.313 m/s, the
    // same paper's speed to reach the far rim at all; faster balls still drop.
    assert!(MAHONEY.far_rim_holds(1.4, 0.0));
}

#[test]
fn far_rim_limit_falls_off_center() {
    let rim = Rim {
        cup: R,
        ball: BALL.radius_m(),
        inertia: 0.4,
        gravity: crate::units::G,
    };
    let mut last = f64::INFINITY;
    for (frac, want) in [(0.0, 1.607), (0.25, 1.555), (0.5, 1.385), (0.7, 1.131)] {
        let v = far_rim_limit(&rim, frac * R);
        assert!((v - want).abs() < 0.002, "{frac}: {v}");
        assert!(v < last);
        last = v;
    }
}

#[test]
fn centered_putts() {
    assert_eq!(outcome(0.2, 0.0), Outcome::Holed);
    assert_eq!(outcome(1.0, 0.0), Outcome::Holed);
    assert_eq!(outcome(1.60, 0.0), Outcome::Holed);
    assert_eq!(outcome(1.62, 0.0), Outcome::Over);
    assert_eq!(outcome(1.60, -0.0), Outcome::Holed);
}

#[test]
fn edge_putts() {
    // Slow balls near the edge fall off the rim into the hole.
    assert_eq!(outcome(0.25, 0.98), Outcome::Holed);
    assert_eq!(outcome(0.25, -0.98), Outcome::Holed);
    // At 0.8 R: pitch in below c3 (0.561), lip out up to c1 (0.762), fly
    // across and drop up to the far-rim limit (0.939), over above it.
    assert_eq!(outcome(0.5, 0.8), Outcome::Holed);
    assert_eq!(outcome(0.65, 0.8), Outcome::LipOut);
    assert_eq!(outcome(0.85, 0.8), Outcome::Holed);
    assert_eq!(outcome(1.0, 0.8), Outcome::Over);
    // A path that misses the hole.
    assert_eq!(outcome(0.1, 1.0), Outcome::Over);
    assert_eq!(outcome(0.1, 1.5), Outcome::Over);
}

#[test]
fn fast_balls_skip_over() {
    assert_eq!(outcome(2.0, 0.0), Outcome::Over);
    assert_eq!(outcome(1.5, 0.5), Outcome::Over);
    assert_eq!(outcome(3.0, 0.9), Outcome::Over);
}

fn point(time: f64, x: f64, z: f64, speed: f64, phase: Phase) -> GroundPoint {
    let state = BallState {
        time,
        pos: Vec3::new(x, 0.0, z),
        vel: Vec3::new(speed, 0.0, 0.0),
        spin: Vec3::ZERO,
    };
    GroundPoint::new(&state, phase)
}

#[test]
fn capture_is_timed_at_the_edge() {
    let points = [
        point(0.0, -0.1, 0.0, 1.0, Phase::Roll),
        point(0.2, 0.1, 0.0, 1.0, Phase::Roll),
    ];
    let got = CUP.capture(&points, &BALL).expect("holed");
    assert!((got.time - (0.1 - R)).abs() < 1e-12);
    assert!((got.pos.x + R).abs() < 1e-12);
    assert!((got.speed_mps - 1.0).abs() < 1e-12);
    assert!(got.offset_m.abs() < 1e-12);

    // Off center: the offset is the segment line's distance from the center.
    let shifted = Cup::new(0.0, 0.02);
    let got = shifted.capture(&points, &BALL).expect("holed");
    assert!((got.offset_m - 0.02).abs() < 1e-12);
}

#[test]
fn a_hop_never_captures() {
    let points = [
        point(0.0, -0.1, 0.0, 0.5, Phase::Hop),
        point(0.2, 0.1, 0.0, 0.5, Phase::Roll),
    ];
    assert_eq!(CUP.capture(&points, &BALL), None);
}

#[test]
fn a_path_starting_inside_has_no_crossing() {
    let points = [
        point(0.0, 0.0, 0.0, 0.5, Phase::Roll),
        point(0.2, 0.1, 0.0, 0.5, Phase::Roll),
    ];
    assert_eq!(CUP.capture(&points, &BALL), None);
}

/// A ball rolling along +X on `GREEN`, starting 1 m short of a hole
/// `offset` to the right of its line.
fn roll_at_cup(speed: f64, offset: f64) -> Option<Capture> {
    let start = BallState {
        time: 0.0,
        pos: Vec3::new(-1.0, 0.0, 0.0),
        vel: Vec3::new(speed, 0.0, 0.0),
        spin: Vec3::new(0.0, 0.0, -speed / BALL.radius_m()),
    };
    let ground = Surface::GREEN.settle(&start, &BALL);
    Cup::new(0.0, offset).capture(&ground.points, &BALL)
}

#[test]
fn a_rolling_putt_drops_or_runs_over() {
    let holed = roll_at_cup(1.5, 0.0).expect("holed");
    assert!(holed.time > 0.0 && holed.speed_mps < 1.5);
    assert!((holed.pos.x + R).abs() < 0.01);
    assert_eq!(roll_at_cup(3.0, 0.0), None);
    // A ball that dies at the hole's edge drops.
    assert!(roll_at_cup(1.12, 0.0).is_some());
}

/// Horizontal speed where a ball rolling along +X from `start` m/s, 1 m
/// short, crosses the hole's edge on a line `offset` from the center.
fn speed_at_edge(start: f64, offset: f64) -> f64 {
    let state = BallState {
        time: 0.0,
        pos: Vec3::new(-1.0, 0.0, 0.0),
        vel: Vec3::new(start, 0.0, 0.0),
        spin: Vec3::new(0.0, 0.0, -start / BALL.radius_m()),
    };
    let edge = -(R * R - offset * offset).sqrt();
    let points = Surface::GREEN.settle(&state, &BALL).points;
    let pair = points
        .windows(2)
        .find(|p| p[0].pos.x <= edge && p[1].pos.x > edge)
        .expect("reaches the hole");
    let t = (edge - pair[0].pos.x) / (pair[1].pos.x - pair[0].pos.x);
    pair[0].vel.x + (pair[1].vel.x - pair[0].vel.x) * t
}

#[test]
fn a_rolling_putt_lips_out_at_speed() {
    // At 0.8 R the ball pitches in when slow, lips out a little faster,
    // flies across and drops faster still, and runs over above that.
    let edge = 0.8 * R;
    for (start, want) in [
        (SLOW, Outcome::Holed),
        (LIP, Outcome::LipOut),
        (ACROSS, Outcome::Holed),
        (FAST, Outcome::Over),
    ] {
        let at_edge = speed_at_edge(start, edge);
        assert_eq!(outcome(at_edge, 0.8), want, "{start} → {at_edge}");
        assert_eq!(roll_at_cup(start, edge).is_some(), want == Outcome::Holed);
    }
}

/// Start speeds reaching the 0.8 R edge at 0.44, 0.67, 0.86 and 1.10 m/s.
const SLOW: f64 = 1.05;
const LIP: f64 = 1.18;
const ACROSS: f64 = 1.31;
const FAST: f64 = 1.5;
