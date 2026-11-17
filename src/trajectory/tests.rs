use super::*;

// ── Helper ──

use crate::aero::Tour;
use crate::air::{Wind, standard_air_density};

const BALL: Tour = Tour::SEA_LEVEL;

/// Denver, CO: 1,609 m (5,280 ft), 17 °C.
fn denver() -> Tour {
    Tour::new(standard_air_density(1609.0, 17.0), Wind::None)
}

/// Sea level at `temperature_c`, still air.
fn sea_level_at(temperature_c: f64) -> Tour {
    Tour::new(standard_air_density(0.0, temperature_c), Wind::None)
}

/// Distances in yards, for comparison with published figures.
trait Yards {
    fn carry_yards(&self) -> f64;
}

impl Yards for FlightResult {
    fn carry_yards(&self) -> f64 {
        units::meters_to_yards(self.carry_m)
    }
}

impl ShotResult {
    fn total_yards(&self) -> f64 {
        units::meters_to_yards(self.total_m)
    }
}

fn flight(speed_mph: f64, launch_deg: f64, backspin: f64, sidespin: f64) -> FlightResult {
    let input = ShotInput {
        ball_speed_mph: speed_mph,
        launch_angle_deg: launch_deg,
        launch_azimuth_deg: 0.0,
        backspin_rpm: backspin,
        sidespin_rpm: sidespin,
    };
    simulate_flight(&input, &BALL)
}

fn flight_at(input: &ShotInput, model: &Tour) -> FlightResult {
    simulate_flight(input, model)
}

// ══════════════════════════════════════════════════════════════════
//  STRUCTURAL TESTS — invariants that must hold regardless of model
// ══════════════════════════════════════════════════════════════════

#[test]
fn trajectory_starts_at_origin() {
    let r = flight(132.0, 16.0, 7000.0, 0.0);
    assert!(r.points[0].pos.length() < 1e-10);
}

#[test]
fn apex_positive() {
    let r = flight(132.0, 16.0, 7000.0, 0.0);
    assert!(r.apex_m > 0.0);
}

#[test]
fn energy_decreasing() {
    let r = flight(132.0, 16.0, 7000.0, 0.0);
    let energy = |p: &TrajectoryPoint| {
        // Per unit mass.
        let ke = 0.5 * p.vel.length_squared();
        let pe = crate::units::G * p.pos.y;
        ke + pe
    };
    for w in r.points.windows(2) {
        assert!(
            energy(&w[1]) <= energy(&w[0]) + 1e-6,
            "energy increased at t={:.3}",
            w[1].time
        );
    }
}

#[test]
fn interpolate_at_start() {
    let r = flight(132.0, 16.0, 7000.0, 0.0);
    let pos = interpolate_trajectory(&r.points, 0.0).unwrap();
    assert!(pos.length() < 1e-6);
}

#[test]
fn interpolate_midpoint_above_ground() {
    let r = flight(132.0, 16.0, 7000.0, 0.0);
    let pos = interpolate_trajectory(&r.points, r.flight_time / 2.0).unwrap();
    assert!(pos.y > 0.0);
}

#[test]
fn interpolate_past_end_returns_last() {
    let r = flight(132.0, 16.0, 7000.0, 0.0);
    let pos = interpolate_trajectory(&r.points, r.flight_time + 10.0).unwrap();
    let last = r.points.last().unwrap().pos;
    assert!((pos - last).length() < 1e-10);
}

// ══════════════════════════════════════════════════════════════════
//  LATERAL SIGN CONVENTION — positive = right of the target line
// ══════════════════════════════════════════════════════════════════

fn lateral_input(azimuth_deg: f64, sidespin_rpm: f64) -> ShotInput {
    ShotInput {
        ball_speed_mph: 150.0,
        launch_angle_deg: 12.0,
        launch_azimuth_deg: azimuth_deg,
        backspin_rpm: 3000.0,
        sidespin_rpm,
    }
}

#[test]
fn positive_azimuth_starts_right() {
    let r = flight_at(&lateral_input(5.0, 0.0), &BALL);
    assert!(
        r.lateral_m > 0.0,
        "lateral {:.2}m should be > 0",
        r.lateral_m
    );
    assert!(r.points[1].pos.z > 0.0, "ball should leave toward +Z");
    assert!(r.points.iter().all(|p| p.pos.z >= 0.0));
}

#[test]
fn negative_azimuth_starts_left() {
    let r = flight_at(&lateral_input(-5.0, 0.0), &BALL);
    assert!(
        r.lateral_m < 0.0,
        "lateral {:.2}m should be < 0",
        r.lateral_m
    );
    assert!(r.points[1].pos.z < 0.0, "ball should leave toward -Z");
}

#[test]
fn slice_spin_curves_right() {
    // Right-handed slice: square start, positive sidespin.
    let r = flight_at(&lateral_input(0.0, 1500.0), &BALL);
    assert!(
        r.lateral_m > 5.0,
        "slice should finish right: {:.2}m",
        r.lateral_m
    );
    // Curve, not a push: lateral velocity grows toward the right.
    let last = r.points.last().unwrap();
    assert!(last.vel.z > 0.0);
}

#[test]
fn hook_spin_curves_left() {
    let r = flight_at(&lateral_input(0.0, -1500.0), &BALL);
    assert!(
        r.lateral_m < -5.0,
        "hook should finish left: {:.2}m",
        r.lateral_m
    );
}

#[test]
fn push_draw_starts_right_and_curves_left() {
    let r = flight_at(&lateral_input(4.0, -800.0), &BALL);
    assert!(r.points[1].vel.z > 0.0, "should start right");
    assert!(
        r.points.last().unwrap().vel.z < 0.0,
        "should be moving left at landing"
    );
}

#[test]
fn lateral_mirror_symmetry() {
    let right = flight_at(&lateral_input(3.0, 1200.0), &BALL);
    let left = flight_at(&lateral_input(-3.0, -1200.0), &BALL);
    assert!((right.lateral_m + left.lateral_m).abs() < 1e-9);
    assert!((right.carry_m - left.carry_m).abs() < 1e-9);
}

#[test]
fn lateral_matches_landing_position() {
    let r = flight_at(&lateral_input(2.0, 900.0), &BALL);
    assert!((r.lateral_m - r.points.last().unwrap().pos.z).abs() < 1e-12);
}

#[test]
fn landing_sidespin_keeps_input_sign() {
    let input = lateral_input(0.0, 1500.0);
    let r = flight_at(&input, &BALL);
    let initial = units::rpm_to_rads(1500.0);
    // Sidespin is about -Y: positive curves right.
    let sidespin = -r.landing_spin.y;
    assert!(sidespin > 0.0 && sidespin < initial);
    let r = flight_at(&lateral_input(0.0, -1500.0), &BALL);
    assert!(-r.landing_spin.y < 0.0);
}

#[test]
fn slice_shot_finishes_right() {
    let r = shot_on_fairway(&lateral_input(0.0, 1500.0));
    assert!(
        r.ground.rest.pos.z > 0.0,
        "rest z={:.2}",
        r.ground.rest.pos.z
    );
    assert!(r.flight.lateral_m > 0.0);
}

// ══════════════════════════════════════════════════════════════════
//  ALTITUDE — Titleist Learning Lab published data
// ══════════════════════════════════════════════════════════════════

#[test]
fn altitude_denver_driver_gain() {
    // Source: Titleist Learning Lab, "Altitude and Golf Ball Flight"
    // Denver (5,280 ft / 1,609m): ~6% more carry.
    // Formula: elevation_ft × 0.00116 = fractional gain.
    // 5280 × 0.00116 = 6.1%
    let driver = ShotInput::driver();
    let sea = flight_at(&driver, &BALL);
    let denver = flight_at(&driver, &denver());
    let gain_pct = (denver.carry_yards() - sea.carry_yards()) / sea.carry_yards() * 100.0;
    // The tour model gives about 7.5% here (Titleist ~6.1%, which assumes a
    // standard day; Denver at 17 °C is warmer than standard).
    assert!(
        gain_pct > 4.0 && gain_pct < 9.0,
        "Denver gain {gain_pct:.1}%, expected ~6.1% (Titleist: elevation_ft × 0.00116)"
    );
}

#[test]
fn altitude_increases_carry() {
    // Directional: higher altitude = more carry (less drag).
    let driver = ShotInput::driver();
    let sea = flight_at(&driver, &BALL);
    let denver = flight_at(&driver, &denver());
    assert!(denver.carry_yards() > sea.carry_yards());
}

#[test]
fn altitude_reduces_curve() {
    // Directional: higher altitude = less Magnus = less lateral deviation.
    let fade = ShotInput {
        ball_speed_mph: 150.0,
        launch_angle_deg: 12.0,
        launch_azimuth_deg: 0.0,
        backspin_rpm: 3000.0,
        sidespin_rpm: 2000.0,
    };
    let sea = flight_at(&fade, &BALL);
    let denver = flight_at(&fade, &denver());
    assert!(
        denver.lateral_m.abs() < sea.lateral_m.abs(),
        "Denver lateral ({:.1}m) should be less than sea level ({:.1}m)",
        denver.lateral_m,
        sea.lateral_m
    );
}

// ══════════════════════════════════════════════════════════════════
//  TEMPERATURE — Andrew Rice Golf / WeatherWorks published data
// ══════════════════════════════════════════════════════════════════

#[test]
fn temperature_hot_vs_cold() {
    // Source: WeatherWorks (weatherworksinc.com), Andrew Rice Golf
    // 90°F vs 50°F → ~8 yards difference for driver.
    // Note: this only captures air density effect, not ball COR change.
    let driver = ShotInput::driver();
    let hot = sea_level_at(32.2); // 90°F
    let cold = sea_level_at(10.0); // 50°F
    let r_hot = flight_at(&driver, &hot);
    let r_cold = flight_at(&driver, &cold);
    let delta = r_hot.carry_yards() - r_cold.carry_yards();
    // Published: ~8 yds per 40°F, which includes ball COR changes. The model
    // responds to air density only (ball COR is unmodeled) and gives ~9 yds.
    // See docs/AERO.md.
    assert!(
        delta > 0.5 && delta < 16.0,
        "hot-cold delta {delta:.1} yds, expected ~8 yds (WeatherWorks)"
    );
}

#[test]
fn temperature_direction() {
    // Hot air is less dense → more carry.
    let driver = ShotInput::driver();
    let hot = sea_level_at(35.0);
    let cold = sea_level_at(5.0);
    assert!(flight_at(&driver, &hot).carry_yards() > flight_at(&driver, &cold).carry_yards());
}

// ══════════════════════════════════════════════════════════════════
//  SPIN DECAY — tutelman.com / TrajectoWare / Trackman published
// ══════════════════════════════════════════════════════════════════

#[test]
fn spin_decay_driver_6s() {
    // Source: tutelman.com/golf/ballflight/spinDecay.php
    // Trackman (Oct 2010): ~4% per second → at 6s: 0.96^6 = 78.3% remains.
    // TrajectoWare Drive: 3.3%/s → at 6s: 0.967^6 = 81.5% remains.
    // Our model: exp(-0.04 × 6) = 78.7% — matches Trackman consensus.
    let mut ball = BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel: Vec3::new(70.0, 0.0, 0.0),
        spin: ShotInput::driver().launch_spin_rads(),
    };
    let initial = ball.spin.length();
    for _ in 0..6000 {
        ball.spin = BALL.spin_after(&ball, 0.001);
    }
    let remaining_pct = ball.spin.length() / initial * 100.0;
    assert!(
        (remaining_pct - 78.0).abs() < 6.0,
        "after 6s: {remaining_pct:.1}% spin remaining, expected ~78% (Trackman)"
    );
}

// ══════════════════════════════════════════════════════════════════
//  LANDING — interpolated ground crossing
// ══════════════════════════════════════════════════════════════════

#[test]
fn flight_lands_exactly_on_the_ground() {
    for input in [
        ShotInput::driver(),
        ShotInput::seven_iron(),
        ShotInput::wedge(),
    ] {
        let r = flight_at(&input, &BALL);
        let last = r.points.last().unwrap();
        assert!(
            last.pos.y.abs() < f64::EPSILON,
            "landing y = {}",
            last.pos.y
        );
        assert!(last.vel.y < 0.0);
        assert!((r.flight_time - last.time).abs() < f64::EPSILON);
        let before = &r.points[r.points.len() - 2];
        assert!(before.pos.y > 0.0 && last.time - before.time <= 0.010 + 1e-12);
    }
}

#[test]
fn shot_flight_matches_simulate_flight_on_flat_ground() {
    let input = lateral_input(2.0, 700.0);
    let f = flight_at(&input, &BALL);
    let s = shot_on_fairway(&input);
    assert!((f.carry_m - s.flight.carry_m).abs() < 1e-12);
    assert!((f.landing_spin - s.flight.landing_spin).length() < 1e-12);
    let c = &s.ground.contacts[0];
    assert!((c.pos - s.flight.points.last().unwrap().pos).length() < 1e-12);
    assert!((c.time - f.flight_time).abs() < 1e-12);
}

// ══════════════════════════════════════════════════════════════════
//  LAUNCH SPIN AXIS — backspin perpendicular to the start line
// ══════════════════════════════════════════════════════════════════

#[test]
fn launch_spin_axis() {
    let straight = lateral_input(0.0, 600.0).launch_spin_rads();
    let expected = Vec3::new(0.0, -units::rpm_to_rads(600.0), units::rpm_to_rads(3000.0));
    assert!((straight - expected).length() < 1e-12);

    // Backspin axis is d × ŷ for launch heading d.
    let az = units::deg_to_rad(7.0);
    let pushed = lateral_input(7.0, 0.0).launch_spin_rads();
    let d = Vec3::new(az.cos(), 0.0, az.sin());
    let axis = d.cross(Vec3::Y);
    assert!((pushed - units::rpm_to_rads(3000.0) * axis).length() < 1e-9);
}

#[test]
fn pushed_straight_shot_stays_on_its_start_line() {
    // A 6° push with no sidespin: the Magnus force stays in the vertical
    // plane through the start line, so the ball lands on that line.
    let input = lateral_input(6.0, 0.0);
    let az = units::deg_to_rad(6.0);
    let across = Vec3::new(-az.sin(), 0.0, az.cos()); // horizontal, ⟂ start line

    let ball = BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel: 60.0 * Vec3::new(az.cos() * 0.9, 0.3, az.sin() * 0.9).normalize(),
        spin: input.launch_spin_rads(),
    };
    let gravity_free = BALL.acceleration(&ball, 10.0);
    assert!(
        gravity_free.dot(across).abs() < 1e-12,
        "aero force has a sideways component {:.3e}",
        gravity_free.dot(across)
    );

    let r = flight_at(&input, &BALL);
    for p in &r.points {
        assert!(
            p.pos.dot(across).abs() < 1e-9,
            "off the start line by {:.3e} m at t={:.2}",
            p.pos.dot(across),
            p.time
        );
    }
    assert!((r.lateral_m - r.carry_m * az.sin()).abs() < 1e-9);
    assert!(r.landing_spin.y.abs() < 1e-12);
    assert!(r.landing_spin.dot(across) > 0.0);
}

#[test]
fn landing_spin_projections_are_exact() {
    let input = lateral_input(-4.0, 900.0);
    let r = flight_at(&input, &BALL);
    let decay = r.landing_spin.length() / input.launch_spin_rads().length();
    assert!(decay > 0.0 && decay < 1.0);
    // Projected on the launch backspin axis and on -Y.
    let az = units::deg_to_rad(-4.0);
    let backspin_axis = Vec3::new(-az.sin(), 0.0, az.cos());
    assert!(
        (r.landing_spin.dot(backspin_axis) - units::rpm_to_rads(3000.0) * decay).abs() < 1e-9
            && (-r.landing_spin.y - units::rpm_to_rads(900.0) * decay).abs() < 1e-9
    );
}

// ══════════════════════════════════════════════════════════════════
//  UNIFIED SHOT SIMULATION — flight + ground
// ══════════════════════════════════════════════════════════════════

use crate::ground::{GroundPhysics, TEST_TURF, Termination};

fn shot_on_fairway(input: &ShotInput) -> ShotResult {
    simulate_shot(input, &BALL, &Surface::FAIRWAY)
}

fn shot_on_green(input: &ShotInput) -> ShotResult {
    simulate_shot(input, &BALL, &Surface::GREEN)
}

#[test]
fn simulate_shot_accepts_dyn_ground_model() {
    let models: [&dyn GroundModel; 2] = [&Surface::GREEN, &GroundPhysics::new(Surface::GREEN)];
    let a = simulate_shot(&ShotInput::driver(), &BALL, models[0]);
    let b = simulate_shot(&ShotInput::driver(), &BALL, models[1]);
    let c = shot_on_green(&ShotInput::driver());
    assert_eq!(a.ground, c.ground);
    assert_eq!(b.ground, c.ground);
}

#[test]
fn shot_total_greater_than_carry_driver() {
    // Driver on fairway: low spin, should roll forward significantly.
    let r = shot_on_fairway(&ShotInput::driver());
    assert!(
        r.total_yards() > r.flight.carry_yards(),
        "driver total ({:.1}) should exceed carry ({:.1})",
        r.total_yards(),
        r.flight.carry_yards()
    );
}

#[test]
fn shot_driver_total_distance() {
    // Fairway drives landing near 40° run about 21 yd: the mean PGA TOUR
    // roll-out line of 2016–18 gives 21.7 yd at 39° (events 10–60 yd;
    // USGA/R&A 2021). The preset launches slower than the PGA Tour average
    // (Trackman: 282 yd carry), carries about 283 yd, lands at about 40°
    // and runs about 21 yd on FAIRWAY.
    let r = shot_on_fairway(&ShotInput::driver());
    let run = r.total_yards() - r.flight.carry_yards();
    assert!(
        (260.0..290.0).contains(&r.flight.carry_yards()),
        "driver carry should be 260-290 yds, got {:.1}",
        r.flight.carry_yards()
    );
    assert!(
        (14.0..22.0).contains(&run),
        "driver fairway run should be 14-22 yds, got {run:.1}"
    );
    assert!(
        (276.0..310.0).contains(&r.total_yards()),
        "driver total should be 276-310 yds, got {:.1}",
        r.total_yards()
    );
}

#[test]
fn shot_7iron_carry_and_run() {
    // PGA Tour 7-iron (Trackman): about 120–123 mph, 16°, 7100 rpm,
    // 172–176 yd carry. The model carries about 171 yd and runs about
    // 5.8 yd on FAIRWAY.
    let input = ShotInput {
        ball_speed_mph: 123.0,
        launch_angle_deg: 16.3,
        launch_azimuth_deg: 0.0,
        backspin_rpm: 7097.0,
        sidespin_rpm: 0.0,
    };
    let r = shot_on_fairway(&input);
    let run = r.total_yards() - r.flight.carry_yards();
    assert!(
        (165.0..185.0).contains(&r.flight.carry_yards()),
        "7-iron carry should be 165-185 yds, got {:.1}",
        r.flight.carry_yards()
    );
    assert!(
        (3.0..12.0).contains(&run),
        "7-iron fairway run should be 3-12 yds, got {run:.1}"
    );
}

#[test]
fn shot_wedge_minimal_rollout_on_green() {
    // Wedge on green: high spin, should check up with minimal rollout.
    let r = shot_on_green(&ShotInput::wedge());
    let rollout_yards = r.total_yards() - r.flight.carry_yards();
    assert!(
        rollout_yards.abs() < 15.0,
        "wedge on green rollout should be minimal, got {rollout_yards:.1} yds"
    );
}

#[test]
fn shot_driver_first_bounce_height_in_penner_range() {
    // Penner (2002) §3.1 reports first-bounce heights from 2.82 m
    // (landing 68.0° from vertical) to 1.17 m (46.2°) for drives on his
    // green model. The driver preset lands about 50° from vertical.
    for r in [
        shot_on_green(&ShotInput::driver()),
        shot_on_fairway(&ShotInput::driver()),
    ] {
        let v = r.ground.contacts[0].vel_out;
        let peak = v.y * v.y / (2.0 * units::G);
        assert!(
            (1.17..=2.82).contains(&peak),
            "driver first bounce peak {peak:.2} m outside Penner's 1.17-2.82 m"
        );
    }
}

#[test]
fn shot_bounces_end_by_hop_height_rule() {
    // Every preset shot stops bouncing by Penner's 5 mm hop-height rule,
    // well inside the safety bound on impacts, and comes to rest.
    let physics = GroundPhysics::new(Surface::FAIRWAY);
    let hop = |v: Vec3| v.y * v.y / (2.0 * units::G);
    for input in [
        ShotInput::driver(),
        ShotInput::seven_iron(),
        ShotInput::wedge(),
    ] {
        for r in [shot_on_fairway(&input), shot_on_green(&input)] {
            assert_eq!(r.ground.termination, Termination::Rest);
            let contacts = &r.ground.contacts;
            assert!(contacts.len() < physics.max_impacts);
            let (last, rest) = contacts.split_last().unwrap();
            assert!(
                rest.iter()
                    .all(|c| hop(c.vel_out) >= physics.min_hop_height)
            );
            assert!(
                last.vel_out.y.abs() < 1e-12,
                "last of {} contacts leaves along the ground",
                contacts.len()
            );
        }
    }
}

#[test]
fn shot_has_bounces() {
    let r = shot_on_fairway(&ShotInput::driver());
    assert!(
        r.ground.contacts.len() > 1,
        "a drive should bounce more than once"
    );
}

#[test]
fn shot_rest_on_ground() {
    let r = shot_on_fairway(&ShotInput::seven_iron());
    assert!(
        r.ground.rest.pos.y.abs() < 1e-12,
        "final position should be on ground: y={:.6}",
        r.ground.rest.pos.y
    );
    assert_eq!(r.ground.rest.vel, Vec3::ZERO);
}

#[test]
fn shot_rest_downrange() {
    let r = shot_on_fairway(&ShotInput::driver());
    assert!(
        r.ground.rest.pos.x > 0.0,
        "driver should finish downrange: x={:.1}",
        r.ground.rest.pos.x
    );
}

#[test]
fn shot_path_times_are_continuous() {
    // Ground points carry absolute times: flight points then ground
    // points form one time-ordered path.
    let r = shot_on_fairway(&ShotInput::seven_iron());
    let flight_end = r.flight.points.last().unwrap();
    let first = r.ground.points.first().unwrap();
    assert!((first.time - flight_end.time).abs() < 1e-12);
    assert!((first.pos - flight_end.pos).length() < 1e-12);
    for w in r.ground.points.windows(2) {
        assert!(w[1].time >= w[0].time);
    }
    assert!((r.ground.points.last().unwrap().pos - r.ground.rest.pos).length() < 1e-12);
}

#[test]
fn shot_landing_spin_decayed() {
    // Landing spin should be less than initial spin (decay during flight).
    let flight = simulate_flight(&ShotInput::driver(), &BALL);
    let initial_backspin = units::rpm_to_rads(2700.0);
    assert!(
        flight.landing_spin.z < initial_backspin,
        "landing backspin ({:.1}) should be less than initial ({:.1})",
        flight.landing_spin.z,
        initial_backspin
    );
    assert!(
        flight.landing_spin.z > 0.0,
        "should still have some backspin"
    );
}

#[test]
fn shot_fairway_less_rollout_than_green() {
    // Same drive on fairway and green: both presets bounce alike, so the
    // hops are the same, and the fairway's higher rolling resistance
    // gives it the shorter run (about 17 yd against 23).
    let input = ShotInput::driver();
    let fairway = shot_on_fairway(&input);
    let green = shot_on_green(&input);
    assert_eq!(fairway.ground.contacts.len(), green.ground.contacts.len());
    let rollout_fairway = fairway.total_yards() - fairway.flight.carry_yards();
    let rollout_green = green.total_yards() - green.flight.carry_yards();
    assert!(
        rollout_fairway < rollout_green - 3.0,
        "fairway rollout ({rollout_fairway:.1}) should be well short of green ({rollout_green:.1})"
    );
}

// ══════════════════════════════════════════════════════════════════
//  RUN-OUT — the hop train dies out; roll by club
// ══════════════════════════════════════════════════════════════════

/// Trackman's normalization, 25 °C at sea level.
fn trackman_air() -> Tour {
    sea_level_at(25.0)
}

fn launch(mph: f64, launch_deg: f64, rpm: f64) -> ShotInput {
    ShotInput {
        ball_speed_mph: mph,
        launch_angle_deg: launch_deg,
        launch_azimuth_deg: 0.0,
        backspin_rpm: rpm,
        sidespin_rpm: 0.0,
    }
}

#[test]
fn shot_hop_train_dies_out() {
    // After the first bounce each hop rises at most about 0.6× as high as
    // the one before (asserted: 0.7×), so a shot makes a handful of
    // impacts and then rolls. No train of hops a few millimeters high
    // skims on for dozens of impacts: at most two hops below 1 cm, from wedges to drives, on
    // both presets.
    let hop = |v: Vec3| v.y * v.y / (2.0 * units::G);
    for input in [
        launch(65.0, 32.0, 6500.0),  // lob wedge
        launch(52.0, 32.0, 4000.0),  // pitch
        launch(30.0, 18.0, 2500.0),  // low chip
        launch(120.0, 16.3, 7097.0), // 7-iron
        launch(171.0, 10.4, 2545.0), // drive
        ShotInput::driver(),
        ShotInput::seven_iron(),
        ShotInput::wedge(),
    ] {
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let r = simulate_shot(&input, &BALL, &surface);
            let contacts = &r.ground.contacts;
            let label = format!("{input:?} on {surface:?}");
            assert!(contacts.len() <= 10, "{label}: {} impacts", contacts.len());
            let heights: Vec<f64> = contacts.iter().map(|c| hop(c.vel_out)).collect();
            for pair in heights[1..].windows(2) {
                assert!(pair[1] <= 0.7 * pair[0], "{label}: hops {heights:?}");
            }
            let low = heights[..heights.len() - 1]
                .iter()
                .filter(|&&h| h < 0.01)
                .count();
            assert!(low <= 2, "{label}: {low} hops below 1 cm");
        }
    }
}

#[test]
fn shot_roll_by_club_on_fairway() {
    // FAIRWAY roll and total for three launches in Trackman's air.
    // A lob wedge (65 mph, 32°, 6500 rpm) lands at about 44° and rolls
    // about 9.3 yd; a 50-yard pitch (52 mph, 32°, 4000 rpm) lands at
    // about 39° and rolls about 14.8 yd. The PGA Tour average drive
    // (Trackman: 171 mph, 10.4°, 2545 rpm) lands at about 39° and rolls
    // about 23 yd, against 21.7 yd at 39° on the mean PGA TOUR roll-out
    // line of 2016–18 (USGA/R&A 2021). Wedge and pitch landings at 39–44°
    // and 17–19 m/s are inside the measured impacts, but no published
    // roll-out covers them on a fairway, so their windows guard against
    // drift rather than validate.
    for (name, input, roll, total) in [
        (
            "lob wedge",
            launch(65.0, 32.0, 6500.0),
            8.0..10.5,
            76.0..87.0,
        ),
        ("pitch", launch(52.0, 32.0, 4000.0), 13.0..16.5, 60.0..71.0),
        (
            "drive",
            launch(171.0, 10.4, 2545.0),
            19.5..25.5,
            305.0..330.0,
        ),
    ] {
        let r = simulate_shot(&input, &trackman_air(), &Surface::FAIRWAY);
        let run = r.total_yards() - r.flight.carry_yards();
        assert!(roll.contains(&run), "{name}: roll {run:.2} yd");
        assert!(
            total.contains(&r.total_yards()),
            "{name}: total {:.1} yd",
            r.total_yards()
        );
    }
}

#[test]
fn shot_total_is_the_rest_point_and_descent_is_reported() {
    let r = shot_on_fairway(&ShotInput::driver());
    let rest = r.ground.rest.pos;
    assert!((r.total_m - rest.x.hypot(rest.z)).abs() < 1e-12);
    assert!(r.total_m > r.flight.carry_m);
    let last = r.flight.points.last().unwrap().vel;
    let descent = units::rad_to_deg((-last.y).atan2(last.x.hypot(last.z)));
    assert!((r.flight.descent_angle_deg - descent).abs() < 1e-12);
    assert!((30.0..50.0).contains(&r.flight.descent_angle_deg));
}

// ══════════════════════════════════════════════════════════════════
//  GROUND LATERAL SIGNS — fades rest right, draws rest left
// ══════════════════════════════════════════════════════════════════

fn driver_with_sidespin(sidespin_rpm: f64) -> ShotInput {
    ShotInput {
        sidespin_rpm,
        ..ShotInput::driver()
    }
}

#[test]
fn fade_rests_right_draw_rests_left() {
    for (surface, name) in [(Surface::FAIRWAY, "fairway"), (Surface::GREEN, "green")] {
        let fade = simulate_shot(&driver_with_sidespin(600.0), &BALL, &surface);
        let draw = simulate_shot(&driver_with_sidespin(-600.0), &BALL, &surface);
        assert!(
            fade.ground.rest.pos.z > fade.flight.lateral_m && fade.flight.lateral_m > 0.0,
            "{name}: fade lands right and keeps moving right: land {:.2}, rest {:.2}",
            fade.flight.lateral_m,
            fade.ground.rest.pos.z
        );
        assert!(
            draw.ground.rest.pos.z < draw.flight.lateral_m && draw.flight.lateral_m < 0.0,
            "{name}: draw lands left and keeps moving left: land {:.2}, rest {:.2}",
            draw.flight.lateral_m,
            draw.ground.rest.pos.z
        );
    }
}

#[test]
#[allow(clippy::float_cmp, reason = "no sidespin, no lateral motion at all")]
fn sidespin_hop_curves_right() {
    // A hop through air carrying positive sidespin (angular velocity
    // along -Y, a fade's) curves right, like the flight; the mirrored
    // spin curves it left by the same amount.
    let hop = |sidespin: f64| {
        let leave = BallState {
            time: 0.0,
            pos: Vec3::ZERO,
            vel: Vec3::new(6.0, 2.0, 0.0),
            spin: Vec3::new(0.0, -sidespin, 0.0),
        };
        let g = Surface::FAIRWAY.settle(&leave, &BALL);
        g.contacts[1].pos
    };
    let (right, left) = (hop(200.0), hop(-200.0));
    assert!(right.z > 1e-3, "fade spin: z {:.6} m", right.z);
    assert!(
        (left.z + right.z).abs() < 1e-12,
        "mirror: {left:?} vs {right:?}"
    );
    assert_eq!(hop(0.0).z, 0.0);
}

#[test]
fn fade_first_bounce_kicks_right_of_its_heading() {
    // A fade lands heading right of its start line while its backspin
    // axis still points across the start line. Friction removes more
    // velocity along the start line than across it, so the rebound's
    // heading is further right than the arrival heading (a draw: left).
    for (side, sign) in [(600.0, 1.0), (-600.0, -1.0)] {
        let r = shot_on_fairway(&driver_with_sidespin(side));
        let c = &r.ground.contacts[0];
        let heading_in = c.vel_in.z.atan2(c.vel_in.x);
        let heading_out = c.vel_out.z.atan2(c.vel_out.x);
        assert!(
            sign * heading_in > 0.0 && sign * (heading_out - heading_in) > 0.0,
            "sidespin {side}: heading in {:.3}°, out {:.3}°",
            units::rad_to_deg(heading_in),
            units::rad_to_deg(heading_out)
        );
    }
}

#[test]
fn pivot_does_not_steer_on_level_ground() {
    // Pivot friction acts about the true normal, so it never changes an
    // impact's velocity or a slide or roll on level ground. Spin about
    // the normal reaches the path only through a later impact's tilted
    // contact plane and through the Magnus force it gives the hops in
    // between. Every lever that clears that spin at each impact, from the
    // presets' 1 mm up, gives the same path, on rigid ground (no tilt)
    // and on compliant turf.
    let rigid = Surface {
        compliance: 0.0,
        mu_roll: 0.12,
        plow: 0.36,
        ..TEST_TURF
    };
    for side in [600.0, -600.0] {
        let input = driver_with_sidespin(side);
        for (base, levers) in [
            (rigid, [0.001, 0.002, 0.005]),
            (Surface::FAIRWAY, [0.001, 0.002, 0.005]),
            (Surface::GREEN, [0.001, 0.002, 0.005]),
        ] {
            let shots: Vec<ShotResult> = levers
                .iter()
                .map(|&pivot| simulate_shot(&input, &BALL, &Surface { pivot, ..base }))
                .collect();
            for (shot, pivot) in shots.iter().zip(levers).skip(1) {
                assert!(
                    (shot.ground.rest.pos - shots[0].ground.rest.pos).length() < 1e-9,
                    "sidespin {side}, pivot {pivot}: rest {:?} vs {:?}",
                    shot.ground.rest.pos,
                    shots[0].ground.rest.pos
                );
                assert_eq!(shot.ground.contacts.len(), shots[0].ground.contacts.len());
            }
        }
        // With no pivot friction on rigid ground the sidespin survives
        // into the hops, and their Magnus force curves the ball further
        // the way the sidespin curves it in flight.
        let kept = simulate_shot(
            &input,
            &BALL,
            &Surface {
                pivot: 0.0,
                ..rigid
            },
        );
        let cleared = simulate_shot(&input, &BALL, &rigid);
        let curve = side.signum() * (kept.ground.rest.pos.z - cleared.ground.rest.pos.z);
        assert!(
            curve > 0.0 && curve < 0.2,
            "sidespin {side}: curve {curve:.3} m"
        );
        // With no pivot friction on compliant turf the spin survives to
        // the later impacts, whose tilt throws move the rest slightly.
        let free = simulate_shot(
            &input,
            &BALL,
            &Surface {
                pivot: 0.0,
                ..Surface::FAIRWAY
            },
        );
        let preset = shot_on_fairway(&input);
        let shift = (free.ground.rest.pos - preset.ground.rest.pos).length();
        assert!(shift > 0.01 && shift < 0.2, "shift {shift:.3} m");
    }
}

// ══════════════════════════════════════════════════════════════════
//  SYMMETRY — mirror and frame invariance of the full shot
// ══════════════════════════════════════════════════════════════════

fn mirror_z(v: Vec3) -> Vec3 {
    Vec3::new(v.x, v.y, -v.z)
}

/// Mirror of an angular velocity (a pseudovector) in the z = 0 plane.
fn mirror_spin_z(w: Vec3) -> Vec3 {
    Vec3::new(-w.x, -w.y, w.z)
}

#[test]
fn shot_mirror_symmetry_on_the_ground() {
    // ±sidespin and ±azimuth give mirrored z at every contact, every
    // ground point and at rest.
    for surface in [Surface::FAIRWAY, Surface::GREEN] {
        let right = simulate_shot(&lateral_input(3.0, 900.0), &BALL, &surface);
        let left = simulate_shot(&lateral_input(-3.0, -900.0), &BALL, &surface);
        assert!(right.ground.rest.pos.z > 1.0);
        assert_eq!(right.ground.contacts.len(), left.ground.contacts.len());
        assert_eq!(right.ground.points.len(), left.ground.points.len());
        for (r, l) in right.ground.contacts.iter().zip(&left.ground.contacts) {
            assert!((mirror_z(r.pos) - l.pos).length() < 1e-9);
            assert!((mirror_z(r.vel_out) - l.vel_out).length() < 1e-9);
            assert!((mirror_spin_z(r.spin_out) - l.spin_out).length() < 1e-7);
        }
        for (r, l) in right.ground.points.iter().zip(&left.ground.points) {
            assert!((mirror_z(r.pos) - l.pos).length() < 1e-9);
            assert_eq!(r.phase, l.phase);
        }
        assert!((mirror_z(right.ground.rest.pos) - left.ground.rest.pos).length() < 1e-9);
    }
}

#[test]
fn shot_frame_invariance_about_y() {
    // Rotating the launch azimuth (and with it the spin vector) about Y
    // on flat ground rotates the whole result: every contact and the rest
    // position, velocity and spin.
    let az = units::deg_to_rad(10.0);
    let rot = |v: Vec3| {
        Vec3::new(
            v.x * az.cos() - v.z * az.sin(),
            v.y,
            v.x * az.sin() + v.z * az.cos(),
        )
    };
    for side in [0.0, 800.0, -800.0] {
        let base = shot_on_fairway(&lateral_input(0.0, side));
        let turned = shot_on_fairway(&lateral_input(10.0, side));
        assert_eq!(base.ground.contacts.len(), turned.ground.contacts.len());
        for (b, t) in base.ground.contacts.iter().zip(&turned.ground.contacts) {
            assert!((rot(b.pos) - t.pos).length() < 1e-7);
            assert!((rot(b.vel_out) - t.vel_out).length() < 1e-7);
            assert!((rot(b.spin_out) - t.spin_out).length() < 1e-5);
        }
        assert!(
            (rot(base.ground.rest.pos) - turned.ground.rest.pos).length() < 1e-6,
            "sidespin {side}: rest {:?} vs rotated {:?}",
            turned.ground.rest.pos,
            rot(base.ground.rest.pos)
        );
    }
}

// ══════════════════════════════════════════════════════════════════
//  PINNED OUTPUTS — exact f64 results of three full shots
// ══════════════════════════════════════════════════════════════════

/// `[carry_m, apex_m, total_m, rest.x, rest.y, rest.z]` of
/// a driver on `FAIRWAY`, a 7-iron fade on `GREEN`, and a wedge draw onto
/// a tilted plane, bit for bit.
fn pinned_shots() -> [[f64; 6]; 3] {
    use crate::ground::{FnTerrain, GroundPhysics, GroundSample};

    let shot = |speed, launch, az, back, side| ShotInput {
        ball_speed_mph: speed,
        launch_angle_deg: launch,
        launch_azimuth_deg: az,
        backspin_rpm: back,
        sidespin_rpm: side,
    };
    let slope = GroundPhysics::new(FnTerrain(|x: f64, z: f64| GroundSample {
        height: -0.03 * x + 0.02 * z,
        normal: Vec3::new(0.03, 1.0, -0.02).normalize(),
        surface: Surface::GREEN,
    }));
    let results = [
        simulate_shot(
            &shot(167.0, 10.9, 0.0, 2686.0, 0.0),
            &BALL,
            &Surface::FAIRWAY,
        ),
        simulate_shot(
            &shot(120.0, 16.3, 1.5, 7097.0, 900.0),
            &BALL,
            &Surface::GREEN,
        ),
        simulate_shot(&shot(102.0, 24.2, -2.0, 9304.0, -700.0), &BALL, &slope),
    ];
    results.map(|r| {
        [
            r.flight.carry_m,
            r.flight.apex_m,
            r.total_m,
            r.ground.rest.pos.x,
            r.ground.rest.pos.y,
            r.ground.rest.pos.z,
        ]
    })
}

/// Pinned to 1e-9 m rather than bit for bit: the platform's `exp` may differ in the last bit.
#[test]
#[allow(
    clippy::unreadable_literal,
    reason = "f64 values as `{:?}` prints them"
)]
fn full_shot_outputs_are_pinned() {
    let expected = [
        [
            259.2803652763695,
            34.69155184333128,
            278.1377613735037,
            278.1377613735037,
            0.0,
            0.0,
        ],
        [
            151.63377813048913,
            28.141717089577387,
            159.26766027887996,
            158.44059205247387,
            0.0,
            16.210070967465583,
        ],
        [
            123.56988682054241,
            27.660157785411993,
            128.42230212531018,
            127.96555720954358,
            -4.055395782001129,
            -10.821453285741082,
        ],
    ];
    for (got, want) in pinned_shots()
        .iter()
        .flatten()
        .zip(expected.iter().flatten())
    {
        assert!((got - want).abs() < 1e-9, "{got} != {want}");
    }
}

// ── Wind ──

/// Sea-level air moving uniformly at `wind` (m/s).
fn windy(wind: Vec3) -> Tour {
    Tour::new(1.225, Wind::Uniform(wind))
}

fn driver() -> ShotInput {
    ShotInput {
        ball_speed_mph: 150.0,
        launch_angle_deg: 13.0,
        launch_azimuth_deg: 0.0,
        backspin_rpm: 2700.0,
        sidespin_rpm: 0.0,
    }
}

#[test]
#[allow(
    clippy::float_cmp,
    reason = "the flight and the ground phase are the same computation either way"
)]
fn ground_model_gets_the_shot_model() {
    // A custom ground model receives the shot's model, wind included: one
    // that forwards it settles exactly as the built-in surface does, and
    // one that settles with its own still-air model settles as in still
    // air.
    struct Forwarding;
    impl GroundModel for Forwarding {
        fn height(&self, x: f64, z: f64) -> f64 {
            Surface::FAIRWAY.height(x, z)
        }
        fn settle(&self, landing: &BallState, model: &dyn AeroModel) -> GroundResult {
            Surface::FAIRWAY.settle(landing, model)
        }
    }
    struct StillAir;
    impl GroundModel for StillAir {
        fn height(&self, x: f64, z: f64) -> f64 {
            Surface::FAIRWAY.height(x, z)
        }
        fn settle(&self, landing: &BallState, _model: &dyn AeroModel) -> GroundResult {
            Surface::FAIRWAY.settle(landing, &BALL)
        }
    }

    let input = driver();
    let model = windy(Vec3::new(-5.0, 0.0, 3.0));
    let windy = simulate_shot(&input, &model, &Surface::FAIRWAY);
    let forwarded = simulate_shot(&input, &model, &Forwarding);
    let still_air = simulate_shot(&input, &model, &StillAir);

    assert_eq!(forwarded.ground.rest.pos, windy.ground.rest.pos);
    // The flight is dimple's own, so the wind acts on it either way.
    assert_eq!(still_air.flight.carry_m, windy.flight.carry_m);
    assert_eq!(still_air.flight.lateral_m, windy.flight.lateral_m);
    assert_ne!(still_air.ground.rest.pos, windy.ground.rest.pos);
}

#[test]
fn headwind_shortens_carry_and_tailwind_lengthens_it() {
    let input = driver();
    let still = simulate_flight(&input, &BALL).carry_m;
    // +x is downrange, so a tailwind blows along +x and a headwind along −x.
    let head = simulate_flight(&input, &windy(Vec3::new(-5.0, 0.0, 0.0))).carry_m;
    let tail = simulate_flight(&input, &windy(Vec3::new(5.0, 0.0, 0.0))).carry_m;
    assert!(
        head < still,
        "headwind {head} should fall short of still {still}"
    );
    assert!(tail > still, "tailwind {tail} should pass still {still}");
    // Moving air is not symmetric: the headwind costs more than the same
    // tailwind gains, because drag grows with the square of the airspeed.
    assert!(
        still - head > tail - still,
        "head {head} tail {tail} still {still}"
    );
}

#[test]
fn crosswind_drifts_the_ball_downwind() {
    let input = driver();
    for (wind, sign) in [
        (Vec3::new(0.0, 0.0, 5.0), 1.0),
        (Vec3::new(0.0, 0.0, -5.0), -1.0),
    ] {
        let f = simulate_flight(&input, &windy(wind));
        assert!(
            f.lateral_m * sign > 1.0,
            "wind {wind:?} should drift the ball with it, lateral {}",
            f.lateral_m
        );
    }
}

#[test]
fn wind_carries_through_the_ground_phase() {
    let input = driver();
    let still = simulate_shot(&input, &BALL, &Surface::FAIRWAY).total_m;
    let head = simulate_shot(&input, &windy(Vec3::new(-5.0, 0.0, 0.0)), &Surface::FAIRWAY).total_m;
    let tail = simulate_shot(&input, &windy(Vec3::new(5.0, 0.0, 0.0)), &Surface::FAIRWAY).total_m;
    assert!(head < still, "head {head} still {still}");
    assert!(tail > still, "tail {tail} still {still}");
}

/// Wind on the bounce hops must not keep the ball bouncing: a fast landing
/// in a strong wind must not settle into a train of hops that the wind
/// feeds back as fast as the impacts drain it, which would carry the ball
/// on to the impact limit, hundreds of meters past the landing. (A wind
/// this strong at the turf may still blow the ball a long way along the
/// ground.)
#[test]
fn strong_wind_does_not_keep_the_ball_bouncing() {
    for (speed_mph, launch_deg, spin, wind) in [
        (67.1, 8.0, 1500.0, Vec3::new(20.0, 0.0, 0.0)),
        (67.1, 8.0, 1500.0, Vec3::new(14.0, 0.0, 14.0)),
        (100.7, 18.0, 10000.0, Vec3::new(-20.0, 0.0, 0.0)),
        (157.0, 12.0, 2500.0, Vec3::new(20.0, 0.0, 0.0)),
    ] {
        let input = ShotInput {
            ball_speed_mph: speed_mph,
            launch_angle_deg: launch_deg,
            launch_azimuth_deg: 0.0,
            backspin_rpm: spin,
            sidespin_rpm: 0.0,
        };
        let r = simulate_shot(&input, &windy(wind), &Surface::FAIRWAY);
        let hops = r.ground.contacts.len();
        assert_eq!(
            r.ground.termination,
            Termination::Rest,
            "{wind:?} {input:?}"
        );
        assert!(hops < 60, "{wind:?} {input:?}: {hops} impacts");
    }
}

/// A straight shot (mph, degrees, rpm).
fn straight(speed_mph: f64, launch_deg: f64, backspin_rpm: f64) -> ShotInput {
    ShotInput {
        ball_speed_mph: speed_mph,
        launch_angle_deg: launch_deg,
        launch_azimuth_deg: 0.0,
        backspin_rpm,
        sidespin_rpm: 0.0,
    }
}

#[test]
fn topped_shot_skims_forward_on_fairway() {
    // 130 mph at 2° with 1500 rpm lands at about 46 m/s and 3°. At so flat a
    // landing the contact angle is nearly zero, so the ball
    // skims on and rests ahead of its landing, never hopping above 0.5 m.
    let r = shot_on_fairway(&straight(130.0, 2.0, 1500.0));
    let landing = r.ground.contacts[0].pos;
    assert!(
        r.ground.rest.pos.x > landing.x,
        "rests at {:.1} m, landed at {:.1} m",
        r.ground.rest.pos.x,
        landing.x
    );
    let top = r
        .ground
        .points
        .iter()
        .map(|p| p.pos.y)
        .fold(f64::NEG_INFINITY, f64::max);
    assert!(top < 0.5, "ground phase reaches {top:.2} m");
}
