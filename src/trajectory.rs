use glam::DVec3;

use crate::aero::{AeroParams, BallModel};
use crate::environment::Environment;
use crate::ground::{BallState, GroundModel, GroundResult, Surface};
use crate::rk4::{self, State};
use crate::units;

/// High-level shot input in common golf units.
///
/// Lateral conventions match launch-monitor and Flight Relay Protocol usage:
/// positive azimuth and positive sidespin both mean **right** of the target
/// line, from the golfer's view.
#[derive(Debug, Clone, Copy)]
pub struct ShotInput {
    /// Ball speed at launch (mph).
    pub ball_speed_mph: f64,
    /// Vertical launch angle above horizontal (degrees).
    pub launch_angle_deg: f64,
    /// Horizontal launch direction relative to the target line (degrees).
    /// Positive = right of target, negative = left.
    pub launch_azimuth_deg: f64,
    /// Backspin (rpm). Positive = backspin, negative = topspin.
    pub backspin_rpm: f64,
    /// Sidespin (rpm). Positive = ball curves right (fade/slice for a
    /// right-handed golfer), negative = curves left (draw/hook).
    pub sidespin_rpm: f64,
}

impl ShotInput {
    /// Launch angular velocity (rad/s, crate frame).
    ///
    /// With the launch heading `d = (cos az, 0, sin az)`, backspin is about
    /// the horizontal axis perpendicular to the start line, `d × ŷ =
    /// (−sin az, 0, cos az)`, and sidespin is about −Y:
    ///
    /// `ω = backspin · (−sin az, 0, cos az) − sidespin · ŷ`
    ///
    /// For a straight launch this is `(0, −sidespin, backspin)`. Backspin
    /// lift then stays in the vertical plane of the start line at any
    /// azimuth, so a pushed or pulled shot with no sidespin flies straight
    /// along its start line.
    #[must_use]
    pub fn launch_spin_rads(&self) -> DVec3 {
        let backspin = units::rpm_to_rads(self.backspin_rpm);
        let sidespin = units::rpm_to_rads(self.sidespin_rpm);
        let azimuth = units::deg_to_rad(self.launch_azimuth_deg);
        DVec3::new(
            -backspin * azimuth.sin(),
            -sidespin,
            backspin * azimuth.cos(),
        )
    }

    /// Typical 7-iron shot.
    #[must_use]
    pub fn seven_iron() -> Self {
        Self {
            ball_speed_mph: 132.0,
            launch_angle_deg: 16.0,
            launch_azimuth_deg: 0.0,
            backspin_rpm: 7000.0,
            sidespin_rpm: 0.0,
        }
    }

    /// Typical driver shot.
    #[must_use]
    pub fn driver() -> Self {
        Self {
            ball_speed_mph: 167.0,
            launch_angle_deg: 10.5,
            launch_azimuth_deg: 0.0,
            backspin_rpm: 2700.0,
            sidespin_rpm: 0.0,
        }
    }

    /// Typical wedge shot.
    #[must_use]
    pub fn wedge() -> Self {
        Self {
            ball_speed_mph: 100.0,
            launch_angle_deg: 28.0,
            launch_azimuth_deg: 0.0,
            backspin_rpm: 9500.0,
            sidespin_rpm: 0.0,
        }
    }
}

/// A single point along the flight trajectory.
///
/// Vectors use the crate frame: +X downrange, +Y up, +Z right of the target
/// line, origin at the tee.
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryPoint {
    /// Time since launch (s).
    pub time: f64,
    /// Position (m).
    pub pos: DVec3,
    /// Velocity (m/s).
    pub vel: DVec3,
}

/// Result of a flight simulation.
#[derive(Debug, Clone)]
pub struct FlightResult {
    /// Sampled trajectory from launch (t = 0) to landing (last point). The
    /// last point is the interpolated ground crossing, exactly on the ground.
    pub points: Vec<TrajectoryPoint>,
    /// Time from launch to landing (s).
    pub flight_time: f64,
    /// Straight-line ground distance from the tee to the landing point (m).
    pub carry_m: f64,
    /// [`carry_m`](Self::carry_m) in yards.
    pub carry_yards: f64,
    /// Maximum height above the tee (m).
    pub apex_m: f64,
    /// Lateral offset of the landing point from the target line (m).
    /// Positive = right, negative = left.
    pub lateral_m: f64,
    /// Angular velocity at landing (rad/s, crate frame).
    ///
    /// The spin axis is fixed during flight (only the rate decays), so this
    /// is the launch spin vector ([`ShotInput::launch_spin_rads`]) scaled by
    /// the decay.
    pub landing_spin: DVec3,
    /// Backspin at landing (rad/s): `landing_spin` projected on the launch
    /// backspin axis `(−sin az, 0, cos az)`, positive = backspin. Exact,
    /// because the axis does not move in flight.
    pub landing_backspin_rads: f64,
    /// Sidespin at landing (rad/s): `−landing_spin.y`. Same sign convention
    /// as [`ShotInput::sidespin_rpm`]: positive = curves right.
    pub landing_sidespin_rads: f64,
}

/// Simulate a ball flight from launch to landing on flat ground (y = 0).
///
/// Uses RK4 integration with a fixed 1 ms timestep. Records a trajectory
/// point every 10 ms, plus the landing point, which is interpolated to the
/// exact ground crossing. Spin decays exponentially over the flight.
#[must_use]
pub fn simulate_flight(input: &ShotInput, ball: &BallModel, env: &Environment) -> FlightResult {
    fly(input, ball, env, &Surface::FAIRWAY).0
}

/// Launch velocity (m/s): ball speed along the launch direction.
fn launch_velocity(input: &ShotInput) -> DVec3 {
    let speed = units::mph_to_mps(input.ball_speed_mph);
    let launch = units::deg_to_rad(input.launch_angle_deg);
    let azimuth = units::deg_to_rad(input.launch_azimuth_deg);
    DVec3::new(
        speed * launch.cos() * azimuth.cos(),
        speed * launch.sin(),
        speed * launch.cos() * azimuth.sin(),
    )
}

/// Flight from the tee at `(0, ground.height(0, 0), 0)` until
/// `ground.landing` reports a crossing (or 30 s). Returns the flight and the
/// full landing state.
fn fly<G: GroundModel + ?Sized>(
    input: &ShotInput,
    ball: &BallModel,
    env: &Environment,
    ground: &G,
) -> (FlightResult, BallState) {
    const DT: f64 = 0.001; // 1 ms timestep
    const SAMPLE_EVERY: u32 = 10; // record every 10 ms
    const MAX_STEPS: u32 = 30_000; // 30 seconds max

    let spin = input.launch_spin_rads();
    let mut aero = AeroParams::from_spin_vector(spin);
    let air_density = env.air_density();

    let tee = DVec3::new(0.0, ground.height(0.0, 0.0), 0.0);
    let mut prev = BallState {
        time: 0.0,
        pos: tee,
        vel: launch_velocity(input),
        spin,
    };
    let mut apex = 0.0_f64;
    let mut points = Vec::with_capacity((MAX_STEPS / SAMPLE_EVERY + 2) as usize);
    points.push(trajectory_point(&prev));

    let mut landing = None;
    for step_num in 1..=MAX_STEPS {
        let s = rk4::step(
            State {
                pos: prev.pos,
                vel: prev.vel,
            },
            prev.time,
            DT,
            |_pos, vel, _t| ball.acceleration(vel, &aero, air_density),
        );
        aero.decay_spin(DT);
        let next = BallState {
            time: f64::from(step_num) * DT,
            pos: s.pos,
            vel: s.vel,
            spin: aero.spin_vector(),
        };
        apex = apex.max(next.pos.y - tee.y);

        if let Some(l) = ground.landing(&prev, &next) {
            landing = Some(l);
            break;
        }
        if step_num % SAMPLE_EVERY == 0 {
            points.push(trajectory_point(&next));
        }
        prev = next;
    }

    let end = match landing {
        Some(l) => {
            points.push(trajectory_point(&l));
            l
        }
        None => {
            if points.last().is_none_or(|p| p.time < prev.time) {
                points.push(trajectory_point(&prev));
            }
            prev
        }
    };

    let carry_m = end.pos.x.hypot(end.pos.z);
    let azimuth = units::deg_to_rad(input.launch_azimuth_deg);
    let backspin_axis = DVec3::new(-azimuth.sin(), 0.0, azimuth.cos());

    let flight = FlightResult {
        points,
        flight_time: end.time,
        carry_m,
        carry_yards: units::meters_to_yards(carry_m),
        apex_m: apex,
        lateral_m: end.pos.z,
        landing_spin: end.spin,
        landing_backspin_rads: end.spin.dot(backspin_axis),
        landing_sidespin_rads: -end.spin.y,
    };
    (flight, end)
}

fn trajectory_point(s: &BallState) -> TrajectoryPoint {
    TrajectoryPoint {
        time: s.time,
        pos: s.pos,
        vel: s.vel,
    }
}

/// Result of a full shot simulation: flight, then everything on the ground.
#[derive(Debug, Clone)]
pub struct ShotResult {
    /// Airborne phase, launch to first ground contact.
    pub flight: FlightResult,
    /// Impacts, hops, slide and roll from the landing point to rest. Its
    /// times are absolute, so `flight.points` followed by `ground.points` is
    /// the whole path.
    pub ground: GroundResult,
    /// Horizontal distance from the tee to the final position (m).
    pub total_m: f64,
    /// [`total_m`](Self::total_m) in yards.
    pub total_yards: f64,
    /// Where the ball came to rest (m, crate frame: `final_pos.z` > 0 is
    /// right of the target line). Equal to `ground.rest.pos`.
    pub final_pos: DVec3,
}

/// Simulate a complete shot: flight, landing, and the ground phase.
///
/// `ground` decides where the flight lands ([`GroundModel::landing`]) and
/// what happens afterwards ([`GroundModel::settle`]). Pass a [`Surface`]
/// for flat ground with the default physics, a
/// [`GroundPhysics`](crate::ground::GroundPhysics) for per-point surfaces or
/// custom tunables, or your own [`GroundModel`] (also as `&dyn GroundModel`).
#[must_use]
pub fn simulate_shot<G: GroundModel + ?Sized>(
    input: &ShotInput,
    ball: &BallModel,
    env: &Environment,
    ground: &G,
) -> ShotResult {
    let (flight, landing) = fly(input, ball, env, ground);
    let settled = ground.settle(&landing, ball);
    let final_pos = settled.rest.pos;
    let total_m = final_pos.x.hypot(final_pos.z);

    ShotResult {
        flight,
        ground: settled,
        total_m,
        total_yards: units::meters_to_yards(total_m),
        final_pos,
    }
}

/// Interpolate trajectory position at a given time.
///
/// Returns `None` if `t` is outside the trajectory time range.
#[must_use]
pub fn interpolate_trajectory(points: &[TrajectoryPoint], t: f64) -> Option<DVec3> {
    if points.is_empty() || t < points[0].time {
        return None;
    }

    let last = points.last()?;
    if t >= last.time {
        return Some(last.pos);
    }

    // Binary search for the interval containing t
    let idx = points.partition_point(|p| p.time <= t);
    if idx == 0 {
        return Some(points[0].pos);
    }

    let a = &points[idx - 1];
    let b = &points[idx];
    let frac = (t - a.time) / (b.time - a.time);

    Some(a.pos.lerp(b.pos, frac))
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── Helper ──

    use crate::aero::BallModel;

    const BALL: BallModel = BallModel::TOUR;

    fn flight(speed_mph: f64, launch_deg: f64, backspin: f64, sidespin: f64) -> FlightResult {
        let input = ShotInput {
            ball_speed_mph: speed_mph,
            launch_angle_deg: launch_deg,
            launch_azimuth_deg: 0.0,
            backspin_rpm: backspin,
            sidespin_rpm: sidespin,
        };
        simulate_flight(&input, &BALL, &Environment::SEA_LEVEL)
    }

    fn flight_at(input: &ShotInput, env: &Environment) -> FlightResult {
        simulate_flight(input, &BALL, env)
    }

    /// Assert carry within ±tolerance yards of expected, with source citation.
    fn assert_carry(result: &FlightResult, expected: f64, tol: f64, label: &str) {
        let delta = (result.carry_yards - expected).abs();
        assert!(
            delta <= tol,
            "{label}: model {:.1} yds, expected {expected:.0}±{tol:.0} yds (Δ={delta:+.1})",
            result.carry_yards
        );
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
            let ke = 0.5 * BALL.mass_kg * p.vel.length_squared();
            let pe = BALL.mass_kg * crate::units::G * p.pos.y;
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
        let r = flight_at(&lateral_input(5.0, 0.0), &Environment::SEA_LEVEL);
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
        let r = flight_at(&lateral_input(-5.0, 0.0), &Environment::SEA_LEVEL);
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
        let r = flight_at(&lateral_input(0.0, 1500.0), &Environment::SEA_LEVEL);
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
        let r = flight_at(&lateral_input(0.0, -1500.0), &Environment::SEA_LEVEL);
        assert!(
            r.lateral_m < -5.0,
            "hook should finish left: {:.2}m",
            r.lateral_m
        );
    }

    #[test]
    fn push_draw_starts_right_and_curves_left() {
        let r = flight_at(&lateral_input(4.0, -800.0), &Environment::SEA_LEVEL);
        assert!(r.points[1].vel.z > 0.0, "should start right");
        assert!(
            r.points.last().unwrap().vel.z < 0.0,
            "should be moving left at landing"
        );
    }

    #[test]
    fn lateral_mirror_symmetry() {
        let right = flight_at(&lateral_input(3.0, 1200.0), &Environment::SEA_LEVEL);
        let left = flight_at(&lateral_input(-3.0, -1200.0), &Environment::SEA_LEVEL);
        assert!((right.lateral_m + left.lateral_m).abs() < 1e-9);
        assert!((right.carry_m - left.carry_m).abs() < 1e-9);
    }

    #[test]
    fn lateral_matches_landing_position() {
        let r = flight_at(&lateral_input(2.0, 900.0), &Environment::SEA_LEVEL);
        assert!((r.lateral_m - r.points.last().unwrap().pos.z).abs() < 1e-12);
    }

    #[test]
    fn landing_sidespin_keeps_input_sign() {
        let input = lateral_input(0.0, 1500.0);
        let r = flight_at(&input, &Environment::SEA_LEVEL);
        let initial = units::rpm_to_rads(1500.0);
        assert!(r.landing_sidespin_rads > 0.0 && r.landing_sidespin_rads < initial);
        let r = flight_at(&lateral_input(0.0, -1500.0), &Environment::SEA_LEVEL);
        assert!(r.landing_sidespin_rads < 0.0);
    }

    #[test]
    fn slice_shot_finishes_right() {
        let r = shot_on_fairway(&lateral_input(0.0, 1500.0));
        assert!(r.final_pos.z > 0.0, "final_pos.z={:.2}", r.final_pos.z);
        let stats = crate::stats::ShotStats::from_shot(&r);
        assert!(stats.lateral_yards > 0.0);
    }

    // ══════════════════════════════════════════════════════════════════
    //  CARRY DISTANCE — published Trackman / PING data (sea level)
    // ══════════════════════════════════════════════════════════════════
    //
    // Tolerance: ±10 yards. Wider than measurement error (~±2 yds) to absorb
    // differences in temperature normalization and source data across
    // publications; the per-club fit in docs/AERO.md is tighter.

    #[test]
    fn carry_pga_driver() {
        // Source: Trackman PGA Tour Averages 2024
        // 171 mph ball speed, 10.9° launch, 2686 rpm backspin → 282 yds carry
        let r = flight(171.0, 10.9, 2686.0, 0.0);
        assert_carry(&r, 282.0, 10.0, "PGA Tour driver (Trackman 2024)");
    }

    #[test]
    fn carry_trackman_classic_driver() {
        // Source: Trackman PGA Tour Averages (2009-2014), normalized 25°C sea level
        // 165 mph ball speed, 11.2° launch, 2685 rpm → 269 yds carry
        // Modern ball model overshoots classic-era data by design (~10 yds for driver).
        let r = flight(165.0, 11.2, 2685.0, 0.0);
        assert_carry(&r, 269.0, 12.0, "Trackman classic driver");
    }

    #[test]
    fn carry_trackman_classic_5iron() {
        // Source: Trackman PGA Tour Averages (2009-2014)
        // 132 mph ball speed, 12.1° launch, 5361 rpm → 194 yds carry
        let r = flight(132.0, 12.1, 5361.0, 0.0);
        assert_carry(&r, 194.0, 10.0, "Trackman classic 5-iron");
    }

    #[test]
    fn carry_trackman_classic_8iron() {
        // Source: Trackman PGA Tour Averages (2009-2014)
        // 115 mph ball speed, 18.1° launch, 7998 rpm → 160 yds carry
        let r = flight(115.0, 18.1, 7998.0, 0.0);
        assert_carry(&r, 160.0, 10.0, "Trackman classic 8-iron");
    }

    #[test]
    fn carry_trackman_classic_9iron() {
        // Source: Trackman PGA Tour Averages (2009-2014)
        // 109 mph ball speed, 20.4° launch, 8647 rpm → 148 yds carry
        let r = flight(109.0, 20.4, 8647.0, 0.0);
        assert_carry(&r, 148.0, 10.0, "Trackman classic 9-iron");
    }

    #[test]
    fn carry_trackman_7iron() {
        // Source: Trackman PGA Tour Averages (pre-2024, multiple secondary sources)
        // ~123 mph ball speed, 16.3° launch, 7097 rpm backspin → 176 yds carry
        let r = flight(123.0, 16.3, 7097.0, 0.0);
        assert_carry(&r, 176.0, 10.0, "Trackman PGA 7-iron");
    }

    #[test]
    fn carry_trackman_pw() {
        // Source: Trackman PGA Tour Averages (pre-2024, multiple secondary sources)
        // ~102 mph ball speed, ~25° launch, ~9300 rpm backspin → 142 yds carry
        // Wider tolerance: PW is the model's coldest club, a residual limit of
        // the single-sigmoid Cd/Cl model. See docs/AERO.md.
        let r = flight(102.0, 25.0, 9300.0, 0.0);
        assert_carry(&r, 142.0, 12.0, "Trackman PGA PW");
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
        let sea = flight_at(&driver, &Environment::SEA_LEVEL);
        let denver = flight_at(&driver, &Environment::DENVER);
        let gain_pct = (denver.carry_yards - sea.carry_yards) / sea.carry_yards * 100.0;
        // The TOUR model (saturating spin drag, cd_super=0.186) reproduces the
        // published altitude response: model ~7% vs Titleist ~6.1%.
        assert!(
            gain_pct > 4.0 && gain_pct < 9.0,
            "Denver gain {gain_pct:.1}%, expected ~6.1% (Titleist: elevation_ft × 0.00116)"
        );
    }

    #[test]
    fn altitude_increases_carry() {
        // Directional: higher altitude = more carry (less drag).
        let driver = ShotInput::driver();
        let sea = flight_at(&driver, &Environment::SEA_LEVEL);
        let denver = flight_at(&driver, &Environment::DENVER);
        assert!(denver.carry_yards > sea.carry_yards);
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
        let sea = flight_at(&fade, &Environment::SEA_LEVEL);
        let denver = flight_at(&fade, &Environment::DENVER);
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
        let hot = Environment {
            altitude_m: 0.0,
            temperature_c: 32.2,
        }; // 90°F
        let cold = Environment {
            altitude_m: 0.0,
            temperature_c: 10.0,
        }; // 50°F
        let r_hot = flight_at(&driver, &hot);
        let r_cold = flight_at(&driver, &cold);
        let delta = r_hot.carry_yards - r_cold.carry_yards;
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
        let hot = Environment {
            altitude_m: 0.0,
            temperature_c: 35.0,
        };
        let cold = Environment {
            altitude_m: 0.0,
            temperature_c: 5.0,
        };
        assert!(flight_at(&driver, &hot).carry_yards > flight_at(&driver, &cold).carry_yards);
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
        let mut aero = crate::aero::AeroParams::from_spin(2700.0, 0.0);
        let initial = aero.spin_rate;
        for _ in 0..6000 {
            aero.decay_spin(0.001);
        }
        let remaining_pct = aero.spin_rate / initial * 100.0;
        assert!(
            (remaining_pct - 78.0).abs() < 6.0,
            "after 6s: {remaining_pct:.1}% spin remaining, expected ~78% (Trackman)"
        );
    }

    // ══════════════════════════════════════════════════════════════════
    //  BALL ERA COMPARISON — same launch, different ball technology
    // ══════════════════════════════════════════════════════════════════
    //
    // These tests document how much ball technology affects carry distance.
    // Same launch conditions (speed, angle, spin), same environment (sea level),
    // but BallModel::TOUR (2020s) vs BallModel::PATENT_1997 (1990s).
    //
    // Modern tour balls carry ~15-25% farther than 1990s balls at the same
    // launch conditions. The primary driver is ~30% lower drag coefficient
    // from optimized dimple geometry. Secondary: spin-sensitive lift that
    // rewards the modern "high launch, low spin" strategy.
    //
    // Historical context: the USGA distance insight report (2020) attributes
    // roughly 30 yards of driver carry increase since 1995 to ball technology.
    // The model reproduces this: at 15°C sea level, TOUR driver ≈ 277 yds,
    // PATENT ≈ 248 yds, Δ ≈ 29 yards.

    const PATENT_BALL: BallModel = BallModel::PATENT_1997;

    fn flight_1997(speed_mph: f64, launch_deg: f64, backspin: f64) -> FlightResult {
        let input = ShotInput {
            ball_speed_mph: speed_mph,
            launch_angle_deg: launch_deg,
            launch_azimuth_deg: 0.0,
            backspin_rpm: backspin,
            sidespin_rpm: 0.0,
        };
        simulate_flight(&input, &PATENT_BALL, &Environment::SEA_LEVEL)
    }

    #[test]
    fn era_driver_carry_gap() {
        // Driver: 171 mph, 10.9°, 2686 rpm (Trackman PGA 2024 conditions).
        //
        // Modern tour ball: ~277 yds at 15°C sea level (validated above against
        // Trackman data).
        // 1990s patent ball: expected ~245-255 yds — the ~30 yard gap matches
        // USGA Distance Insights (2020) findings on ball technology contribution.
        //
        // The gap comes almost entirely from lower Cd (0.22 vs 0.30 at launch).
        // The Cl difference is small at driver SR≈0.08: modern Cl≈0.19, patent
        // Cl≈0.28. Higher patent Cl actually helps lift, partially offsetting
        // the drag penalty — but drag dominates.
        let modern = flight(171.0, 10.9, 2686.0, 0.0);
        let vintage = flight_1997(171.0, 10.9, 2686.0);
        let gap = modern.carry_yards - vintage.carry_yards;
        assert!(
            gap > 20.0 && gap < 50.0,
            "driver era gap should be 20-50 yds, got {gap:.1} (modern={:.1}, vintage={:.1})",
            modern.carry_yards,
            vintage.carry_yards
        );
    }

    #[test]
    fn era_7iron_carry_gap() {
        // 7-iron: 123 mph, 16.3°, 7097 rpm (Trackman PGA conditions).
        //
        // The gap is smaller for irons because high spin partially compensates:
        // 1990s balls had higher Cl (nearly saturated at all SR values), giving
        // more lift. But the drag penalty still dominates.
        //
        // Expected: ~5-25 yard gap — proportionally much less than driver.
        // Historically, iron distance gains came mostly from loft changes, not
        // ball aero; the modern ball's saturating spin drag (high late-flight
        // drag at iron SR) narrows the aero-only gap further.
        let modern = flight(123.0, 16.3, 7097.0, 0.0);
        let vintage = flight_1997(123.0, 16.3, 7097.0);
        let gap = modern.carry_yards - vintage.carry_yards;
        assert!(
            gap > 5.0 && gap < 40.0,
            "7-iron era gap should be 5-40 yds, got {gap:.1} (modern={:.1}, vintage={:.1})",
            modern.carry_yards,
            vintage.carry_yards
        );
    }

    #[test]
    fn era_wedge_carry_gap() {
        // PW: 102 mph, 25°, 9300 rpm (Trackman PGA conditions).
        //
        // Wedges have very high spin (SR≈0.46), so the patent ball's Cl advantage
        // is fully saturated (both models give Cl≈0.28). The gap is pure drag.
        // At this SR, spin-dependent drag is large: CD_SPIN × 0.46 adds ~0.07-0.09
        // to Cd in both models.
        //
        // With cl_sub=0.203 (modern) vs cl_sub=0.12 (patent), the modern ball has
        // more subcritical lift than the 1990s ball, but also more subcritical
        // drag (cd_sub=0.295 vs 0.19). At PW landing speeds (Re≈45-60k) the two
        // effects nearly cancel, so the era gap is very small for wedges.
        //
        // Expected: ~0-15 yard gap (may be near zero).
        let modern = flight(102.0, 25.0, 9300.0, 0.0);
        let vintage = flight_1997(102.0, 25.0, 9300.0);
        let gap = modern.carry_yards - vintage.carry_yards;
        assert!(
            gap > -1.0 && gap < 15.0,
            "PW era gap should be -1 to 15 yds, got {gap:.1} (modern={:.1}, vintage={:.1})",
            modern.carry_yards,
            vintage.carry_yards
        );
    }

    #[test]
    fn era_modern_always_farther() {
        // Sanity check: modern tour ball should carry farther than 1990s ball
        // at every club condition. The aerodynamic improvements help universally.
        let cases = [
            (171.0, 10.9, 2686.0, "driver"),
            (132.0, 12.1, 5361.0, "5-iron"),
            (123.0, 16.3, 7097.0, "7-iron"),
            (109.0, 20.4, 8647.0, "9-iron"),
            (102.0, 25.0, 9300.0, "PW"),
        ];
        for (speed, launch, spin, label) in cases {
            let modern = flight(speed, launch, spin, 0.0);
            let vintage = flight_1997(speed, launch, spin);
            assert!(
                modern.carry_yards > vintage.carry_yards,
                "{label}: modern ({:.1}) should carry farther than vintage ({:.1})",
                modern.carry_yards,
                vintage.carry_yards
            );
        }
    }

    #[test]
    fn era_gap_consistent_across_clubs() {
        // The era gap (modern minus vintage carry) should be positive for all clubs,
        // largest for driver (low spin, drag-dominated) and smaller for wedges
        // (high spin, subcritical lift similarity between eras).
        //
        // Source: USGA Distance Insights (2020) — ball technology contributes ~20-30
        // yards for driver over the past 25 years.
        let driver_modern = flight(171.0, 10.9, 2686.0, 0.0);
        let driver_vintage = flight_1997(171.0, 10.9, 2686.0);
        let iron_modern = flight(123.0, 16.3, 7097.0, 0.0);
        let iron_vintage = flight_1997(123.0, 16.3, 7097.0);

        let gap_driver = driver_modern.carry_yards - driver_vintage.carry_yards;
        let gap_iron = iron_modern.carry_yards - iron_vintage.carry_yards;

        // Driver gap should be larger than iron gap (drag reduction benefits low-spin more).
        assert!(
            gap_driver > gap_iron,
            "driver era gap ({gap_driver:.1}) should exceed iron ({gap_iron:.1})"
        );
        // Both gaps should be substantial.
        assert!(
            gap_driver > 15.0 && gap_iron > 5.0,
            "era gaps should be substantial: driver={gap_driver:.1}, iron={gap_iron:.1}"
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
            let r = flight_at(&input, &Environment::SEA_LEVEL);
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
        let f = flight_at(&input, &Environment::SEA_LEVEL);
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
        let expected = DVec3::new(0.0, -units::rpm_to_rads(600.0), units::rpm_to_rads(3000.0));
        assert!((straight - expected).length() < 1e-12);

        // Backspin axis is d × ŷ for launch heading d.
        let az = units::deg_to_rad(7.0);
        let pushed = lateral_input(7.0, 0.0).launch_spin_rads();
        let d = DVec3::new(az.cos(), 0.0, az.sin());
        let axis = d.cross(DVec3::Y);
        assert!((pushed - units::rpm_to_rads(3000.0) * axis).length() < 1e-9);
    }

    #[test]
    fn pushed_straight_shot_stays_on_its_start_line() {
        // A 6° push with no sidespin: the Magnus force stays in the vertical
        // plane through the start line, so the ball lands on that line.
        let input = lateral_input(6.0, 0.0);
        let az = units::deg_to_rad(6.0);
        let across = DVec3::new(-az.sin(), 0.0, az.cos()); // horizontal, ⟂ start line

        let aero = crate::aero::AeroParams::from_spin_vector(input.launch_spin_rads());
        let v = 60.0 * DVec3::new(az.cos() * 0.9, 0.3, az.sin() * 0.9).normalize();
        let gravity_free = BALL.acceleration(v, &aero, Environment::SEA_LEVEL.air_density())
            - DVec3::NEG_Y * units::G;
        assert!(
            gravity_free.dot(across).abs() < 1e-12,
            "aero force has a sideways component {:.3e}",
            gravity_free.dot(across)
        );

        let r = flight_at(&input, &Environment::SEA_LEVEL);
        for p in &r.points {
            assert!(
                p.pos.dot(across).abs() < 1e-9,
                "off the start line by {:.3e} m at t={:.2}",
                p.pos.dot(across),
                p.time
            );
        }
        assert!((r.lateral_m - r.carry_m * az.sin()).abs() < 1e-9);
        assert!(r.landing_sidespin_rads.abs() < 1e-12);
        assert!(r.landing_spin.dot(across) > 0.0);
        assert!((r.landing_spin.dot(across) - r.landing_backspin_rads).abs() < 1e-12);
    }

    #[test]
    fn landing_spin_projections_are_exact() {
        let input = lateral_input(-4.0, 900.0);
        let r = flight_at(&input, &Environment::SEA_LEVEL);
        let decay = r.landing_spin.length() / input.launch_spin_rads().length();
        assert!(decay > 0.0 && decay < 1.0);
        assert!(
            (r.landing_backspin_rads - units::rpm_to_rads(3000.0) * decay).abs() < 1e-9
                && (r.landing_sidespin_rads - units::rpm_to_rads(900.0) * decay).abs() < 1e-9
        );
    }

    // ══════════════════════════════════════════════════════════════════
    //  UNIFIED SHOT SIMULATION — flight + ground
    // ══════════════════════════════════════════════════════════════════

    use crate::ground::{GroundPhysics, Termination};

    fn shot_on_fairway(input: &ShotInput) -> ShotResult {
        simulate_shot(input, &BALL, &Environment::SEA_LEVEL, &Surface::FAIRWAY)
    }

    fn shot_on_green(input: &ShotInput) -> ShotResult {
        simulate_shot(input, &BALL, &Environment::SEA_LEVEL, &Surface::GREEN)
    }

    #[test]
    fn simulate_shot_accepts_dyn_ground_model() {
        let models: [&dyn GroundModel; 2] = [&Surface::GREEN, &GroundPhysics::new(Surface::GREEN)];
        let a = simulate_shot(
            &ShotInput::driver(),
            &BALL,
            &Environment::SEA_LEVEL,
            models[0],
        );
        let b = simulate_shot(
            &ShotInput::driver(),
            &BALL,
            &Environment::SEA_LEVEL,
            models[1],
        );
        let c = shot_on_green(&ShotInput::driver());
        assert_eq!(a.ground, c.ground);
        assert_eq!(b.ground, c.ground);
    }

    #[test]
    fn shot_total_greater_than_carry_driver() {
        // Driver on fairway: low spin, should roll forward significantly.
        let r = shot_on_fairway(&ShotInput::driver());
        assert!(
            r.total_yards > r.flight.carry_yards,
            "driver total ({:.1}) should exceed carry ({:.1})",
            r.total_yards,
            r.flight.carry_yards
        );
    }

    #[test]
    fn shot_driver_total_distance() {
        // PGA Tour driver: ~282 carry, ~305-315 total on fairway.
        let r = shot_on_fairway(&ShotInput::driver());
        assert!(
            r.total_yards > 290.0 && r.total_yards < 330.0,
            "driver total should be 290-330 yds, got {:.1} (carry={:.1}, rollout={:.1})",
            r.total_yards,
            r.flight.carry_yards,
            r.total_yards - r.flight.carry_yards
        );
    }

    #[test]
    fn shot_7iron_total_distance() {
        // PGA Tour 7-iron: 123 mph, 16.3°, 7097 rpm → ~176 carry, ~185-195 total.
        let input = ShotInput {
            ball_speed_mph: 123.0,
            launch_angle_deg: 16.3,
            launch_azimuth_deg: 0.0,
            backspin_rpm: 7097.0,
            sidespin_rpm: 0.0,
        };
        let r = shot_on_fairway(&input);
        assert!(
            r.total_yards > 178.0 && r.total_yards < 225.0,
            "7-iron total should be 178-225 yds, got {:.1} (carry={:.1}, rollout={:.1})",
            r.total_yards,
            r.flight.carry_yards,
            r.total_yards - r.flight.carry_yards
        );
    }

    #[test]
    fn shot_wedge_minimal_rollout_on_green() {
        // Wedge on green: high spin, should check up with minimal rollout.
        let r = shot_on_green(&ShotInput::wedge());
        let rollout_yards = r.total_yards - r.flight.carry_yards;
        assert!(
            rollout_yards.abs() < 15.0,
            "wedge on green rollout should be minimal, got {rollout_yards:.1} yds"
        );
    }

    #[test]
    fn shot_driver_first_bounce_height_in_penner_range() {
        // Penner (2002) §3.1 reports first-bounce heights from 2.82 m
        // (landing 68.0° from vertical) to 1.17 m (46.2°) for drives on his
        // green model. The driver preset lands about 49° from vertical.
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
        let hop = |v: DVec3| v.y * v.y / (2.0 * units::G);
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
    fn shot_final_pos_on_ground() {
        let r = shot_on_fairway(&ShotInput::seven_iron());
        assert!(
            r.final_pos.y.abs() < 1e-12,
            "final position should be on ground: y={:.6}",
            r.final_pos.y
        );
        assert_eq!(r.ground.rest.vel, DVec3::ZERO);
    }

    #[test]
    fn shot_final_pos_downrange() {
        let r = shot_on_fairway(&ShotInput::driver());
        assert!(
            r.final_pos.x > 0.0,
            "driver should finish downrange: x={:.1}",
            r.final_pos.x
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
        assert!((r.ground.points.last().unwrap().pos - r.final_pos).length() < 1e-12);
    }

    #[test]
    fn shot_landing_spin_decayed() {
        // Landing spin should be less than initial spin (decay during flight).
        let flight = simulate_flight(&ShotInput::driver(), &BALL, &Environment::SEA_LEVEL);
        let initial_backspin = units::rpm_to_rads(2700.0);
        assert!(
            flight.landing_backspin_rads < initial_backspin,
            "landing backspin ({:.1}) should be less than initial ({:.1})",
            flight.landing_backspin_rads,
            initial_backspin
        );
        assert!(
            flight.landing_backspin_rads > 0.0,
            "should still have some backspin"
        );
    }

    #[test]
    fn shot_green_less_rollout_than_fairway() {
        // Same drive on green vs fairway: green should have less total
        // rollout because its more compliant turf steepens each rebound and
        // takes more forward speed out of it, which outweighs the green's
        // lower rolling friction.
        let input = ShotInput::driver();
        let fairway = shot_on_fairway(&input);
        let green = shot_on_green(&input);
        let rollout_fairway = fairway.total_yards - fairway.flight.carry_yards;
        let rollout_green = green.total_yards - green.flight.carry_yards;
        assert!(
            rollout_green < rollout_fairway,
            "green rollout ({rollout_green:.1}) should be less than fairway ({rollout_fairway:.1})"
        );
    }

    #[test]
    fn shot_stats_from_shot() {
        let r = shot_on_fairway(&ShotInput::driver());
        let stats = crate::stats::ShotStats::from_shot(&r);
        assert!(stats.total_yards > stats.carry_yards);
        assert!(stats.rollout_yards > 0.0);
        assert!((stats.carry_yards + stats.rollout_yards - stats.total_yards).abs() < 0.1);
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
            let fade = simulate_shot(
                &driver_with_sidespin(600.0),
                &BALL,
                &Environment::SEA_LEVEL,
                &surface,
            );
            let draw = simulate_shot(
                &driver_with_sidespin(-600.0),
                &BALL,
                &Environment::SEA_LEVEL,
                &surface,
            );
            assert!(
                fade.final_pos.z > fade.flight.lateral_m && fade.flight.lateral_m > 0.0,
                "{name}: fade lands right and keeps moving right: land {:.2}, rest {:.2}",
                fade.flight.lateral_m,
                fade.final_pos.z
            );
            assert!(
                draw.final_pos.z < draw.flight.lateral_m && draw.flight.lateral_m < 0.0,
                "{name}: draw lands left and keeps moving left: land {:.2}, rest {:.2}",
                draw.flight.lateral_m,
                draw.final_pos.z
            );
            let stats = crate::stats::ShotStats::from_shot(&fade);
            assert!(stats.lateral_yards > 0.0);
        }
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
        // contact plane. On rigid ground (no tilt) the lever has no effect at
        // all; on compliant turf every lever that clears that spin at each
        // impact, from the presets' 1 mm up, gives the same path.
        let env = Environment::SEA_LEVEL;
        let rigid = Surface {
            compliance: 0.0,
            ..Surface::FAIRWAY
        };
        for side in [600.0, -600.0] {
            let input = driver_with_sidespin(side);
            for (base, levers) in [
                (rigid, [0.0, 0.001, 0.005]),
                (Surface::FAIRWAY, [0.001, 0.002, 0.005]),
                (Surface::GREEN, [0.001, 0.002, 0.005]),
            ] {
                let shots: Vec<ShotResult> = levers
                    .iter()
                    .map(|&pivot| simulate_shot(&input, &BALL, &env, &Surface { pivot, ..base }))
                    .collect();
                for (shot, pivot) in shots.iter().zip(levers).skip(1) {
                    assert!(
                        (shot.final_pos - shots[0].final_pos).length() < 1e-9,
                        "sidespin {side}, pivot {pivot}: rest {:?} vs {:?}",
                        shot.final_pos,
                        shots[0].final_pos
                    );
                    assert_eq!(shot.ground.contacts.len(), shots[0].ground.contacts.len());
                }
            }
            // With no pivot friction on compliant turf the spin survives to
            // the later impacts, whose tilt throws move the rest slightly.
            let free = simulate_shot(
                &input,
                &BALL,
                &env,
                &Surface {
                    pivot: 0.0,
                    ..Surface::FAIRWAY
                },
            );
            let preset = shot_on_fairway(&input);
            let shift = (free.final_pos - preset.final_pos).length();
            assert!(shift > 0.01 && shift < 0.2, "shift {shift:.3} m");
        }
    }

    // ══════════════════════════════════════════════════════════════════
    //  SYMMETRY — mirror and frame invariance of the full shot
    // ══════════════════════════════════════════════════════════════════

    fn mirror_z(v: DVec3) -> DVec3 {
        DVec3::new(v.x, v.y, -v.z)
    }

    /// Mirror of an angular velocity (a pseudovector) in the z = 0 plane.
    fn mirror_spin_z(w: DVec3) -> DVec3 {
        DVec3::new(-w.x, -w.y, w.z)
    }

    #[test]
    fn shot_mirror_symmetry_on_the_ground() {
        // ±sidespin and ±azimuth give mirrored z at every contact, every
        // ground point and at rest.
        for surface in [Surface::FAIRWAY, Surface::GREEN] {
            let right = simulate_shot(
                &lateral_input(3.0, 900.0),
                &BALL,
                &Environment::SEA_LEVEL,
                &surface,
            );
            let left = simulate_shot(
                &lateral_input(-3.0, -900.0),
                &BALL,
                &Environment::SEA_LEVEL,
                &surface,
            );
            assert!(right.final_pos.z > 1.0);
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
            assert!((mirror_z(right.final_pos) - left.final_pos).length() < 1e-9);
        }
    }

    #[test]
    fn shot_frame_invariance_about_y() {
        // Rotating the launch azimuth (and with it the spin vector) about Y
        // on flat ground rotates the whole result: every contact and the rest
        // position, velocity and spin.
        let az = units::deg_to_rad(10.0);
        let rot = |v: DVec3| {
            DVec3::new(
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
                (rot(base.final_pos) - turned.final_pos).length() < 1e-6,
                "sidespin {side}: rest {:?} vs rotated {:?}",
                turned.final_pos,
                rot(base.final_pos)
            );
        }
    }
}
