use crate::Vec3;

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
    pub fn launch_spin_rads(&self) -> Vec3 {
        let backspin = units::rpm_to_rads(self.backspin_rpm);
        let sidespin = units::rpm_to_rads(self.sidespin_rpm);
        let azimuth = units::deg_to_rad(self.launch_azimuth_deg);
        Vec3::new(
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
    pub pos: Vec3,
    /// Velocity (m/s).
    pub vel: Vec3,
}

/// Result of a flight simulation.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct FlightResult {
    /// Sampled trajectory from launch (t = 0) to landing (last point). The
    /// last point is the interpolated ground crossing, exactly on the ground.
    pub points: Vec<TrajectoryPoint>,
    /// Time from launch to landing (s).
    pub flight_time: f64,
    /// Horizontal distance from the tee to the landing point (m).
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
    pub landing_spin: Vec3,
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
///
/// `env` gives the air: its density, viscosity and wind
/// ([`Environment::wind`]). Drag and Magnus lift act on the ball's velocity
/// relative to the air, so a headwind makes the ball fly as if the air were
/// faster and a crosswind curves it downwind. A uniform wind is not a frame
/// shift of the still-air trajectory: the relative airspeed at launch
/// changes, and with it the drag and the spin ratio, so the whole flight is
/// integrated through the moving air.
#[must_use]
pub fn simulate_flight(input: &ShotInput, ball: &BallModel, env: &Environment) -> FlightResult {
    fly(input, ball, env, &Surface::FAIRWAY).0
}

/// Launch velocity (m/s): ball speed along the launch direction.
fn launch_velocity(input: &ShotInput) -> Vec3 {
    let speed = units::mph_to_mps(input.ball_speed_mph);
    let launch = units::deg_to_rad(input.launch_angle_deg);
    let azimuth = units::deg_to_rad(input.launch_azimuth_deg);
    Vec3::new(
        speed * launch.cos() * azimuth.cos(),
        speed * launch.sin(),
        speed * launch.cos() * azimuth.sin(),
    )
}

/// Flight from the tee at `(0, ground.height(0, 0), 0)` until
/// `ground.landing` reports a crossing (or 30 s). Returns the flight and the
/// full landing state.
///
/// The aerodynamic forces act on the ball's velocity relative to the air,
/// with the wind sampled at the ball's height above `ground`.
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
    let viscosity = env.dynamic_viscosity();
    let wind = env.wind;

    let tee = Vec3::new(0.0, ground.height(0.0, 0.0), 0.0);
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
            |pos, vel, _t| {
                ball.acceleration_with_viscosity(
                    wind.airspeed(vel, || pos.y - ground.height(pos.x, pos.z)),
                    &aero,
                    air_density,
                    viscosity,
                )
            },
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

    let end = if let Some(l) = landing {
        points.push(trajectory_point(&l));
        l
    } else {
        if points.last().is_none_or(|p| p.time < prev.time) {
            points.push(trajectory_point(&prev));
        }
        prev
    };

    let carry_m = end.pos.x.hypot(end.pos.z);
    let azimuth = units::deg_to_rad(input.launch_azimuth_deg);
    let backspin_axis = Vec3::new(-azimuth.sin(), 0.0, azimuth.cos());

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
#[non_exhaustive]
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
    pub final_pos: Vec3,
}

/// Simulate a complete shot: flight, landing, and the ground phase.
///
/// The ball launches from the ground under the tee,
/// `(0, ground.height(0, 0), 0)`. `ground` decides where the flight lands
/// ([`GroundModel::landing`]: on terrain, where the flight crosses the
/// surface, so elevated greens and downhill landings change the carry) and
/// what happens afterward ([`GroundModel::settle`]). Pass a [`Surface`]
/// for flat ground with the default physics, a
/// [`GroundPhysics`](crate::ground::GroundPhysics) for terrain (slopes,
/// heights, per-point surfaces) or custom tunables, or your own
/// [`GroundModel`] (also as `&dyn GroundModel`).
///
/// `env` gives the air for the flight (see [`simulate_flight`]) and is
/// handed to [`GroundModel::settle`], so dimple's ground models fly the hops
/// between impacts through the same air and wind. The ground contacts
/// themselves are unchanged by the wind, because the turf does not move.
#[must_use]
pub fn simulate_shot<G: GroundModel + ?Sized>(
    input: &ShotInput,
    ball: &BallModel,
    env: &Environment,
    ground: &G,
) -> ShotResult {
    let (flight, landing) = fly(input, ball, env, ground);
    let settled = ground.settle(&landing, ball, env);
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
/// Returns `None` for an empty trajectory or a `t` before its first point; a
/// `t` past the last point returns the last position.
#[must_use]
pub fn interpolate_trajectory(points: &[TrajectoryPoint], t: f64) -> Option<Vec3> {
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
mod tests;
