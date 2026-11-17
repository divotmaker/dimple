use crate::Vec3;

use crate::aero::AeroModel;
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
    /// Backspin (rpm). Positive = backspin, negative = topspin. A measured
    /// value: a ball with no spin gets no lift, so decide what an
    /// unmeasured spin should be before calling.
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
    /// Shot time (s): 0 at launch, and continuing from the start state's
    /// time for a flight from [`simulate_flight_from`].
    pub time: f64,
    /// Position (m).
    pub pos: Vec3,
    /// Velocity (m/s).
    pub vel: Vec3,
}

/// Result of a flight simulation.
///
/// The distances are measured from where the flight starts: the tee for
/// [`simulate_flight`], the start state for [`simulate_flight_from`].
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct FlightResult {
    /// Sampled trajectory from the start (the launch, at t = 0, for
    /// [`simulate_flight`]) to landing (last point), with absolute times.
    /// The last point is the interpolated ground crossing, exactly on the
    /// ground.
    pub points: Vec<TrajectoryPoint>,
    /// Time in the air (s): from launch, or from the start state's time, to
    /// landing.
    pub flight_time: f64,
    /// Horizontal distance from the start (the tee) to the landing point
    /// (m).
    pub carry_m: f64,
    /// Maximum height above the start (the tee) (m).
    pub apex_m: f64,
    /// Lateral offset (`z`) of the landing point from the start (from the
    /// target line, for a shot from the tee) (m). Positive = right,
    /// negative = left.
    pub lateral_m: f64,
    /// Angle of the landing velocity below horizontal (degrees); 90 for a
    /// ball falling straight down.
    pub descent_angle_deg: f64,
    /// Angular velocity at landing (rad/s, crate frame), as the model's
    /// [`spin_after`](AeroModel::spin_after) left it. With [`Tour`](crate::Tour)
    /// the spin axis is fixed during flight and only the rate decays, so
    /// this is the launch spin vector ([`ShotInput::launch_spin_rads`])
    /// scaled by the decay.
    pub landing_spin: Vec3,
}

/// Simulate a ball flight from launch to landing on flat ground (y = 0).
///
/// Uses RK4 integration with a fixed 1 ms timestep. Records a trajectory
/// point every 10 ms, plus the landing point, which is interpolated to the
/// exact ground crossing.
///
/// Gravity acts along −Y with the magnitude [`AeroModel::gravity`] of `model`
/// (standard gravity, [`units::G`], unless the model overrides it), added to
/// [`AeroModel::acceleration`] at every RK4 stage; the spin is advanced by
/// [`AeroModel::spin_after`] once per step. A flight that has not landed
/// after 120 s ends there, in the air. [`Tour`](crate::Tour) is dimple's model; its wind acts on the
/// ball's velocity relative to the air, so a headwind makes the ball fly as
/// if the air were faster and a crosswind curves it downwind. A uniform wind
/// is not a frame shift of the still-air trajectory: the relative airspeed
/// at launch changes, and with it the drag and the spin ratio, so the whole
/// flight is integrated through the moving air.
///
/// The same as [`simulate_flight_from`] with the launch state: time 0, at
/// the origin, with the launch velocity and
/// [`launch_spin_rads`](ShotInput::launch_spin_rads).
#[must_use]
pub fn simulate_flight(input: &ShotInput, model: &dyn AeroModel) -> FlightResult {
    simulate_flight_from(&launch_state(input, 0.0), model)
}

/// Simulate a ball flight from any state to landing on flat ground
/// (y = 0): a ball that hit a tree, a wall or a flagstick and flies on.
///
/// Integrates exactly as [`simulate_flight`] does, starting from `start`:
///
/// - Time continues from `start.time`; the points carry absolute times and
///   [`flight_time`](FlightResult::flight_time) is the time in the air.
/// - Positions stay in the crate frame (origin at the tee).
/// - [`carry_m`](FlightResult::carry_m) and
///   [`lateral_m`](FlightResult::lateral_m) are measured from the start's
///   horizontal position, [`apex_m`](FlightResult::apex_m) from its height.
/// - A start on the ground flies like a launch: if it moves into the ground
///   it lands within the first step, at its start.
/// - A start below the ground lands at once: the flight is one point, the
///   start placed on the ground (`y` set to the ground height) with its
///   velocity and spin, at `start.time`, so its flight time, carry and
///   apex are zero.
///
/// ```
/// use dimple::{BallState, ShotInput, Tour, Vec3, simulate_flight, simulate_flight_from};
///
/// // A ball 5 m up, 1.5 s into the shot, moving 20 m/s downrange.
/// let start = BallState {
///     time: 1.5,
///     pos: Vec3::new(40.0, 5.0, 0.0),
///     vel: Vec3::new(20.0, 2.0, 0.0),
///     spin: Vec3::new(0.0, 0.0, 300.0),
/// };
/// let flight = simulate_flight_from(&start, &Tour::SEA_LEVEL);
/// assert_eq!(flight.points[0].time, 1.5);
/// assert!(flight.carry_m > 0.0);
/// let landing = flight.points.last().unwrap();
/// assert_eq!(landing.pos.y, 0.0);
/// assert!((landing.pos.x - 40.0 - flight.carry_m).abs() < 1e-9);
/// ```
#[must_use]
pub fn simulate_flight_from(start: &BallState, model: &dyn AeroModel) -> FlightResult {
    fly(start, model, &Surface::FAIRWAY).0
}

/// The ball at launch: time 0 at `(0, tee_height, 0)`, with the launch
/// velocity and spin of `input`.
fn launch_state(input: &ShotInput, tee_height: f64) -> BallState {
    BallState {
        time: 0.0,
        pos: Vec3::new(0.0, tee_height, 0.0),
        vel: launch_velocity(input),
        spin: input.launch_spin_rads(),
    }
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

/// Flight from `start` until `ground.landing` reports a crossing (or 120 s).
/// Returns the flight and the full landing state.
///
/// The model is handed the ball's height above `ground` at every stage.
fn fly<G: GroundModel + ?Sized>(
    start: &BallState,
    model: &dyn AeroModel,
    ground: &G,
) -> (FlightResult, BallState) {
    const DT: f64 = 0.001; // 1 ms timestep
    const SAMPLE_EVERY: u32 = 10; // record every 10 ms
    const MAX_STEPS: u32 = 120_000; // 120 seconds max

    let floor = ground.height(start.pos.x, start.pos.z);
    if start.pos.y < floor {
        let mut landed = *start;
        landed.pos.y = floor;
        return (
            summarize(start, vec![trajectory_point(&landed)], 0.0, &landed),
            landed,
        );
    }

    let gravity = Vec3::new(0.0, -model.gravity(), 0.0);
    let mut prev = *start;
    let mut apex = 0.0_f64;
    let mut points = Vec::with_capacity(1024);
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
            |pos, vel, time| {
                let ball = BallState {
                    time,
                    pos,
                    vel,
                    spin: prev.spin,
                };
                gravity + model.acceleration(&ball, pos.y - ground.height(pos.x, pos.z))
            },
        );
        let next = BallState {
            time: start.time + f64::from(step_num) * DT,
            pos: s.pos,
            vel: s.vel,
            spin: model.spin_after(&prev, DT),
        };
        apex = apex.max(next.pos.y - start.pos.y);

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

    (summarize(start, points, apex, &end), end)
}

/// The flight from `start` to `end` along `points`.
fn summarize(
    start: &BallState,
    points: Vec<TrajectoryPoint>,
    apex: f64,
    end: &BallState,
) -> FlightResult {
    let horizontal_speed = end.vel.x.hypot(end.vel.z);
    let descent_angle_deg = if horizontal_speed > 1e-10 {
        units::rad_to_deg((-end.vel.y).atan2(horizontal_speed))
    } else {
        90.0
    };
    let lateral_m = end.pos.z - start.pos.z;
    FlightResult {
        points,
        flight_time: end.time - start.time,
        carry_m: (end.pos.x - start.pos.x).hypot(lateral_m),
        apex_m: apex,
        lateral_m,
        descent_angle_deg,
        landing_spin: end.spin,
    }
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
    /// Airborne phase, from the start (launch) to first ground contact.
    pub flight: FlightResult,
    /// Impacts, hops, slide and roll from the landing point to rest. Its
    /// times are absolute, so `flight.points` followed by `ground.points` is
    /// the whole path.
    pub ground: GroundResult,
    /// Horizontal distance from the start (the tee) to where the ball came
    /// to rest, `ground.rest.pos` (m).
    pub total_m: f64,
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
/// `model` flies the shot (see [`simulate_flight`]) and is handed to
/// [`GroundModel::settle`], so dimple's ground models fly the hops between
/// impacts through the same air and wind, and use the model's radius,
/// inertia factor and gravity for the contacts. The ground contacts
/// themselves are unchanged by the wind, because the turf does not move.
///
/// The same as [`simulate_shot_from`] with the launch state: time 0, at
/// `(0, ground.height(0, 0), 0)`, with the launch velocity and
/// [`launch_spin_rads`](ShotInput::launch_spin_rads).
#[must_use]
pub fn simulate_shot<G: GroundModel + ?Sized>(
    input: &ShotInput,
    model: &dyn AeroModel,
    ground: &G,
) -> ShotResult {
    simulate_shot_from(&launch_state(input, ground.height(0.0, 0.0)), model, ground)
}

/// Simulate a complete shot from any state: flight to `ground`, landing,
/// and the ground phase, as [`simulate_shot`] does from the tee.
///
/// The flight follows the rules of [`simulate_flight_from`] on `ground`
/// instead of the plane y = 0: time continues from `start.time`, positions
/// stay in the crate frame, carry and lateral offset are measured from the
/// start's horizontal position, and a start below `ground.height` lands at
/// once, placed on the ground. [`total_m`](ShotResult::total_m) is measured
/// from the start's horizontal position too.
///
/// A [`GroundModel`] whose [`landing`](GroundModel::landing) reports
/// obstacles (a tree trunk, a wall, a flagstick) can bounce the ball with
/// [`impact`](crate::ground::impact) and the obstacle's outward normal in
/// its [`settle`](GroundModel::settle), then continue with
/// `simulate_shot_from` from the rebound.
#[must_use]
pub fn simulate_shot_from<G: GroundModel + ?Sized>(
    start: &BallState,
    model: &dyn AeroModel,
    ground: &G,
) -> ShotResult {
    let (flight, landing) = fly(start, model, ground);
    let settled = ground.settle(&landing, model);
    let rest = settled.rest.pos;

    ShotResult {
        flight,
        total_m: (rest.x - start.pos.x).hypot(rest.z - start.pos.z),
        ground: settled,
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
