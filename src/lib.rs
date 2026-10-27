//! Golf ball-flight physics.
//!
//! `dimple` turns launch-monitor data (ball speed, launch angle, azimuth,
//! backspin, sidespin) into a full shot: an RK4-integrated flight with drag,
//! Magnus lift and spin decay, followed by the ground phase (impacts, hops,
//! slide and roll) with the full spin vector carried through every contact.
//! Impacts are checked against a published fit to bounce measurements; the
//! run and total distances of full shots are beta, compared with measured
//! data only for drives on a fairway, and may change between releases.
//! Pure math: no I/O, no allocation beyond the returned trajectories, no
//! async runtime.
//!
//! # Example
//!
//! ```
//! use dimple::aero::BallModel;
//! use dimple::environment::Environment;
//! use dimple::ground::Surface;
//! use dimple::stats::ShotStats;
//! use dimple::trajectory::{ShotInput, simulate_flight, simulate_shot};
//!
//! // A right-handed fade: starts 1° right, curves further right.
//! let input = ShotInput {
//!     ball_speed_mph: 150.0,
//!     launch_angle_deg: 12.5,
//!     launch_azimuth_deg: 1.0,
//!     backspin_rpm: 3200.0,
//!     sidespin_rpm: 600.0,
//! };
//!
//! let flight = simulate_flight(&input, &BallModel::TOUR, &Environment::SEA_LEVEL);
//! assert!(flight.carry_yards > 200.0);
//! assert!(flight.lateral_m > 0.0); // positive = right of the target line
//!
//! let shot = simulate_shot(
//!     &input,
//!     &BallModel::TOUR,
//!     &Environment::SEA_LEVEL,
//!     &Surface::FAIRWAY,
//! );
//! let stats = ShotStats::from_shot(&shot);
//! assert!(stats.total_yards > stats.carry_yards);
//! assert!(stats.lateral_yards > 0.0);
//! assert!(shot.final_pos.z > flight.lateral_m); // keeps running right
//! ```
//!
//! # Ground
//!
//! The ground API has three layers; use the simplest one that fits:
//!
//! 1. [`ground::Surface`]: flat ground made of one material. Presets
//!    `FAIRWAY`, `GREEN` and `FIRM_GREEN`; set a green's speed from a
//!    stimpmeter reading with [`ground::Surface::with_stimp`], e.g.
//!    `Surface::GREEN.with_stimp(12.0)`, and adjust any field, e.g.
//!    `Surface { mu_roll: 0.09, ..Surface::FAIRWAY }`.
//! 2. [`ground::GroundPhysics`] over a [`ground::Terrain`]: dimple's contact
//!    physics with tunable steps and limits, on ground you describe point by
//!    point: a heightfield giving the height, normal and surface at each
//!    (x, z), for example with [`ground::FnTerrain`]. Slopes, elevated
//!    greens and per-point surfaces all go here.
//! 3. [`ground::GroundModel`]: replace the physics entirely. It is
//!    object-safe; [`ground::impact`] and
//!    [`ground::GroundPhysics::step`] are public so a custom model can reuse
//!    pieces, and the result types have constructors.
//!
//! [`trajectory::simulate_shot`] launches from the terrain under the tee and
//! lands where the flight crosses the surface;
//! [`trajectory::simulate_flight`] always lands on the plane y = 0.
//!
//! # Conventions
//!
//! **Units.** All internal computation is SI (meters, seconds, kilograms,
//! radians, rad/s). [`trajectory::ShotInput`] takes common launch-monitor
//! units (mph, degrees, rpm); [`stats::ShotStats`] reports yards, feet and
//! degrees. Conversions live in [`units`].
//!
//! **Coordinate frame.** Right-handed, origin at the tee:
//!
//! | Axis | Direction |
//! |---|---|
//! | +X | downrange, along the target line |
//! | +Y | up |
//! | +Z | right of the target line, from the golfer's view |
//!
//! Every position, velocity and angular velocity vector in the public API
//! (trajectory points, contacts, ground points, final rest position) uses
//! this frame. Positions are of the ball's centre, and the ball touches
//! the ground when its centre reaches the ground height: a ball on ground
//! level with the tee has `y = 0`.
//!
//! **Lateral signs: positive = right.** Every lateral quantity is positive to
//! the right of the target line as seen by the golfer standing behind the
//! ball, matching the Flight Relay Protocol and common launch-monitor usage:
//!
//! | Quantity | Positive means |
//! |---|---|
//! | [`ShotInput::launch_azimuth_deg`](trajectory::ShotInput::launch_azimuth_deg) | ball starts right of target |
//! | [`ShotInput::sidespin_rpm`](trajectory::ShotInput::sidespin_rpm) | ball curves right (fade/slice for a right-handed golfer, draw/hook for a left-handed golfer) |
//! | [`FlightResult::lateral_m`](trajectory::FlightResult::lateral_m), [`ShotStats::lateral_yards`](stats::ShotStats::lateral_yards) | ball lands right of target |
//! | [`FlightResult::landing_sidespin_rads`](trajectory::FlightResult::landing_sidespin_rads) | same as input sidespin |
//! | `pos.z`, `vel.z` of any vector | right of / moving right of the target line |
//!
//! **Spin axis.** Backspin and sidespin combine into one angular velocity
//! vector (right-hand rule). Positive backspin points along the horizontal
//! axis perpendicular to the start line, `(−sin az, 0, cos az)` for launch
//! azimuth `az` (+Z for a straight shot): the top of the ball moves back
//! toward the golfer, producing lift. Positive sidespin is a clockwise
//! rotation viewed from above, i.e. along -Y, which curves the ball toward
//! +Z. Equivalently, a spin axis tilted to the right produces a rightward
//! curve. See [`trajectory::ShotInput::launch_spin_rads`].
//!
//! **Integration.** Flight uses fixed-step RK4 at 1 ms with exponential spin
//! decay; the ground phase uses RK4 at 1 ms within each smooth phase. See
//! `docs/AERO.md` in the repository for the aerodynamic and ground models,
//! their sources and their validation.

/// Ball models, Cd/Cl coefficients, and the drag/Magnus/gravity force model.
pub mod aero;
/// Atmosphere: altitude and temperature to air density.
pub mod environment;
/// Everything after landing: surfaces, heightfield terrain, 3D impacts,
/// hops, slide and roll on slopes, with the spin vector carried through
/// every contact.
///
/// Spin about the ground normal (launch-monitor sidespin, on flat ground)
/// cannot move a rigid ball's contact point, so it does not steer the roll;
/// pivot friction decays it. On a slope the normal is tilted, so part of
/// that vertical spin lies in the contact plane and does slip the contact
/// point. Lateral movement on the ground comes from the
/// spin components that do slip the contact point: a backspin axis that no
/// longer matches the heading of a curving ball (fades kick right, draws
/// left), spin along the direction of travel, and, at impact, the part of
/// any spin lying in Penner's tilted contact plane. Between impacts the spin
/// also curves the hops through Magnus lift, as in flight. The effects are
/// modest.
pub mod ground;
/// Fixed-step fourth-order Runge-Kutta integrators.
pub mod rk4;
/// Per-shot and per-session statistics in display units.
pub mod stats;
/// Shot input, flight simulation, and full flight + ground shots.
pub mod trajectory;
/// Physical constants and unit conversions.
pub mod units;

#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;
