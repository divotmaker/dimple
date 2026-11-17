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
//! use dimple::units::meters_to_yards;
//! use dimple::{ShotInput, Surface, Tour, simulate_flight, simulate_shot};
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
//! let flight = simulate_flight(&input, &Tour::SEA_LEVEL);
//! assert!(meters_to_yards(flight.carry_m) > 200.0);
//! assert!(flight.lateral_m > 0.0); // positive = right of the target line
//!
//! let shot = simulate_shot(&input, &Tour::SEA_LEVEL, &Surface::FAIRWAY);
//! assert!(shot.total_m > shot.flight.carry_m);
//! assert!(shot.ground.rest.pos.z > flight.lateral_m); // keeps running right
//! ```
//!
//! # The air
//!
//! Every shot is flown by one [`AeroModel`]: the ball and the air it flies
//! through. [`Tour`] is dimple's model of a modern tour ball, built from an
//! air density and a [`Wind`]:
//!
//! | You know | Model |
//! |---|---|
//! | nothing (sea level, standard day, still air) | [`Tour::SEA_LEVEL`] |
//! | altitude and temperature | `Tour::new(standard_air_density(altitude_m, temperature_c), wind)` ([`standard_air_density`]) |
//! | a measured density (weather station: ρ = p / (287.05 · T)) | `Tour::new(air_density, wind)` |
//! | a wind | [`Wind::Uniform`], [`Wind::PowerLaw`] (a 10 m reading scaled with height), or [`Wind::custom`] (any function of time, position and height) |
//! | another ball, air or gravity | implement [`AeroModel`] yourself, or wrap [`Tour`] |
//!
//! [`Tour`]'s drag and lift coefficients, measured density and wind are
//! **beta**. A ball model of your own needs
//! four methods; see [`AeroModel`] for a drag-only ball and for a wrapper
//! that changes gravity in flight and on the ground alike.
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
//! **Parameter structs grow.** [`Surface`] and [`ground::GroundPhysics`]
//! gain fields in minor releases, so build them from a preset (or
//! [`ground::GroundPhysics::new`]) with `..`, e.g.
//! `Surface { mu_roll: 0.09, ..Surface::FAIRWAY }`.
//!
//! # Conventions
//!
//! **Units.** All internal computation is SI (meters, seconds, kilograms,
//! radians, rad/s). [`trajectory::ShotInput`] takes common launch-monitor
//! units (mph, degrees, rpm); results are SI, with the descent angle in
//! degrees. Conversions to display units live in [`units`].
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
//! (trajectory points, contacts, ground points, final rest position) is a
//! [`Vec3`] in this frame; convert it to your own math library's vector
//! through `[f64; 3]`. Positions are of the ball's center, and the ball touches
//! the ground when its center reaches the ground height: a ball on ground
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
//! | [`FlightResult::lateral_m`](trajectory::FlightResult::lateral_m) | ball lands right of target |
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
//! **Integration.** Flight uses fixed-step RK4 at 1 ms, with the spin
//! advanced by the model once per step ([`Tour`]: exponential decay); the
//! ground phase uses RK4 at 1 ms within each smooth phase. See
//! `docs/AERO.md` in the repository for the aerodynamic and ground models,
//! their sources and their validation.

/// The [`AeroModel`] trait and dimple's tour-ball model, [`Tour`].
pub mod aero;
/// Air density from altitude and temperature ([`standard_air_density`]) and
/// the wind ([`Wind`]).
pub mod air;
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
mod rk4;
/// Shot input, flight simulation, and full flight + ground shots.
pub mod trajectory;
/// Physical constants and unit conversions.
pub mod units;
mod vec3;

pub use aero::{AeroModel, Tour};
pub use air::{Wind, standard_air_density};
pub use ground::{BallState, Surface};
pub use trajectory::{
    ShotInput, simulate_flight, simulate_flight_from, simulate_shot, simulate_shot_from,
};
pub use vec3::Vec3;
