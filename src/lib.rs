//! Golf ball-flight physics.
//!
//! `dimple` turns launch-monitor data (ball speed, launch angle, azimuth,
//! backspin, sidespin) into a full shot: an RK4-integrated flight with drag,
//! Magnus lift and spin decay, followed by bounce and rollout on a turf
//! surface. Pure math: no I/O, no allocation beyond the returned trajectories,
//! no async runtime.
//!
//! # Example
//!
//! ```
//! use dimple::aero::BallModel;
//! use dimple::bounce::BounceSurface;
//! use dimple::environment::Environment;
//! use dimple::rollout::RolloutSurface;
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
//!     &BounceSurface::FAIRWAY,
//!     &RolloutSurface::FAIRWAY,
//! );
//! let stats = ShotStats::from_shot(&shot);
//! assert!(stats.total_yards > stats.carry_yards);
//! assert!(stats.lateral_yards > 0.0);
//! ```
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
//! Every position and velocity vector in the public API (trajectory points,
//! bounce events, rollout points, final rest position) uses this frame.
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
//! vector (right-hand rule). Positive backspin points along +Z (the top of the
//! ball moves back toward the golfer, producing lift). Positive sidespin is a
//! clockwise rotation viewed from above, i.e. along -Y, which curves the ball
//! toward +Z. Equivalently, a spin axis tilted to the right produces a
//! rightward curve. See [`aero::AeroParams::from_spin`].
//!
//! **Integration.** Flight uses fixed-step RK4 at 1 ms with exponential spin
//! decay; see `docs/AERO.md` in the repository for the aerodynamic model,
//! its sources and its validation against published Trackman data.

/// Ball models, Cd/Cl coefficients, and the drag/Magnus/gravity force model.
pub mod aero;
/// Ground impact: oblique bounce with turf compliance.
pub mod bounce;
/// Atmosphere: altitude and temperature to air density.
pub mod environment;
/// Fixed-step fourth-order Runge-Kutta integrator.
pub mod rk4;
/// Ground roll after the last bounce: sliding then pure rolling.
pub mod rollout;
/// Per-shot and per-session statistics in display units.
pub mod stats;
/// Shot input, flight simulation, and full flight + bounce + rollout shots.
pub mod trajectory;
/// Physical constants and unit conversions.
pub mod units;

#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;
