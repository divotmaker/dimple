use std::fmt;

use glam::DVec3;

use crate::aero::BallModel;
use crate::rk4::{self, Integrable};
use crate::units::G;

/// Penner (2002) eq. 8 reference contact angle θc (rad): 15.4°, his fit to
/// Haake's green impact.
const PENNER_THETA_C: f64 = 15.4 * std::f64::consts::PI / 180.0;
/// Impact speed (m/s) of Haake's measured green impact.
const PENNER_REF_SPEED: f64 = 18.6;
/// Incidence angle from the normal (rad) of Haake's measured green impact:
/// 44.4°.
const PENNER_REF_ANGLE: f64 = 44.4 * std::f64::consts::PI / 180.0;

/// Tangential speed (m/s) below which an impact has no tilt direction, so
/// Penner's contact angle is zero.
const MIN_TILT_SPEED: f64 = 0.01;

/// Contact-point slip speed (m/s) treated as zero. Slides end exactly at
/// zero slip; this only absorbs rounding.
const SLIP_EPS: f64 = 1e-9;

/// Spacing (s) of the sampled [`GroundPoint`]s.
const SAMPLE_INTERVAL: f64 = 0.01;

/// Time tolerance (s) for reaching sample times and the time limit.
const TIME_EPS: f64 = 1e-9;

/// Gravity vector (m/s²).
const GRAVITY: DVec3 = DVec3::new(0.0, -G, 0.0);

/// Normal coefficient of restitution between ball and turf.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Restitution {
    /// Constant normal coefficient of restitution (0..1).
    Constant(f64),
    /// Penner (2002) eq. 5, a speed-dependent fit to measured ball–turf
    /// restitution:
    /// e = 0.510 − 0.0375·v + 0.000903·v² for v ≤ 20 m/s, 0.120 above
    /// (v = normal impact speed against the effective contact plane, m/s).
    Penner2002,
}

impl Restitution {
    /// Coefficient of restitution at the given normal impact speed (m/s).
    ///
    /// The speed is taken as a magnitude; its sign is ignored.
    #[must_use]
    pub fn coefficient_at(&self, normal_speed: f64) -> f64 {
        match *self {
            Self::Constant(e) => e,
            Self::Penner2002 => {
                let v = normal_speed.abs();
                if v <= 20.0 {
                    0.510 - 0.0375 * v + 0.000_903 * v * v
                } else {
                    0.120
                }
            }
        }
    }
}

/// Material properties of a patch of ground.
///
/// Plain public fields with `const` presets; adjust one with struct update
/// syntax, e.g. `Surface { mu_roll: 0.09, ..Surface::FAIRWAY }`. A `Surface`
/// is also a flat [`Terrain`] (the plane y = 0 made of this surface) and a
/// complete [`GroundModel`] with the default [`GroundPhysics`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Surface {
    /// Normal coefficient of restitution at impact.
    pub restitution: Restitution,
    /// Coulomb friction coefficient between the ball and the turf, used at
    /// impact and while the ball slides.
    pub mu_slide: f64,
    /// Base rolling-resistance coefficient: the effective friction of a
    /// slowly rolling ball (stimpmeter speeds).
    pub mu_roll: f64,
    /// Speed scale (m/s) of the speed-dependent rolling resistance,
    /// μ_eff = `mu_roll` · (1 + (v / `roll_v_ref`)²) (Kensrud et al. 2024).
    pub roll_v_ref: f64,
    /// Turf compliance: dimensionless multiplier on Penner's (2002) eq. 8
    /// effective contact angle, θc = 15.4° · (v / 18.6 m/s) · (φ / 44.4°),
    /// with v the impact speed and φ the incidence angle from the normal.
    /// `1.0` is Penner's green; `0.0` is a rigid flat surface.
    pub compliance: f64,
    /// Torsional friction lever (m): the friction moment about the surface
    /// normal is `pivot` × normal force. It opposes spin about the normal,
    /// both at impact and on the ground, saturating at zero.
    ///
    /// It never changes the velocity of the impact it acts in, nor a slide
    /// or roll on level ground. The path depends on it only through spin
    /// about the normal carried into a later impact, where Penner's tilted
    /// contact plane turns part of that spin into slip. With the presets'
    /// 1 mm a full-shot landing removes all of it at once: a driver landing
    /// on `FAIRWAY` can remove up to about 155 rad/s (~1500 rpm), against
    /// the 47 rad/s a 600 rpm-sidespin drive lands with. Any lever that
    /// clears the spin at every impact gives the same path.
    ///
    /// The presets use 1 mm, the full-slip torsional couple of a Hertzian
    /// contact, (3π/16)·μ·a (Johnson 1985), with μ = 0.4 and a contact
    /// radius a ≈ 4 mm (a = 3–5 mm gives 0.7–1.2 mm). Uncalibrated.
    pub pivot: f64,
}

impl Surface {
    /// Standard fairway: Penner's restitution and friction with compliance
    /// 0.9, slightly below the green's so rebounds come off slightly firmer.
    /// The compliance is chosen, not measured; the run is sensitive to it
    /// (see `docs/AERO.md`).
    pub const FAIRWAY: Self = Self {
        restitution: Restitution::Penner2002,
        mu_slide: 0.40,
        mu_roll: 0.12,
        roll_v_ref: 2.5,
        compliance: 0.9,
        pivot: 0.001,
    };

    /// Standard green: Penner's (2002) green impact model with a rolling
    /// resistance of stimp ≈ 11.4 ft.
    pub const GREEN: Self = Self {
        restitution: Restitution::Penner2002,
        mu_slide: 0.40,
        mu_roll: 0.055,
        roll_v_ref: 2.5,
        compliance: 1.0,
        pivot: 0.001,
    };

    /// Effective rolling-resistance coefficient at `speed` (m/s):
    /// `mu_roll · (1 + (speed / roll_v_ref)²)` (Kensrud et al. 2024). Turf
    /// deformation grows with speed, so a fast-rolling ball decelerates much
    /// harder than a stimpmeter roll.
    #[must_use]
    pub fn rolling_friction(&self, speed: f64) -> f64 {
        let s = speed / self.roll_v_ref;
        self.mu_roll * (1.0 + s * s)
    }

    /// Effective contact angle β (rad) for an impact with tangential speed
    /// `v_t` and normal approach speed `v_n` (both ≥ 0).
    ///
    /// `compliance` times Penner's eq. 8 θc, capped at the incidence angle
    /// from the normal: at the cap the ball strikes the effective plane
    /// head-on, and a larger tilt would turn its tangential velocity
    /// backward. Zero when there is (almost) no tangential motion, since the
    /// tilt direction is undefined.
    fn contact_angle(&self, v_t: f64, v_n: f64) -> f64 {
        if v_t <= MIN_TILT_SPEED {
            return 0.0;
        }
        let speed = v_t.hypot(v_n);
        let incidence = v_t.atan2(v_n);
        let theta_c = PENNER_THETA_C * (speed / PENNER_REF_SPEED) * (incidence / PENNER_REF_ANGLE);
        (self.compliance * theta_c).min(incidence)
    }
}

/// The ground at one horizontal location.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundSample {
    /// Height (m) of the ground surface at (x, z). dimple treats the ball as
    /// touching when its centre reaches this height (the ball radius is not
    /// added), consistent with the flight's landing; on ground level with
    /// the tee this is 0 and a ball on the ground has `pos.y == 0`.
    pub height: f64,
    /// Unit outward surface normal, with `normal.y > 0`. [`GroundPhysics`]
    /// currently requires +Y (see [`Terrain`]).
    pub normal: DVec3,
    /// Material at this location.
    pub surface: Surface,
}

/// Layer 1 of the ground API: "what is the ground like here?"
///
/// A heightfield point query. dimple never needs to know how the ground is
/// stored: an analytic plane, a grid, a mesh lookup or an engine query all
/// fit behind `sample`. [`GroundPhysics`] samples the terrain under the
/// ball's contact point every step, so surfaces can change mid-roll.
///
/// [`GroundPhysics`] currently requires level terrain: every sample must
/// have normal +Y and the same height (0 for ground level with the tee).
/// Other terrain is not rejected but is unsupported until slopes are
/// modeled: the contact physics would use the sampled normal while the
/// slide and roll keep the ball at its landing height.
pub trait Terrain {
    /// The ground at horizontal position (`x`, `z`) in the crate frame.
    #[must_use]
    fn sample(&self, x: f64, z: f64) -> GroundSample;
}

/// A surface on its own is flat ground: the plane y = 0, normal +Y, made of
/// this surface everywhere.
impl Terrain for Surface {
    fn sample(&self, _x: f64, _z: f64) -> GroundSample {
        GroundSample {
            height: 0.0,
            normal: DVec3::Y,
            surface: *self,
        }
    }
}

/// A [`Terrain`] from a closure `|x, z| -> GroundSample`.
///
/// ```
/// use dimple::ground::{FnTerrain, GroundPhysics, GroundSample, Surface};
/// use glam::DVec3;
///
/// // A green from 150 m to 170 m downrange, fairway everywhere else.
/// let range = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| GroundSample {
///     height: 0.0,
///     normal: DVec3::Y,
///     surface: if (150.0..170.0).contains(&x) { Surface::GREEN } else { Surface::FAIRWAY },
/// }));
/// # let _ = range;
/// ```
#[derive(Clone, Copy)]
pub struct FnTerrain<F>(pub F);

impl<F: Fn(f64, f64) -> GroundSample> Terrain for FnTerrain<F> {
    fn sample(&self, x: f64, z: f64) -> GroundSample {
        (self.0)(x, z)
    }
}

impl<F> fmt::Debug for FnTerrain<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("FnTerrain").finish_non_exhaustive()
    }
}

/// Full kinematic state of the ball.
///
/// Complete by design: time, position, velocity and spin are the whole
/// rigid-body state, so this is a plain struct to build and match on
/// freely.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BallState {
    /// Absolute shot time (s since launch).
    pub time: f64,
    /// Position of the ball's centre (m, crate frame).
    pub pos: DVec3,
    /// Velocity (m/s).
    pub vel: DVec3,
    /// Angular velocity (rad/s, right-hand rule, crate frame).
    pub spin: DVec3,
}

/// What the ball is doing between ground points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    /// In the air between impacts (gravity only).
    Hop,
    /// On the ground with the contact point slipping: sliding friction.
    Slide,
    /// On the ground, rolling without slipping: rolling resistance.
    Roll,
}

/// Why the ground phase ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Termination {
    /// The ball came to rest.
    Rest,
    /// The time limit ([`GroundPhysics::max_time`]) was reached first.
    MaxTime,
    /// The impact limit ([`GroundPhysics::max_impacts`]) was reached while
    /// the ball was still bouncing.
    MaxImpacts,
}

/// One impact with the ground.
///
/// Non-exhaustive so it can gain fields; build one with [`Contact::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct Contact {
    /// Absolute shot time (s).
    pub time: f64,
    /// Ball centre at impact (m).
    pub pos: DVec3,
    /// Unit surface normal at the contact.
    pub normal: DVec3,
    /// Surface at the contact point.
    pub surface: Surface,
    /// Velocity arriving at the ground (m/s).
    pub vel_in: DVec3,
    /// Velocity leaving the contact (m/s). For the last contact before the
    /// ball stays on the ground this is the velocity it slides or rolls
    /// away with: the rebound projected onto the ground plane, with a slip
    /// at or below [`GroundPhysics::roll_threshold`] removed.
    pub vel_out: DVec3,
    /// Angular velocity arriving (rad/s).
    pub spin_in: DVec3,
    /// Angular velocity leaving (rad/s); for the last contact, the spin the
    /// ball slides or rolls away with.
    pub spin_out: DVec3,
}

impl Contact {
    /// A contact from the ball arriving at the ground (`arrival`: time,
    /// position, incoming velocity and spin), the ground there (`at`:
    /// normal and surface) and the ball leaving it (`departure`: outgoing
    /// velocity and spin; its time and position are not used).
    #[must_use]
    pub fn new(arrival: &BallState, at: &GroundSample, departure: &BallState) -> Self {
        Self {
            time: arrival.time,
            pos: arrival.pos,
            normal: at.normal,
            surface: at.surface,
            vel_in: arrival.vel,
            vel_out: departure.vel,
            spin_in: arrival.spin,
            spin_out: departure.spin,
        }
    }
}

/// A sampled point of the ground phase.
///
/// Non-exhaustive so it can gain fields; build one with
/// [`GroundPoint::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct GroundPoint {
    /// Absolute shot time (s).
    pub time: f64,
    /// Ball centre (m).
    pub pos: DVec3,
    /// Velocity (m/s).
    pub vel: DVec3,
    /// Angular velocity (rad/s).
    pub spin: DVec3,
    /// What the ball does from this point until the next one.
    pub phase: Phase,
}

impl GroundPoint {
    /// The ball in `state`, doing `phase` until the next point.
    #[must_use]
    pub fn new(state: &BallState, phase: Phase) -> Self {
        Self {
            time: state.time,
            pos: state.pos,
            vel: state.vel,
            spin: state.spin,
            phase,
        }
    }
}

/// Everything from landing to rest.
///
/// Non-exhaustive so it can gain fields; build one with
/// [`GroundResult::new`].
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct GroundResult {
    /// Impacts in order, starting with the landing.
    pub contacts: Vec<Contact>,
    /// Sampled path (every 10 ms, at each contact, at each phase change, and
    /// at the end), in time order. Times are absolute shot times, so a full
    /// shot's path is `flight.points` followed by these.
    pub points: Vec<GroundPoint>,
    /// Where the simulation stopped. A true rest state only when
    /// [`termination`](Self::termination) is [`Termination::Rest`]; spin
    /// about the normal that the pivot friction has not yet removed is kept.
    pub rest: BallState,
    /// Why the simulation stopped.
    pub termination: Termination,
}

impl GroundResult {
    /// A result from its parts; see the fields for what each holds.
    #[must_use]
    pub fn new(
        contacts: Vec<Contact>,
        points: Vec<GroundPoint>,
        rest: BallState,
        termination: Termination,
    ) -> Self {
        Self {
            contacts,
            points,
            rest,
            termination,
        }
    }
}

/// Layer 2 of the ground API: dimple's contact physics over any
/// [`Terrain`].
///
/// Each impact uses [`impact`]; the ball keeps hopping (gravity-only flight,
/// RK4) while its rebound would rise at least
/// [`min_hop_height`](Self::min_hop_height), then slides and rolls
/// ([`step`](Self::step)) until it stops. The tunables are public fields; start
/// from [`GroundPhysics::new`] and override with struct update syntax.
#[derive(Debug, Clone, Copy)]
pub struct GroundPhysics<T: Terrain> {
    /// Where the ball is and what it lands on.
    pub terrain: T,
    /// Integration step (s) for hops, slides and rolls. Default 1 ms.
    pub dt: f64,
    /// A rolling ball slower than this (m/s) comes to rest. It also stops
    /// once its speed is below one step's deceleration. Either way it
    /// finishes at constant deceleration, so this barely moves the rest
    /// position. Default 0.01 m/s.
    pub stop_speed: f64,
    /// Contact-point slip speed (m/s) at or below which a ball that has
    /// stopped bouncing starts out rolling instead of sliding. The residual
    /// slip is removed by the friction impulse that would have finished the
    /// slide, so the rolling velocity is the same either way. Slides are
    /// otherwise stepped to exactly zero slip. Default 0.01 m/s.
    pub roll_threshold: f64,
    /// The ball keeps bouncing while its rebound would rise at least this
    /// high (m) above the ground: Penner's (2002) bounce-to-roll transition,
    /// 5 mm. Penner tested 1–10 mm and found only a minor effect on the run.
    pub min_hop_height: f64,
    /// Safety bound on impacts. Every hop rises at least
    /// [`min_hop_height`](Self::min_hop_height) and loses energy, so real
    /// shots stop bouncing far below this; it guards degenerate surfaces
    /// (e.g. `Restitution::Constant(1.0)` with no friction). Default 200.
    pub max_impacts: usize,
    /// Time limit (s) on the ground phase, from landing. Default 30 s.
    pub max_time: f64,
}

impl<T: Terrain> GroundPhysics<T> {
    /// Default physics over `terrain`: 1 ms steps, 5 mm hop transition,
    /// 0.01 m/s stop and roll thresholds, at most 200 impacts and 30 s.
    #[must_use]
    pub fn new(terrain: T) -> Self {
        Self {
            terrain,
            dt: 0.001,
            stop_speed: DEFAULT_STOP_SPEED,
            roll_threshold: DEFAULT_ROLL_THRESHOLD,
            min_hop_height: 0.005,
            max_impacts: 200,
            max_time: 30.0,
        }
    }

    /// Advance a ball that is on the ground by one step of
    /// [`dt`](Self::dt), on the terrain under its contact point.
    ///
    /// Returns the new state and the phase the ball was in over the step.
    /// RK4 within each smooth phase:
    ///
    /// - **Slide** (the contact point slips): friction `−μ_slide·m·g·û` and
    ///   its torque `(−R n) × F`. On flat ground the slip speed falls
    ///   linearly at `((1+k)/k)·μ_slide·g`, so a slide that would end within
    ///   the step stops exactly at the slide-to-roll time (the step is
    ///   shortened) and the ball leaves rolling.
    /// - **Roll**: spin tangential to the ground is locked to the velocity;
    ///   deceleration `μ_eff·g/(1+k)` with `μ_eff` from
    ///   [`Surface::rolling_friction`]. A ball whose speed is below one
    ///   step's deceleration (or [`stop_speed`](Self::stop_speed)) comes to
    ///   rest: it finishes at constant deceleration, the returned velocity is
    ///   exactly zero, and the time advances by the stopping time only.
    /// - In both, pivot friction decays the spin about the normal at
    ///   `pivot·g/(k R²)`, saturating at zero.
    ///
    /// The velocity is kept in the ground plane. Unlike
    /// [`settle`](GroundModel::settle), this does not apply
    /// [`max_time`](Self::max_time).
    #[must_use]
    pub fn step(&self, state: &BallState, ball: &BallModel) -> (BallState, Phase) {
        let at = self.sample_contact(state.pos, ball.radius_m());
        contact_step(state, &at, ball, self.dt, self.stop_speed)
    }

    /// Terrain under the contact point of a ball centred at `pos`.
    fn sample_contact(&self, pos: DVec3, radius: f64) -> GroundSample {
        let below = self.terrain.sample(pos.x, pos.z);
        if below.normal == DVec3::Y {
            return below;
        }
        let contact = pos - radius * below.normal;
        self.terrain.sample(contact.x, contact.z)
    }

    /// Gravity-only flight from a rebound until the ball meets the ground
    /// again: `Ok(arrival)`, or `Err(state)` if the time limit `t_end` comes
    /// first. Sampled points are appended to `points`.
    fn hop(
        &self,
        start: &BallState,
        t_end: f64,
        points: &mut Vec<GroundPoint>,
    ) -> Result<BallState, BallState> {
        let mut s = rk4::State {
            pos: start.pos,
            vel: start.vel,
        };
        let mut t = start.time;
        let mut next_sample = t + SAMPLE_INTERVAL;
        loop {
            if t >= t_end - TIME_EPS {
                return Err(BallState {
                    time: t,
                    pos: s.pos,
                    vel: s.vel,
                    spin: start.spin,
                });
            }
            let h = self.dt.min(t_end - t);
            let b = rk4::step(s, t, h, |_, _, _| GRAVITY);
            let below = self.terrain.sample(b.pos.x, b.pos.z).height;
            if b.pos.y <= below {
                return Ok(self.hop_arrival(s, t, h, start.spin));
            }
            s = b;
            t += h;
            if t >= next_sample - TIME_EPS {
                let sample = BallState {
                    time: t,
                    pos: s.pos,
                    vel: s.vel,
                    spin: start.spin,
                };
                points.push(GroundPoint::new(&sample, Phase::Hop));
                next_sample += SAMPLE_INTERVAL;
            }
        }
    }

    /// Exact ground crossing within a gravity-only step of length `h` from
    /// `a` (at time `t`), which ends at or below the ground.
    ///
    /// Within the step the motion is the exact parabola RK4 reproduces, so
    /// the crossing time solves `y(τ) = height(x(τ), z(τ))`; the height is
    /// iterated to a fixed point (one pass on level ground).
    fn hop_arrival(&self, a: rk4::State, t: f64, h: f64, spin: DVec3) -> BallState {
        let at = |tau: f64| a.pos + a.vel * tau + 0.5 * GRAVITY * tau * tau;
        let end = at(h);
        let mut ground = self.terrain.sample(end.x, end.z).height;
        let mut tau = h;
        for _ in 0..8 {
            tau = drop_time(a.pos.y - ground, a.vel.y).clamp(0.0, h);
            let p = at(tau);
            let next = self.terrain.sample(p.x, p.z).height;
            let converged = (next - ground).abs() <= 1e-12;
            ground = next;
            if converged {
                break;
            }
        }
        let mut pos = at(tau);
        pos.y = ground;
        BallState {
            time: t + tau,
            pos,
            vel: a.vel + GRAVITY * tau,
            spin,
        }
    }
}

/// Time (s) for a ball `height` above the ground moving up at `vy` to fall
/// back to it under gravity alone. The positive root of
/// `height + vy·τ − g·τ²/2 = 0`, written to avoid cancellation.
fn drop_time(height: f64, vy: f64) -> f64 {
    if height <= 0.0 && vy <= 0.0 {
        return 0.0;
    }
    let disc = (vy * vy + 2.0 * G * height.max(0.0)).sqrt();
    if vy >= 0.0 {
        (vy + disc) / G
    } else {
        2.0 * height.max(0.0) / (disc - vy)
    }
}

const DEFAULT_STOP_SPEED: f64 = 0.01;
const DEFAULT_ROLL_THRESHOLD: f64 = 0.01;

/// Layer 3 of the ground API: the whole ground model.
///
/// [`Surface`] (flat ground) and [`GroundPhysics`] (any [`Terrain`])
/// implement it; implement it yourself to replace the physics entirely.
/// Object-safe, so `&dyn GroundModel` works with
/// [`simulate_shot`](crate::trajectory::simulate_shot).
pub trait GroundModel {
    /// Height (m) of the ground surface at (`x`, `z`); the ball touches
    /// when its centre reaches it (see [`GroundSample::height`]). Used by
    /// the default [`landing`](Self::landing) and to place the ball on the
    /// tee at (0, `height(0, 0)`, 0).
    #[must_use]
    fn height(&self, x: f64, z: f64) -> f64;

    /// Where the flight segment from `a` to `b` (one integration step) first
    /// meets the ground, if it does.
    ///
    /// The default interpolates linearly in `y − height` between the two
    /// states and returns the interpolated time, velocity and spin, with the
    /// position placed exactly on the ground. A segment that starts on the
    /// ground and ends above it does not land. An engine can override this
    /// with its own sweep test.
    #[must_use]
    fn landing(&self, a: &BallState, b: &BallState) -> Option<BallState> {
        let fa = a.pos.y - self.height(a.pos.x, a.pos.z);
        let fb = b.pos.y - self.height(b.pos.x, b.pos.z);
        if fb > 0.0 || fa < 0.0 {
            return None;
        }
        let frac = if fa > fb { fa / (fa - fb) } else { 0.0 };
        let mut pos = a.pos.lerp(b.pos, frac);
        pos.y = self.height(pos.x, pos.z);
        Some(BallState {
            time: a.time + (b.time - a.time) * frac,
            pos,
            vel: a.vel.lerp(b.vel, frac),
            spin: a.spin.lerp(b.spin, frac),
        })
    }

    /// Everything from the landing state to rest.
    #[must_use]
    fn settle(&self, landing: &BallState, ball: &BallModel) -> GroundResult;
}

/// Flat ground with the default [`GroundPhysics`].
impl GroundModel for Surface {
    fn height(&self, _x: f64, _z: f64) -> f64 {
        0.0
    }

    fn settle(&self, landing: &BallState, ball: &BallModel) -> GroundResult {
        GroundPhysics::new(*self).settle(landing, ball)
    }
}

impl<T: Terrain> GroundModel for GroundPhysics<T> {
    fn height(&self, x: f64, z: f64) -> f64 {
        self.terrain.sample(x, z).height
    }

    fn settle(&self, landing: &BallState, ball: &BallModel) -> GroundResult {
        let t_end = landing.time + self.max_time;
        let radius = ball.radius_m();
        let mut contacts = Vec::new();
        let mut points = Vec::new();
        let mut state = *landing;

        // Impacts and the hops between them.
        loop {
            let at = self.sample_contact(state.pos, radius);
            let n = at.normal;
            let out = impact(&state, &at, ball);
            let v_n = out.vel.dot(n);
            let hops = v_n > 0.0 && v_n * v_n >= 2.0 * G * n.y * self.min_hop_height;
            if !hops {
                // The ball stays on the ground: drop the (small or inward)
                // normal velocity and slide/roll from here.
                let mut leave = BallState {
                    vel: out.vel - v_n * n,
                    ..out
                };
                if slip(leave.vel, leave.spin, n, radius).length() <= self.roll_threshold {
                    (leave.vel, leave.spin) =
                        lock_rolling(leave.vel, leave.spin, n, radius, ball.inertia_factor);
                }
                contacts.push(Contact::new(&state, &at, &leave));
                state = leave;
                break;
            }
            contacts.push(Contact::new(&state, &at, &out));
            points.push(GroundPoint::new(&out, Phase::Hop));
            if contacts.len() >= self.max_impacts {
                return GroundResult::new(contacts, points, out, Termination::MaxImpacts);
            }
            match self.hop(&out, t_end, &mut points) {
                Ok(arrival) => state = arrival,
                Err(stopped) => {
                    points.push(GroundPoint::new(&stopped, Phase::Hop));
                    return GroundResult::new(contacts, points, stopped, Termination::MaxTime);
                }
            }
        }

        // Slide and roll.
        let mut next_sample = state.time;
        let mut last_phase = None;
        loop {
            if state.time >= t_end - TIME_EPS {
                let phase = last_phase.unwrap_or(Phase::Roll);
                points.push(GroundPoint::new(&state, phase));
                return GroundResult::new(contacts, points, state, Termination::MaxTime);
            }
            let at = self.sample_contact(state.pos, radius);
            let h = self.dt.min(t_end - state.time);
            let (next, phase) = contact_step(&state, &at, ball, h, self.stop_speed);
            if state.time >= next_sample - TIME_EPS || last_phase != Some(phase) {
                points.push(GroundPoint::new(&state, phase));
                next_sample = state.time + SAMPLE_INTERVAL;
            }
            last_phase = Some(phase);
            state = next;
            if phase == Phase::Roll && state.vel == DVec3::ZERO {
                points.push(GroundPoint::new(&state, Phase::Roll));
                return GroundResult::new(contacts, points, state, Termination::Rest);
            }
        }
    }
}

/// One impact of the ball with the ground: the post-impact state (same time
/// and position, new velocity and spin).
///
/// A rigid-sphere point-contact impulse against Penner's (2002) effective
/// contact plane, in vector form. With `n` the surface normal, `t̂` the
/// direction of the incoming tangential velocity, `R` the radius and `k` the
/// inertia factor:
///
/// 1. Tilt: β = `compliance` · θc (Penner eq. 8), capped at the incidence
///    angle; effective normal `n' = n cos β − t̂ sin β`, leaning back toward
///    the incoming ball.
/// 2. Normal impulse (per unit mass) `j_n = (1 + e)·|v·n'|`, with `e` from
///    [`Surface::restitution`] at the normal speed `|v·n'|`.
/// 3. Contact-point slip `u = (v − (v·n')n') + ω × (−R n')`. Friction
///    `j_t = min(|u|·k/(1+k), μ·j_n)` acts along `−û`: the ball either grips
///    (leaves with zero slip) or slides throughout (Coulomb limit).
/// 4. `v += j_n n' − j_t û`, `ω += (−R n') × (−j_t û) / (k R²)`.
/// 5. Pivot friction on the spin about the true normal `n` (not `n'`, so it
///    never removes spin along the heading):
///    `Δω_n = −sign(ω_n) · min(|ω_n|, pivot · j_n / (k R²))`.
///
/// No impulse is applied if the ball is not moving into the ground
/// (`v·n ≥ 0`, grazing). Whether the ball then hops or stays on the ground
/// is up to the caller ([`GroundPhysics`] uses the 5 mm hop rule).
#[must_use]
pub fn impact(state: &BallState, at: &GroundSample, ball: &BallModel) -> BallState {
    let n = at.normal.normalize();
    let v = state.vel;
    let v_dot_n = v.dot(n);
    if v_dot_n >= 0.0 {
        return *state;
    }
    let surface = &at.surface;
    let radius = ball.radius_m();
    let k = ball.inertia_factor;

    // Penner's effective contact plane.
    let v_t_vec = v - v_dot_n * n;
    let v_t = v_t_vec.length();
    let beta = surface.contact_angle(v_t, -v_dot_n);
    let t_hat = if v_t > 0.0 {
        v_t_vec / v_t
    } else {
        DVec3::ZERO
    };
    let n_eff = n * beta.cos() - t_hat * beta.sin();
    let v_n_eff = v.dot(n_eff);
    if v_n_eff >= 0.0 {
        return *state;
    }

    let e = surface.restitution.coefficient_at(v_n_eff);
    let j_n = (1.0 + e) * -v_n_eff;

    let r = -radius * n_eff;
    let slip = (v - v_n_eff * n_eff) + state.spin.cross(r);
    let slip_speed = slip.length();
    let j_t = (slip_speed * k / (1.0 + k)).min(surface.mu_slide * j_n);
    let j_t_vec = if slip_speed > 0.0 {
        slip * (-j_t / slip_speed)
    } else {
        DVec3::ZERO
    };

    let inertia = k * radius * radius; // per unit mass
    let vel = v + j_n * n_eff + j_t_vec;
    let mut spin = state.spin + r.cross(j_t_vec) / inertia;

    let w_n = spin.dot(n);
    let dw = w_n.abs().min(surface.pivot * j_n / inertia);
    spin -= n * (w_n.signum() * dw);

    BallState {
        vel,
        spin,
        ..*state
    }
}

/// Position, velocity and spin of a ball on the ground, for RK4.
#[derive(Debug, Clone, Copy)]
struct Motion {
    pos: DVec3,
    vel: DVec3,
    spin: DVec3,
}

impl Integrable for Motion {
    fn add_scaled(&self, rate: &Self, h: f64) -> Self {
        Self {
            pos: self.pos + rate.pos * h,
            vel: self.vel + rate.vel * h,
            spin: self.spin + rate.spin * h,
        }
    }
}

/// Contact-point slip velocity of a ball on ground with unit normal `n`.
fn slip(vel: DVec3, spin: DVec3, n: DVec3, radius: f64) -> DVec3 {
    (vel - vel.dot(n) * n) + spin.cross(-radius * n)
}

/// Remove the contact slip with the friction impulse that brings the ball
/// to rolling (angular momentum about the contact point is conserved).
fn lock_rolling(vel: DVec3, spin: DVec3, n: DVec3, radius: f64, k: f64) -> (DVec3, DVec3) {
    let u = slip(vel, spin, n, radius);
    (
        vel - u * (k / (1.0 + k)),
        spin + n.cross(u) / ((1.0 + k) * radius),
    )
}

/// Decay the spin component about `n` from `start` to `end` without letting
/// it cross zero (pivot friction saturates).
fn saturate_pivot(start: DVec3, end: DVec3, n: DVec3) -> DVec3 {
    let w0 = start.dot(n);
    let w1 = end.dot(n);
    if w0 == 0.0 || w0.signum() != w1.signum() {
        end - w1 * n
    } else {
        end
    }
}

/// One ground step of at most `dt` on the ground sample `at`; see
/// [`GroundPhysics::step`].
fn contact_step(
    state: &BallState,
    at: &GroundSample,
    ball: &BallModel,
    dt: f64,
    stop_speed: f64,
) -> (BallState, Phase) {
    let n = at.normal.normalize();
    let surface = &at.surface;
    let radius = ball.radius_m();
    let k = ball.inertia_factor;
    let inertia = k * radius * radius; // per unit mass

    // Normal load and tangential gravity per unit mass.
    let g_n = G * n.y;
    let g_t = GRAVITY - GRAVITY.dot(n) * n;

    let vel = state.vel - state.vel.dot(n) * n;
    let spin = state.spin;
    let w_n = spin.dot(n);
    let pivot_accel = -w_n.signum() * surface.pivot * g_n / inertia * n;

    let u = slip(vel, spin, n, radius);
    let slip_speed = u.length();

    if slip_speed > SLIP_EPS {
        let u_hat = u / slip_speed;
        let friction = surface.mu_slide * g_n;
        let slip_decay = (1.0 + k) / k * friction;
        let to_roll = if slip_decay > 0.0 {
            slip_speed / slip_decay
        } else {
            f64::INFINITY
        };
        let h = dt.min(to_roll);
        let accel = -friction * u_hat + g_t;
        let spin_accel = friction * n.cross(u_hat) / (k * radius) + pivot_accel;
        let m = rk4::step_state(
            Motion {
                pos: state.pos,
                vel,
                spin,
            },
            state.time,
            h,
            |s, _| Motion {
                pos: s.vel,
                vel: accel,
                spin: spin_accel,
            },
        );
        let mut spin_out = saturate_pivot(spin, m.spin, n);
        let mut vel_out = m.vel;
        if to_roll <= dt {
            (vel_out, spin_out) = lock_rolling(vel_out, spin_out, n, radius, k);
        }
        let next = BallState {
            time: state.time + h,
            pos: m.pos,
            vel: vel_out,
            spin: spin_out,
        };
        return (next, Phase::Slide);
    }

    // Rolling: lock the tangential spin to the velocity.
    let (vel, spin) = lock_rolling(vel, spin, n, radius, k);
    let speed = vel.length();
    let decel = |v: f64| surface.rolling_friction(v) * g_n / (1.0 + k);
    let a0 = decel(speed);
    if speed <= stop_speed.max(a0 * dt) {
        let t_stop = if a0 > 0.0 { speed / a0 } else { 0.0 };
        let spin_n = w_n * n;
        let spin_rest = saturate_pivot(spin_n, spin_n + pivot_accel * t_stop, n);
        let rest = BallState {
            time: state.time + t_stop,
            pos: state.pos + vel * (0.5 * t_stop),
            vel: DVec3::ZERO,
            spin: spin_rest,
        };
        return (rest, Phase::Roll);
    }

    let slope = g_t / (1.0 + k);
    let m = rk4::step_state(
        Motion {
            pos: state.pos,
            vel,
            spin,
        },
        state.time,
        dt,
        |s, _| {
            let v = s.vel.length();
            let dv = -decel(v) * s.vel / v + slope;
            Motion {
                pos: s.vel,
                vel: dv,
                spin: n.cross(dv) / radius + pivot_accel,
            }
        },
    );
    let spin_n = saturate_pivot(spin, m.spin, n).dot(n) * n;
    let next = BallState {
        time: state.time + dt,
        pos: m.pos,
        vel: m.vel,
        spin: n.cross(m.vel) / radius + spin_n,
    };
    (next, Phase::Roll)
}

#[cfg(test)]
mod regression_data;
#[cfg(test)]
mod tests;
