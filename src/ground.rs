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

/// Contact-point slip speed, or rolling speed, (m/s) treated as zero.
/// Slides end exactly at zero slip and rolls exactly at zero speed; this
/// only absorbs rounding.
const SLIP_EPS: f64 = 1e-9;

/// Spacing (s) of the sampled [`GroundPoint`]s.
const SAMPLE_INTERVAL: f64 = 0.01;

/// Time tolerance (s) for reaching sample times and the time limit.
const TIME_EPS: f64 = 1e-9;

/// Smallest slip speed (m/s) at which a slide on a slope switches to
/// finishing along a fixed slip direction, whatever
/// [`GroundPhysics::roll_threshold`] is.
const MIN_SWITCH_SLIP: f64 = 1e-6;

/// Normal speed, relative to the ball's speed, below which an impact counts
/// as grazing: rounding of a velocity that lies in the ground plane.
const GRAZING: f64 = 1e-12;

/// Residual accepted by the crossing solves: metres for ground crossings,
/// m/s for the slide-to-roll event.
const ROOT_TOL: f64 = 1e-12;

/// Iteration cap of the crossing solves; they converge in a handful.
const ROOT_MAX_ITER: usize = 200;

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
    /// Unit outward surface normal at (x, z), with `normal.y > 0`. For a
    /// heightfield `y = h(x, z)` this is `(−∂h/∂x, 1, −∂h/∂z)` normalised.
    /// [`GroundPhysics`] uses it for the impact plane, the normal load and
    /// the downhill pull, so it should match the slope of `height`.
    pub normal: DVec3,
    /// Material at this location.
    pub surface: Surface,
}

impl GroundSample {
    /// Whether a ball at rest here stays put: the slope's pull, `g·sin θ`,
    /// does not exceed the rolling resistance at zero speed,
    /// `mu_roll·g·cos θ`, i.e. `tan θ ≤ mu_roll` (θ the slope angle). Always
    /// true on level ground.
    #[must_use]
    pub fn holds_at_rest(&self) -> bool {
        holds_at_rest(self.normal.normalize(), &self.surface)
    }
}

/// Whether rolling resistance at zero speed holds a ball against the slope
/// of ground with unit normal `n`.
fn holds_at_rest(n: DVec3, surface: &Surface) -> bool {
    tangential_gravity(n).length() <= surface.mu_roll * G * n.y
}

/// Layer 1 of the ground API: "what is the ground like here?"
///
/// A heightfield point query: one ground height, normal and surface per
/// horizontal position. dimple never needs to know how the ground is
/// stored: an analytic plane, a grid, a mesh lookup or an engine query all
/// fit behind `sample`.
///
/// [`GroundPhysics`] tracks the ball's centre and samples the terrain at
/// the centre's (x, z): the ball is on the ground when `pos.y` equals the
/// sampled [`height`](GroundSample::height), and the sampled
/// [`normal`](GroundSample::normal) gives the contact plane, the normal load
/// and the downhill pull. It samples every step, so surfaces and slopes can
/// change mid-roll. Each ground step is taken in the tangent plane under
/// the ball, which is then put back on the heightfield with its velocity
/// projected onto the new tangent plane. On a plane this is exact; on
/// curved ground it is a first-order approximation:
///
/// - The ball always follows the heightfield. Over a crest where
///   `v²·κ > g·n_y` (κ the curvature along the path) a real ball would
///   leave the ground; that launch is not modelled.
/// - Into a concave transition (the foot of a slope) the velocity component
///   along the new normal is removed, as in an inelastic contact.
///
/// One height per (x, z) means no walls, lips or cup; a step in `height`
/// is crossed as if the ball were lifted onto it. A ball can roll without
/// slipping only where `tan θ ≤ ((1+k)/k)·mu_slide` (54° for the presets);
/// steeper ground is outside the model.
///
/// See [`FnTerrain`] for examples.
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
/// Per-point surfaces on level ground: a green from 150 m to 170 m
/// downrange, fairway everywhere else.
///
/// ```
/// use dimple::ground::{FnTerrain, GroundPhysics, GroundSample, Surface};
/// use glam::DVec3;
///
/// let range = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| GroundSample {
///     height: 0.0,
///     normal: DVec3::Y,
///     surface: if (150.0..170.0).contains(&x) { Surface::GREEN } else { Surface::FAIRWAY },
/// }));
/// # let _ = range;
/// ```
///
/// A plane: fairway rising 3° downrange and falling 2° to the right. For
/// `height = a·x + b·z` the normal is `(−a, 1, −b)` normalised.
///
/// ```
/// use dimple::aero::BallModel;
/// use dimple::environment::Environment;
/// use dimple::ground::{FnTerrain, GroundPhysics, GroundSample, Surface};
/// use dimple::trajectory::{ShotInput, simulate_shot};
/// use glam::DVec3;
///
/// let (a, b) = (3.0_f64.to_radians().tan(), -(2.0_f64.to_radians().tan()));
/// let normal = DVec3::new(-a, 1.0, -b).normalize();
/// let hill = GroundPhysics::new(FnTerrain(move |x: f64, z: f64| GroundSample {
///     height: a * x + b * z,
///     normal,
///     surface: Surface::FAIRWAY,
/// }));
/// let shot = simulate_shot(
///     &ShotInput::seven_iron(),
///     &BallModel::TOUR,
///     &Environment::SEA_LEVEL,
///     &hill,
/// );
/// let rest = shot.final_pos;
/// assert!((rest.y - (a * rest.x + b * rest.z)).abs() < 1e-9); // on the plane
/// assert!(rest.z > 0.0); // drifted right, down the side slope
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
/// [`min_hop_height`](Self::min_hop_height) and, on a slope, the hops die
/// out rather than being sustained by the slope (see `min_hop_height`), then
/// slides and rolls ([`step`](Self::step)) until it stops. The tunables are public fields; start
/// from [`GroundPhysics::new`] and override with struct update syntax.
#[derive(Debug, Clone, Copy)]
pub struct GroundPhysics<T: Terrain> {
    /// Where the ball is and what it lands on.
    pub terrain: T,
    /// Integration step (s) for hops, slides and rolls. Default 1 ms.
    pub dt: f64,
    /// A rolling ball slower than this (m/s) stops, unless the slope is
    /// speeding it up. It also stops once its speed is below one step's
    /// deceleration. Either way it finishes at constant deceleration, so
    /// this barely moves the stopping point. Where the slope is too steep
    /// to hold it ([`GroundSample::holds_at_rest`]) the ball starts rolling
    /// downhill again from there. Default 0.01 m/s.
    pub stop_speed: f64,
    /// Contact-point slip speed (m/s) at or below which a ball that has
    /// stopped bouncing starts out rolling instead of sliding. The residual
    /// slip is removed by the friction impulse that would have finished the
    /// slide, so the rolling velocity is the same either way. Slides are
    /// otherwise stepped to zero slip; on sloped ground this is also the
    /// slip at which a slide stops following its turning slip direction and
    /// finishes along a fixed one (see [`step`](Self::step)), at no less
    /// than 1e-6 m/s whatever this is set to: below that the slip direction
    /// could flip within one step. Default 0.01 m/s.
    pub roll_threshold: f64,
    /// The ball keeps bouncing while its rebound would rise at least this
    /// high (m) above the ground, measured along the normal:
    /// `(v·n)² / (2·g·n_y)`. Penner's (2002) bounce-to-roll transition,
    /// 5 mm. Penner tested 1–10 mm and found only a minor effect on the run.
    ///
    /// On a slope the ball also needs the hops to die out: the impact must
    /// remove at least the energy of a hop this high (`g·min_hop_height` per
    /// unit mass) more than gravity hands back over the next hop down the
    /// slope. Otherwise the ball goes to the ground at once, as below the
    /// height. Without this a downhill hop train settles at a fixed hop
    /// height just above the transition and never rolls (see
    /// `docs/AERO.md`). Level ground uses the height alone.
    pub min_hop_height: f64,
    /// Safety bound on impacts. Every hop rises at least
    /// [`min_hop_height`](Self::min_hop_height) and the hops die out, so
    /// real shots stop bouncing far below this; it guards degenerate surfaces
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
    /// [`dt`](Self::dt), on the terrain under its centre.
    ///
    /// Returns the new state and the phase the ball was in over the step.
    /// With `n` the sampled normal, normal load `g·n_y` and tangential
    /// gravity `g_t = g − (g·n)n` (zero on level ground), RK4 within each
    /// smooth phase:
    ///
    /// - **Slide** (the contact point slips): friction `−μ_slide·g·n_y·û`,
    ///   its torque `(−R n) × F`, and `g_t`. Where `g_t` is zero the slip
    ///   direction is fixed and the slip speed falls linearly at
    ///   `((1+k)/k)·μ_slide·g·n_y`, so a slide that would end within the
    ///   step stops exactly at the slide-to-roll time (the step is
    ///   shortened) and the ball leaves rolling. On a slope `g_t` turns the
    ///   slip direction, so û is re-evaluated at every RK4 stage; the time
    ///   the slip falls to [`roll_threshold`](Self::roll_threshold) is solved
    ///   within the step, and from there the slide finishes along that slip
    ///   direction, ending at zero slip as on level ground.
    /// - **Roll**: spin tangential to the ground is locked to the velocity;
    ///   deceleration `μ_eff·g·n_y/(1+k)` against the velocity, with
    ///   `μ_eff` from [`Surface::rolling_friction`], plus `g_t/(1+k)`. A ball
    ///   that is slowing and slower than one step's deceleration (or
    ///   [`stop_speed`](Self::stop_speed)) stops: it finishes at constant
    ///   deceleration, the returned velocity is exactly zero, and the time
    ///   advances by the stopping time only. Whether it then stays at rest
    ///   is [`GroundSample::holds_at_rest`]; if not, the next step rolls it
    ///   downhill.
    /// - In both, pivot friction decays the spin about the normal at
    ///   `pivot·g·n_y/(k R²)`, saturating at zero.
    ///
    /// The ball ends on the terrain at its new (x, z), its velocity in the
    /// tangent plane there and, when rolling, its tangential spin locked to
    /// that velocity. Unlike [`settle`](GroundModel::settle), this does not
    /// apply [`max_time`](Self::max_time).
    #[must_use]
    pub fn step(&self, state: &BallState, ball: &BallModel) -> (BallState, Phase) {
        let at = self.sample(state.pos);
        let (next, phase, rolling) = contact_step(
            state,
            &at,
            ball,
            self.dt,
            self.stop_speed,
            self.roll_threshold,
        );
        let under = self.sample(next.pos);
        (
            onto_surface(&next, &under, at.normal, rolling, ball.radius_m()),
            phase,
        )
    }

    /// Terrain under a ball centred at `pos`.
    fn sample(&self, pos: DVec3) -> GroundSample {
        self.terrain.sample(pos.x, pos.z)
    }

    /// Height of a ball centred at `pos` above the terrain under it.
    fn clearance(&self, pos: DVec3) -> f64 {
        pos.y - self.sample(pos).height
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
            if let Some(arrival) = self.hop_arrival(s, t, h, self.clearance(b.pos), start.spin) {
                return Ok(arrival);
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

    /// First ground crossing within a gravity-only step of length `h` from
    /// `a` (at time `t`), if there is one; `end` is the clearance at the end
    /// of the integrated step.
    ///
    /// Within the step the motion is the exact parabola RK4 reproduces. The
    /// clearance `y(τ) − height(x(τ), z(τ))` is checked at
    /// [`CROSSING_SUBSTEPS`] equal sub-steps for the first at or below the
    /// ground, and the crossing is the root inside that sub-step, solved by
    /// [`find_root`] from the drop time onto the height at the end of the
    /// step (exact on level ground). The arrival is placed exactly on the
    /// ground.
    fn hop_arrival(
        &self,
        a: rk4::State,
        t: f64,
        h: f64,
        end: f64,
        spin: DVec3,
    ) -> Option<BallState> {
        let at = |tau: f64| a.pos + a.vel * tau + 0.5 * GRAVITY * tau * tau;
        let clearance = |tau: f64| self.clearance(at(tau));
        let (lo, (hi, f_hi)) = first_crossing(clearance, h, clearance(0.0), end)?;
        let f_hi = if hi == h { clearance(h).min(0.0) } else { f_hi };
        let end_height = self.sample(at(h)).height;
        let guess = drop_time(a.pos.y - end_height, a.vel.y).clamp(lo.0, hi);
        let tau = find_root(clearance, lo, (hi, f_hi), guess);
        let mut pos = at(tau);
        pos.y = self.sample(pos).height;
        Some(BallState {
            time: t + tau,
            pos,
            vel: a.vel + GRAVITY * tau,
            spin,
        })
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

/// Number of equal sub-steps at which a flight or hop step is checked for
/// its first ground crossing.
const CROSSING_SUBSTEPS: u32 = 8;

/// The first of [`CROSSING_SUBSTEPS`] equal sub-steps of `[0, len]` that
/// ends at or below the ground: its bounds, each with the clearance `f`
/// there, or `None` if every sub-step ends above it. `f0` is the clearance
/// at 0 and `f_end` stands in for `f(len)`.
fn first_crossing(
    f: impl Fn(f64) -> f64,
    len: f64,
    f0: f64,
    f_end: f64,
) -> Option<((f64, f64), (f64, f64))> {
    let mut lo = (0.0, f0);
    for i in 1..=CROSSING_SUBSTEPS {
        let hi = if i == CROSSING_SUBSTEPS {
            (len, f_end)
        } else {
            let x = len * f64::from(i) / f64::from(CROSSING_SUBSTEPS);
            (x, f(x))
        };
        if hi.1 <= 0.0 {
            return Some((lo, hi));
        }
        lo = hi;
    }
    None
}

/// A root of `f` on `[lo, hi]` given `f(lo) ≥ 0 ≥ f(hi)` (each bound passed
/// with its value), starting from `guess`.
///
/// The guess is accepted outright if `|f(guess)| ≤ ROOT_TOL`; otherwise the
/// bracket is narrowed by Illinois-modified regula falsi, bisecting whenever
/// an iterate leaves the bracket. If the bracket stops shrinking (a jump in
/// `f`, such as a step in the terrain) the bound with `f ≤ 0` is returned.
fn find_root(f: impl Fn(f64) -> f64, lo: (f64, f64), hi: (f64, f64), guess: f64) -> f64 {
    let (mut lo, mut f_lo) = lo;
    let (mut hi, mut f_hi) = hi;
    let mut x = guess;
    let mut last_side = 0_i8;
    for _ in 0..ROOT_MAX_ITER {
        if !(lo..=hi).contains(&x) {
            x = 0.5 * (lo + hi);
        }
        let fx = f(x);
        if fx.abs() <= ROOT_TOL {
            return x;
        }
        if fx > 0.0 {
            lo = x;
            f_lo = fx;
            if last_side == 1 {
                f_hi *= 0.5;
            }
            last_side = 1;
        } else {
            hi = x;
            f_hi = fx;
            if last_side == -1 {
                f_lo *= 0.5;
            }
            last_side = -1;
        }
        let mid = 0.5 * (lo + hi);
        if !(lo < mid && mid < hi) {
            break;
        }
        x = if f_lo > f_hi {
            lo + (hi - lo) * (f_lo / (f_lo - f_hi))
        } else {
            mid
        };
        if !(lo < x && x < hi) {
            x = mid;
        }
    }
    hi
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
    /// The default treats the segment as straight, with time, velocity and
    /// spin varying linearly along it. It checks the clearance
    /// `y − height(x, z)` at 8 equal sub-steps along the segment for the
    /// first at or below the ground, and solves for the crossing inside that
    /// sub-step: Illinois-modified regula falsi from the linear
    /// interpolation of the end clearances (exact on level ground),
    /// bisecting when needed, to 1e-12 m. It returns the state there with
    /// the position placed exactly on the ground, so elevated greens,
    /// downhill landings, curved terrain and ridges thinner than the
    /// segment all land where the segment first crosses the surface. A
    /// segment that starts below the ground does not land, nor does one
    /// above the ground at every sub-step: a feature narrower than a
    /// sub-step (about 1 cm of a 1 ms flight step at 75 m/s) can be missed.
    /// An engine can override this with its own sweep test.
    #[must_use]
    fn landing(&self, a: &BallState, b: &BallState) -> Option<BallState> {
        let fa = a.pos.y - self.height(a.pos.x, a.pos.z);
        if fa < 0.0 {
            return None;
        }
        let fb = b.pos.y - self.height(b.pos.x, b.pos.z);
        let clearance = |s: f64| {
            let p = a.pos.lerp(b.pos, s);
            p.y - self.height(p.x, p.z)
        };
        let (lo, hi) = first_crossing(clearance, 1.0, fa, fb)?;
        let guess = if fa > fb { fa / (fa - fb) } else { 0.0 };
        let frac = find_root(clearance, lo, hi, guess.clamp(lo.0, hi.0));
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
        let mut at = loop {
            let at = self.sample(state.pos);
            let n = at.normal.normalize();
            let out = impact(&state, &at, ball);
            let v_n = out.vel.dot(n);
            let hops = v_n > 0.0
                && v_n * v_n >= 2.0 * G * n.y * self.min_hop_height
                && hops_die_out(&state, &out, n, ball, self.min_hop_height);
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
                break at;
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
        };

        // Slide and roll, each step on the terrain under the ball.
        let mut next_sample = state.time;
        let mut last_phase = None;
        loop {
            if state.time >= t_end - TIME_EPS {
                let phase = last_phase.unwrap_or(Phase::Roll);
                points.push(GroundPoint::new(&state, phase));
                return GroundResult::new(contacts, points, state, Termination::MaxTime);
            }
            let h = self.dt.min(t_end - state.time);
            let (next, phase, rolling) =
                contact_step(&state, &at, ball, h, self.stop_speed, self.roll_threshold);
            if state.time >= next_sample - TIME_EPS || last_phase != Some(phase) {
                points.push(GroundPoint::new(&state, phase));
                next_sample = state.time + SAMPLE_INTERVAL;
            }
            last_phase = Some(phase);
            let stopped = phase == Phase::Roll && next.vel == DVec3::ZERO;
            let under = self.sample(next.pos);
            let held = stopped && under.holds_at_rest();
            state = onto_surface(&next, &under, at.normal, rolling, radius);
            at = under;
            if held {
                points.push(GroundPoint::new(&state, Phase::Roll));
                return GroundResult::new(contacts, points, state, Termination::Rest);
            }
        }
    }
}

/// Whether the hops that follow a rebound off ground with unit normal `n`
/// die out, so that the ball may keep bouncing; otherwise it goes to the
/// ground at once, as below [`GroundPhysics::min_hop_height`].
///
/// Always true on level ground, where every hop loses energy. On a slope
/// gravity hands energy back over each hop: the potential energy released
/// over the hop on the local plane, `ΔU = −g·Δy` with flight time
/// `T = 2(v·n)/(g·n_y)` and `Δy = v_y·T − g·T²/2`, positive downhill. The
/// impact `arrival → rebound` removed
/// `ΔE = E(arrival) − E(rebound)`, `E = |v|²/2 + k·R²·|ω|²/2` (per unit
/// mass). The train dies out if the slope returns nothing (`ΔU ≤ 0`) or the
/// impact removes at least the energy of a `min_hop_height` hop more than
/// the slope returns (`ΔE − ΔU ≥ g·min_hop_height`). Below that, gravity
/// re-supplies what each impact takes and the hops would settle at a fixed
/// size above the 5 mm transition instead of shrinking through it.
fn hops_die_out(
    arrival: &BallState,
    rebound: &BallState,
    n: DVec3,
    ball: &BallModel,
    min_hop_height: f64,
) -> bool {
    if tangential_gravity(n) == DVec3::ZERO {
        return true;
    }
    let t = 2.0 * rebound.vel.dot(n) / (G * n.y);
    let released = -G * (rebound.vel.y * t - 0.5 * G * t * t);
    if released <= 0.0 {
        return true;
    }
    let spin_inertia = ball.inertia_factor * ball.radius_m() * ball.radius_m();
    let energy =
        |s: &BallState| 0.5 * s.vel.length_squared() + 0.5 * spin_inertia * s.spin.length_squared();
    let removed = energy(arrival) - energy(rebound);
    removed - released >= G * min_hop_height
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
/// (`v·n ≥ 0`, grazing; a normal speed within rounding of zero,
/// `|v·n| ≤ 1e-12·|v|`, counts as grazing, so a ball moving in a sloped
/// ground plane is not struck by its own rounding error). Whether the ball then hops or stays on the ground
/// is up to the caller ([`GroundPhysics`] uses the 5 mm hop rule; see
/// [`GroundPhysics::min_hop_height`]).
#[must_use]
pub fn impact(state: &BallState, at: &GroundSample, ball: &BallModel) -> BallState {
    let n = at.normal.normalize();
    let v = state.vel;
    let v_dot_n = v.dot(n);
    if v_dot_n >= -GRAZING * v.length() {
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

/// Tangential gravity `g − (g·n)n` on ground with unit normal `n`: the
/// downhill pull (m/s²), zero on level ground.
fn tangential_gravity(n: DVec3) -> DVec3 {
    GRAVITY - GRAVITY.dot(n) * n
}

/// Put a ball that has just stepped across the ground back onto it.
///
/// `at` is the terrain under the ball's new (x, z) and `prev_normal` the
/// normal the step was taken on. The centre goes to the terrain height and
/// the velocity into the tangent plane at `at`; a ball that ended the step
/// `rolling` (a roll, or a slide that reached zero slip) also gets its
/// tangential spin re-locked to that velocity, keeping its spin about the
/// normal. Where the normal and height already match (level ground, a
/// plane) the state is returned unchanged.
fn onto_surface(
    state: &BallState,
    at: &GroundSample,
    prev_normal: DVec3,
    rolling: bool,
    radius: f64,
) -> BallState {
    if at.normal == prev_normal && state.pos.y == at.height {
        return *state;
    }
    let n = at.normal.normalize();
    let mut out = *state;
    out.pos.y = at.height;
    out.vel = state.vel - state.vel.dot(n) * n;
    if rolling {
        let w_n = state.spin.dot(prev_normal.normalize());
        out.spin = n.cross(out.vel) / radius + w_n * n;
    }
    out
}

/// The forces on a ball on one ground sample, per unit mass.
struct Patch {
    /// Unit normal.
    n: DVec3,
    /// Tangential gravity (m/s²).
    g_t: DVec3,
    /// Sliding friction magnitude `μ_slide·g·n_y` (m/s²).
    friction: f64,
    /// Ball radius (m).
    radius: f64,
    /// Inertia factor k.
    k: f64,
    /// Pivot-friction angular acceleration about `n` for this step.
    pivot_accel: DVec3,
}

impl Patch {
    /// Rate of change of a sliding ball with slip direction `u_hat`.
    fn slide_rate(&self, s: &Motion, u_hat: DVec3) -> Motion {
        Motion {
            pos: s.vel,
            vel: -self.friction * u_hat + self.g_t,
            spin: self.friction * self.n.cross(u_hat) / (self.k * self.radius) + self.pivot_accel,
        }
    }

    /// One RK4 slide step of length `h` along the fixed slip direction
    /// `u_hat`.
    fn slide_fixed(&self, s: Motion, t: f64, h: f64, u_hat: DVec3) -> Motion {
        let accel = -self.friction * u_hat + self.g_t;
        let spin_accel =
            self.friction * self.n.cross(u_hat) / (self.k * self.radius) + self.pivot_accel;
        rk4::step_state(s, t, h, |s, _| Motion {
            pos: s.vel,
            vel: accel,
            spin: spin_accel,
        })
    }

    /// One RK4 slide step of length `h`, with the slip direction evaluated
    /// at every stage (`fallback` where the slip vanishes).
    fn slide_turning(&self, s: Motion, t: f64, h: f64, fallback: DVec3) -> Motion {
        rk4::step_state(s, t, h, |s, _| {
            let u = slip(s.vel, s.spin, self.n, self.radius);
            let len = u.length();
            let u_hat = if len > SLIP_EPS { u / len } else { fallback };
            self.slide_rate(s, u_hat)
        })
    }

    /// Contact slip speed of `s`.
    fn slip_speed(&self, s: &Motion) -> f64 {
        slip(s.vel, s.spin, self.n, self.radius).length()
    }
}

/// One ground step of at most `dt` on the ground sample `at`; see
/// [`GroundPhysics::step`]. Returns the new state, the phase over the step
/// and whether the ball ends it rolling. The ball is left in the tangent
/// plane at `at`; putting it back on the terrain is up to the caller.
fn contact_step(
    state: &BallState,
    at: &GroundSample,
    ball: &BallModel,
    dt: f64,
    stop_speed: f64,
    roll_threshold: f64,
) -> (BallState, Phase, bool) {
    let n = at.normal.normalize();
    let surface = &at.surface;
    let radius = ball.radius_m();
    let k = ball.inertia_factor;
    let inertia = k * radius * radius; // per unit mass

    // Normal load and tangential gravity per unit mass.
    let g_n = G * n.y;
    let g_t = tangential_gravity(n);

    let vel = state.vel - state.vel.dot(n) * n;
    let spin = state.spin;
    let w_n = spin.dot(n);
    // Pivot friction opposes the spin about the normal; with none (of
    // either sign of zero) there is none.
    let pivot_accel = if w_n == 0.0 {
        DVec3::ZERO
    } else {
        -w_n.signum() * surface.pivot * g_n / inertia * n
    };

    let slip_speed = slip(vel, spin, n, radius).length();

    if slip_speed > SLIP_EPS {
        let patch = Patch {
            n,
            g_t,
            friction: surface.mu_slide * g_n,
            radius,
            k,
            pivot_accel,
        };
        let start = Motion {
            pos: state.pos,
            vel,
            spin,
        };
        let (m, h, rolls) = slide(&patch, start, state.time, dt, roll_threshold);
        let mut spin_out = saturate_pivot(spin, m.spin, n);
        let mut vel_out = m.vel;
        if rolls {
            (vel_out, spin_out) = lock_rolling(vel_out, spin_out, n, radius, k);
        }
        let next = BallState {
            time: state.time + h,
            pos: m.pos,
            vel: vel_out,
            spin: spin_out,
        };
        return (next, Phase::Slide, rolls);
    }

    // Rolling: lock the tangential spin to the velocity.
    let (vel, spin) = lock_rolling(vel, spin, n, radius, k);
    let speed = vel.length();
    let decel = |v: f64| surface.rolling_friction(v) * g_n / (1.0 + k);
    let slope = g_t / (1.0 + k);
    // Direction of motion; from rest (or a rounding-level speed), the way
    // the slope pulls.
    let downhill = slope.normalize_or_zero();
    let heading = if speed > SLIP_EPS {
        vel / speed
    } else {
        downhill
    };
    // Net deceleration along the motion: rolling resistance less the
    // slope's pull along it.
    let a0 = decel(speed) - slope.dot(heading);
    // At (rounding-level) zero speed the ball stops exactly when the slope
    // holds it, the same test the caller uses to end the roll.
    let stops = if speed <= SLIP_EPS {
        holds_at_rest(n, surface)
    } else {
        a0 >= 0.0 && speed <= stop_speed.max(a0 * dt)
    };
    if stops {
        let t_stop = if a0 > 0.0 { speed / a0 } else { 0.0 };
        // The slope's pull across the motion bends the last stretch.
        let across = slope - slope.dot(heading) * heading;
        let spin_n = w_n * n;
        let spin_rest = saturate_pivot(spin_n, spin_n + pivot_accel * t_stop, n);
        let rest = BallState {
            time: state.time + t_stop,
            pos: state.pos + vel * (0.5 * t_stop) + across * (0.5 * t_stop * t_stop),
            vel: DVec3::ZERO,
            spin: spin_rest,
        };
        return (rest, Phase::Roll, true);
    }

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
            let resistance = if v > SLIP_EPS {
                -decel(v) * s.vel / v
            } else {
                -decel(0.0) * downhill
            };
            let dv = resistance + slope;
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
    (next, Phase::Roll, true)
}

/// Slide from `start` (slip above `SLIP_EPS`) for at most `dt`: the state
/// reached, the time taken, and whether the slip reached zero (the caller
/// then locks the ball to rolling).
fn slide(
    patch: &Patch,
    start: Motion,
    t: f64,
    dt: f64,
    roll_threshold: f64,
) -> (Motion, f64, bool) {
    // d|u|/dt = −c + g_t·û: the slip shrinks at `c` less the slope's pull
    // along it. The switch to finishing along a fixed slip direction is kept
    // above zero slip, where û would flip within the RK4 stages.
    let roll_threshold = roll_threshold.max(MIN_SWITCH_SLIP);
    let c = (1.0 + patch.k) / patch.k * patch.friction;
    let finish = |s: Motion, t: f64, left: f64| -> (Motion, f64, bool) {
        // Slide to zero slip along the current slip direction: exact when
        // that direction cannot turn (level ground, or slip along the fall
        // line of a plane).
        let u = slip(s.vel, s.spin, patch.n, patch.radius);
        let len = u.length();
        if len <= SLIP_EPS {
            return (s, 0.0, true);
        }
        let u_hat = u / len;
        let rate = c - patch.g_t.dot(u_hat);
        let to_roll = if rate > 0.0 {
            len / rate
        } else {
            f64::INFINITY
        };
        let h = left.min(to_roll);
        (patch.slide_fixed(s, t, h, u_hat), h, to_roll <= left)
    };

    if patch.g_t == DVec3::ZERO {
        // Level ground: the slip direction is fixed and |u| falls linearly.
        return finish(start, t, dt);
    }

    let u0 = slip(start.vel, start.spin, patch.n, patch.radius);
    let s0 = u0.length();
    let u_hat0 = u0 / s0;
    if s0 <= roll_threshold {
        return finish(start, t, dt);
    }
    // The slip cannot reach the threshold within the step: one RK4 step
    // with the slip direction tracked through the stages.
    let fastest = c + patch.g_t.length();
    if s0 - roll_threshold > fastest * dt {
        return (patch.slide_turning(start, t, dt, u_hat0), dt, false);
    }
    let full = patch.slide_turning(start, t, dt, u_hat0);
    let s_full = patch.slip_speed(&full);
    if s_full > roll_threshold {
        return (full, dt, false);
    }
    // The slip falls to the threshold within the step: solve for that time,
    // then finish along the slip direction there.
    let above =
        |tau: f64| patch.slip_speed(&patch.slide_turning(start, t, tau, u_hat0)) - roll_threshold;
    let rate0 = c - patch.g_t.dot(u_hat0);
    let guess = if rate0 > 0.0 {
        (s0 - roll_threshold) / rate0
    } else {
        dt
    };
    let tau = find_root(
        above,
        (0.0, s0 - roll_threshold),
        (dt, s_full - roll_threshold),
        guess,
    );
    let mid = patch.slide_turning(start, t, tau, u_hat0);
    let (end, h, rolls) = finish(mid, t + tau, dt - tau);
    (end, tau + h, rolls)
}

#[cfg(test)]
mod regression_data;
#[cfg(test)]
mod terrain_tests;
#[cfg(test)]
mod tests;
