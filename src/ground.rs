use std::fmt;

use crate::Vec3;

use crate::aero::AeroModel;
use crate::rk4;
use crate::units::G;

mod contact;
mod surface;

pub use surface::{Restitution, RollLaw, STIMPMETER_SPEED, Surface};

use contact::{contact_step, lock_rolling, onto_surface, slip, tangential_gravity};

/// Spacing (s) of the sampled [`GroundPoint`]s.
const SAMPLE_INTERVAL: f64 = 0.01;

/// Time tolerance (s) for reaching sample times and the time limit.
const TIME_EPS: f64 = 1e-9;

/// Normal speed, relative to the ball's speed, below which an impact counts
/// as grazing: rounding of a velocity that lies in the ground plane.
const GRAZING: f64 = 1e-12;

/// Residual accepted by the crossing solves: meters for ground crossings,
/// m/s for the slide-to-roll event.
const ROOT_TOL: f64 = 1e-12;

/// Iteration cap of the crossing solves; they converge in a handful.
const ROOT_MAX_ITER: usize = 200;

/// The ground at one horizontal location.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GroundSample {
    /// Height (m) of the ground surface at (x, z). dimple treats the ball as
    /// touching when its center reaches this height (the ball radius is not
    /// added), consistent with the flight's landing; on ground level with
    /// the tee this is 0 and a ball on the ground has `pos.y == 0`.
    pub height: f64,
    /// Unit outward surface normal at (x, z). From a [`Terrain`] it has
    /// `normal.y > 0`: for a heightfield `y = h(x, z)` it is
    /// `(−∂h/∂x, 1, −∂h/∂z)` normalized. [`GroundPhysics`] uses it for the
    /// impact plane, the normal load and the downhill pull, so it should
    /// match the slope of `height`. Given to [`impact`] on its own, it may
    /// point any way: a wall, a tree trunk or the underside of an overhang.
    pub normal: Vec3,
    /// Material at this location.
    pub surface: Surface,
}

impl GroundSample {
    /// Whether a ball at rest here stays put: the slope's pull, `g·sin θ`,
    /// does not exceed the rolling resistance at zero speed,
    /// `mu_roll·g·cos θ`, i.e. `tan θ ≤ mu_roll` (θ the slope angle). Always
    /// true on level ground. The answer does not depend on the strength of
    /// gravity.
    #[must_use]
    pub fn holds_at_rest(&self) -> bool {
        holds_at_rest(self.normal.normalize(), &self.surface, G)
    }
}

/// Whether rolling resistance at zero speed holds a ball against the slope
/// of ground with unit normal `n`, under gravity of magnitude `g`.
fn holds_at_rest(n: Vec3, surface: &Surface, g: f64) -> bool {
    tangential_gravity(n, g).length() <= surface.mu_roll * g * n.y
}

/// Layer 1 of the ground API: "what is the ground like here?"
///
/// A heightfield point query: one ground height, normal and surface per
/// horizontal position. dimple never needs to know how the ground is
/// stored: an analytic plane, a grid, a mesh lookup or an engine query all
/// fit behind `sample`.
///
/// [`GroundPhysics`] tracks the ball's center and samples the terrain at
/// the center's (x, z): the ball is on the ground when `pos.y` equals the
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
///   leave the ground; that launch is not modeled.
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
            normal: Vec3::Y,
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
/// use dimple::Vec3;
///
/// let range = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| GroundSample {
///     height: 0.0,
///     normal: Vec3::Y,
///     surface: if (150.0..170.0).contains(&x) { Surface::GREEN } else { Surface::FAIRWAY },
/// }));
/// # let _ = range;
/// ```
///
/// A plane: fairway rising 3° downrange and falling 2° to the right. For
/// `height = a·x + b·z` the normal is `(−a, 1, −b)` normalized.
///
/// ```
/// use dimple::ground::{FnTerrain, GroundPhysics, GroundSample, Surface};
/// use dimple::trajectory::{ShotInput, simulate_shot};
/// use dimple::{Tour, Vec3};
///
/// let (a, b) = (3.0_f64.to_radians().tan(), -(2.0_f64.to_radians().tan()));
/// let normal = Vec3::new(-a, 1.0, -b).normalize();
/// let hill = GroundPhysics::new(FnTerrain(move |x: f64, z: f64| GroundSample {
///     height: a * x + b * z,
///     normal,
///     surface: Surface::FAIRWAY,
/// }));
/// let shot = simulate_shot(&ShotInput::seven_iron(), &Tour::SEA_LEVEL, &hill);
/// let rest = shot.ground.rest.pos;
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
    /// Position of the ball's center (m, crate frame).
    pub pos: Vec3,
    /// Velocity (m/s).
    pub vel: Vec3,
    /// Angular velocity (rad/s, right-hand rule, crate frame).
    pub spin: Vec3,
}

/// What the ball is doing between ground points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Phase {
    /// In the air between impacts: gravity plus the model's
    /// [`acceleration`](AeroModel::acceleration), spin held fixed.
    Hop,
    /// On the ground with the contact point slipping: sliding friction.
    Slide,
    /// On the ground, rolling without slipping: rolling resistance.
    Roll,
}

/// Why the ground phase ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
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
    /// Ball center at impact (m).
    pub pos: Vec3,
    /// Unit surface normal at the contact.
    pub normal: Vec3,
    /// Surface at the contact point.
    pub surface: Surface,
    /// Velocity arriving at the ground (m/s).
    pub vel_in: Vec3,
    /// Velocity leaving the contact (m/s). For the last contact before the
    /// ball stays on the ground this is the velocity it slides or rolls
    /// away with: the rebound projected onto the ground plane, with a slip
    /// at or below [`GroundPhysics::roll_threshold`] removed.
    pub vel_out: Vec3,
    /// Angular velocity arriving (rad/s).
    pub spin_in: Vec3,
    /// Angular velocity leaving (rad/s); for the last contact, the spin the
    /// ball slides or rolls away with.
    pub spin_out: Vec3,
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
    /// Ball center (m).
    pub pos: Vec3,
    /// Velocity (m/s).
    pub vel: Vec3,
    /// Angular velocity (rad/s).
    pub spin: Vec3,
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
/// Each impact uses [`impact`]; the ball keeps hopping while its rebound
/// would rise at least [`min_hop_height`](Self::min_hop_height) and, on a
/// slope, the hops die out rather than being sustained by the slope (see
/// `min_hop_height`), then slides and rolls ([`step`](Self::step)) until it
/// stops. Each hop is RK4 flight under gravity and the acceleration of the
/// model passed to [`settle`](GroundModel::settle), with the rebound's spin
/// held fixed. Gravity is the model's [`gravity`](AeroModel::gravity) here
/// and in every contact. The tunables are
/// public fields; start from [`GroundPhysics::new`] and override with struct
/// update syntax.
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
    /// The ball also needs the hops to die out: the impact must remove at
    /// least the energy of a hop this high (`g·min_hop_height` per unit
    /// mass) more than the next hop hands back. Otherwise the ball goes to
    /// the ground at once, as below the height. Two things hand energy back:
    /// gravity down a slope, and wind pushing a ball slower than the air.
    /// Without this a steep enough downslope, or a strong enough wind, keeps
    /// the ball fast enough that its hops stop shrinking, and it never rolls
    /// (see `docs/AERO.md`). The energy handed back is measured on the hop
    /// actually flown, through whatever air there is. In still air on level
    /// ground every hop loses energy, so the height alone decides.
    pub min_hop_height: f64,
    /// Safety bound on impacts. Every hop rises at least
    /// [`min_hop_height`](Self::min_hop_height) and the hops die out, so
    /// real shots stop bouncing far below this; it guards degenerate surfaces
    /// (e.g. `Restitution::Constant(1.0)` with no friction). Default 200.
    pub max_impacts: usize,
    /// Time limit (s) on the ground phase, from landing. Default 120 s,
    /// enough for a ball to come to rest on a green sloped just inside its
    /// hold angle (a drive takes 82 s on a 3.1° `GREEN`). Where the slope
    /// is too steep to hold the ball ([`GroundSample::holds_at_rest`]) it
    /// never rests, and a long enough slope runs to this limit.
    pub max_time: f64,
}

impl<T: Terrain> GroundPhysics<T> {
    /// Default physics over `terrain`: 1 ms steps, 5 mm hop transition,
    /// 0.01 m/s stop and roll thresholds, at most 200 impacts and 120 s.
    #[must_use]
    pub fn new(terrain: T) -> Self {
        Self {
            terrain,
            dt: 0.001,
            stop_speed: DEFAULT_STOP_SPEED,
            roll_threshold: DEFAULT_ROLL_THRESHOLD,
            min_hop_height: 0.005,
            max_impacts: 200,
            max_time: 120.0,
        }
    }

    /// Advance a ball that is on the ground by one step of
    /// [`dt`](Self::dt), on the terrain under its center.
    ///
    /// Returns the new state and the phase the ball was in over the step.
    /// With `n` the sampled normal, normal load `g·n_y` and tangential
    /// gravity `g_t = g − (g·n)n` (zero on level ground), RK4 within each
    /// smooth phase:
    ///
    /// - **Slide** (the contact point slips): friction `−μ_slide·g·n_y·û`,
    ///   its torque `(−R n) × F`, `g_t`, and plowing `−plow·v`
    ///   ([`Surface::plow`], through the center, so no torque). Where `g_t`
    ///   and `plow` are zero the slip direction is fixed and the slip speed
    ///   falls linearly at `((1+k)/k)·μ_slide·g·n_y`, so a slide that would
    ///   end within the step stops exactly at the slide-to-roll time (the
    ///   step is shortened) and the ball leaves rolling. On a slope `g_t`
    ///   turns the slip direction, and plowing turns it wherever the
    ///   velocity is not along the slip, so û is re-evaluated at every RK4
    ///   stage; the time the slip falls to
    ///   [`roll_threshold`](Self::roll_threshold) is solved within the step,
    ///   and from there the slide finishes along that slip direction,
    ///   ending at zero slip (with plowing, at the solved zero-slip time).
    /// - **Roll**: spin tangential to the ground is locked to the velocity;
    ///   deceleration `(μ_eff·g·n_y + plow·v)/(1+k)` against the velocity,
    ///   with `μ_eff` from [`Surface::rolling_friction`], plus `g_t/(1+k)`. A ball
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
    pub fn step(&self, state: &BallState, model: &dyn AeroModel) -> (BallState, Phase) {
        let at = self.sample(state.pos);
        let (next, phase, rolling) = contact_step(
            state,
            &at,
            model,
            self.dt,
            self.stop_speed,
            self.roll_threshold,
        );
        let under = self.sample(next.pos);
        (
            onto_surface(&next, &under, at.normal, rolling, model.radius_m()),
            phase,
        )
    }

    /// Terrain under a ball centered at `pos`.
    fn sample(&self, pos: Vec3) -> GroundSample {
        self.terrain.sample(pos.x, pos.z)
    }

    /// Height of a ball centered at `pos` above the terrain under it.
    fn clearance(&self, pos: Vec3) -> f64 {
        pos.y - self.sample(pos).height
    }

    /// Flight from a rebound until the ball meets the ground again:
    /// `Ok(arrival)`, or `Err(state)` if the time limit `t_end` comes first.
    /// Sampled points are appended to `points`.
    ///
    /// RK4 in steps of [`dt`](Self::dt) under the model's gravity and
    /// [`acceleration`](AeroModel::acceleration) (handed the ball's height
    /// above the terrain), with the rebound's spin held fixed for the whole
    /// hop: the arrival carries it unchanged.
    fn hop(
        &self,
        start: &BallState,
        t_end: f64,
        points: &mut Vec<GroundPoint>,
        model: &dyn AeroModel,
    ) -> Result<BallState, BallState> {
        let gravity = model.gravity();
        let down = Vec3::new(0.0, -gravity, 0.0);
        let accel = |pos: Vec3, vel: Vec3, time: f64| {
            let ball = BallState {
                time,
                pos,
                vel,
                spin: start.spin,
            };
            down + model.acceleration(&ball, self.clearance(pos))
        };
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
            let b = rk4::step(s, t, h, accel);
            let end = self.clearance(b.pos);
            let arrival =
                self.hop_arrival(t, h, end, start.spin, gravity, |tau| hermite(s, b, h, tau));
            if let Some(arrival) = arrival {
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

    /// First ground crossing within a hop step of length `h` starting at
    /// time `t`, if there is one; `end` is the clearance at the end of the
    /// integrated step, `g` the magnitude of gravity and `arc(τ)` the ball's
    /// state `τ` into the step (`arc(0)` is the step's start).
    ///
    /// [`hop`](Self::hop) passes the cubic Hermite interpolant of the step's
    /// end states (position and velocity at both ends), so `arc(h)` is the
    /// integrated step end; under gravity alone it is the exact parabola,
    /// which RK4 reproduces. The clearance `y(τ) − height(x(τ), z(τ))` is
    /// checked at [`CROSSING_SUBSTEPS`] equal sub-steps for the first at or
    /// below the ground, and the crossing is the root inside that sub-step,
    /// solved by [`find_root`] from the ballistic drop time onto the height
    /// at the end of the step (exact on level ground without air). The
    /// arrival is placed exactly on the ground, with the arc's velocity
    /// there and `spin`.
    fn hop_arrival(
        &self,
        t: f64,
        h: f64,
        end: f64,
        spin: Vec3,
        g: f64,
        arc: impl Fn(f64) -> rk4::State,
    ) -> Option<BallState> {
        let clearance = |tau: f64| self.clearance(arc(tau).pos);
        let (lo, (hi, f_hi)) = first_crossing(clearance, h, clearance(0.0), end)?;
        #[allow(
            clippy::float_cmp,
            reason = "`hi` is `h` itself when the crossing is in the last sub-step"
        )]
        let f_hi = if hi == h { clearance(h).min(0.0) } else { f_hi };
        let a = arc(0.0);
        let end_height = self.sample(arc(h).pos).height;
        let guess = drop_time(a.pos.y - end_height, a.vel.y, g).clamp(lo.0, hi);
        let tau = find_root(clearance, lo, (hi, f_hi), guess);
        let arrival = arc(tau);
        let mut pos = arrival.pos;
        pos.y = self.sample(pos).height;
        Some(BallState {
            time: t + tau,
            pos,
            vel: arrival.vel,
            spin,
        })
    }
}

/// The state `tau` into an integration step of length `h` from `a` to `b`:
/// the cubic Hermite interpolant of the position through both end states
/// (position and velocity), and its derivative for the velocity. Exact at
/// `tau` = 0 and `tau` = `h`; within a 1 ms hop step it is within 1e-14 m
/// of the integrated path.
fn hermite(a: rk4::State, b: rk4::State, h: f64, tau: f64) -> rk4::State {
    let u = tau / h;
    let (u2, u3) = (u * u, u * u * u);
    let chord = b.pos - a.pos;
    let pos = a.pos
        + a.vel * ((u3 - 2.0 * u2 + u) * h)
        + chord * (3.0 * u2 - 2.0 * u3)
        + b.vel * ((u3 - u2) * h);
    let vel = a.vel * (3.0 * u2 - 4.0 * u + 1.0)
        + chord * ((6.0 * u - 6.0 * u2) / h)
        + b.vel * (3.0 * u2 - 2.0 * u);
    rk4::State { pos, vel }
}

/// Time (s) for a ball `height` above the ground moving up at `vy` to fall
/// back to it under gravity `g` alone. The positive root of
/// `height + vy·τ − g·τ²/2 = 0`, written to avoid cancellation.
fn drop_time(height: f64, vy: f64, g: f64) -> f64 {
    if height <= 0.0 && vy <= 0.0 {
        return 0.0;
    }
    let disc = (vy * vy + 2.0 * g * height.max(0.0)).sqrt();
    if vy >= 0.0 {
        (vy + disc) / g
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
    /// when its center reaches it (see [`GroundSample::height`]). Used by
    /// the default [`landing`](Self::landing), to place the ball on the
    /// tee at (0, `height(0, 0)`, 0), to test whether a flight starts below
    /// the ground, and to give the [`AeroModel`] the ball's height above
    /// the ground at every RK4 stage of the flight. Keep it cheap: it is
    /// called several times per millisecond of flight.
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
    ///
    /// `model` is the ball and the air it hops through between impacts
    /// ([`simulate_shot`](crate::trajectory::simulate_shot) passes the
    /// shot's): its [`acceleration`](AeroModel::acceleration) for the hops,
    /// its [`radius_m`](AeroModel::radius_m) and
    /// [`inertia_factor`](AeroModel::inertia_factor) for the contacts, and
    /// its [`gravity`](AeroModel::gravity) throughout. A ground model may
    /// use as much of it as it models, and may pass it on to [`impact`] and
    /// [`GroundPhysics::step`].
    #[must_use]
    fn settle(&self, landing: &BallState, model: &dyn AeroModel) -> GroundResult;
}

/// Flat ground with the default [`GroundPhysics`].
impl GroundModel for Surface {
    fn height(&self, _x: f64, _z: f64) -> f64 {
        0.0
    }

    fn settle(&self, landing: &BallState, model: &dyn AeroModel) -> GroundResult {
        GroundPhysics::new(*self).settle(landing, model)
    }
}

impl<T: Terrain> GroundModel for GroundPhysics<T> {
    fn height(&self, x: f64, z: f64) -> f64 {
        self.terrain.sample(x, z).height
    }

    fn settle(&self, landing: &BallState, model: &dyn AeroModel) -> GroundResult {
        let t_end = landing.time + self.max_time;
        let radius = model.radius_m();
        let g = model.gravity();
        let mut contacts = Vec::new();
        let mut points = Vec::new();
        let mut state = *landing;

        // Impacts and the hops between them.
        let mut at = loop {
            let at = self.sample(state.pos);
            let n = at.normal.normalize();
            let out = impact(&state, &at, model);
            let v_n = out.vel.dot(n);
            // A rebound that would clear the transition height flies its
            // hop. The flown hop must also show the train dying out
            // (`hops_die_out`); otherwise it is discarded and the ball stays
            // on the ground here.
            let mut flown = None;
            if v_n > 0.0 && v_n * v_n >= 2.0 * g * n.y * self.min_hop_height {
                if contacts.len() + 1 >= self.max_impacts {
                    contacts.push(Contact::new(&state, &at, &out));
                    points.push(GroundPoint::new(&out, Phase::Hop));
                    return GroundResult::new(contacts, points, out, Termination::MaxImpacts);
                }
                let mark = points.len();
                points.push(GroundPoint::new(&out, Phase::Hop));
                flown = match self.hop(&out, t_end, &mut points, model) {
                    Ok(arrival)
                        if !hops_die_out(&state, &out, &arrival, model, self.min_hop_height) =>
                    {
                        points.truncate(mark);
                        None
                    }
                    hop => Some(hop),
                };
            }
            let Some(hop) = flown else {
                // The ball stays on the ground: drop the (small or inward)
                // normal velocity and slide/roll from here.
                let mut leave = BallState {
                    vel: out.vel - v_n * n,
                    ..out
                };
                if slip(leave.vel, leave.spin, n, radius).length() <= self.roll_threshold {
                    (leave.vel, leave.spin) =
                        lock_rolling(leave.vel, leave.spin, n, radius, model.inertia_factor());
                }
                contacts.push(Contact::new(&state, &at, &leave));
                state = leave;
                break at;
            };
            contacts.push(Contact::new(&state, &at, &out));
            match hop {
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
                contact_step(&state, &at, model, h, self.stop_speed, self.roll_threshold);
            if state.time >= next_sample - TIME_EPS || last_phase != Some(phase) {
                points.push(GroundPoint::new(&state, phase));
                next_sample = state.time + SAMPLE_INTERVAL;
            }
            last_phase = Some(phase);
            let stopped = phase == Phase::Roll && next.vel == Vec3::ZERO;
            let under = self.sample(next.pos);
            let held = stopped && holds_at_rest(under.normal.normalize(), &under.surface, g);
            state = onto_surface(&next, &under, at.normal, rolling, radius);
            at = under;
            if held {
                points.push(GroundPoint::new(&state, Phase::Roll));
                return GroundResult::new(contacts, points, state, Termination::Rest);
            }
        }
    }
}

/// Whether the hops that follow a rebound die out, so that the ball may keep
/// bouncing; otherwise it goes to the ground at once, as below
/// [`GroundPhysics::min_hop_height`].
///
/// Each hop can hand energy back: gravity down a slope, and wind, whose drag
/// and lift push a ball moving slower than the air (on level ground in still
/// air a hop only loses energy). With `E = |v|²/2 + k·R²·|ω|²/2` (per unit
/// mass), the impact `arrival → rebound` removed `ΔE = E(arrival) −
/// E(rebound)` and the hop `rebound → next`, flown through whatever air
/// there is, returned `ΔU = E(next) − E(rebound)`: the potential energy
/// released down the slope plus the wind's work, less what drag took. The
/// train dies out if the hop returns nothing (`ΔU ≤ 0`) or the impact
/// removes at least the energy of a `min_hop_height` hop more than the hop
/// returns (`ΔE − ΔU ≥ g·min_hop_height`). Below that, the slope or the wind
/// re-supplies what each impact takes and the hops need not shrink through
/// the 5 mm transition: in a strong wind the ball would hop on until the
/// impact limit.
fn hops_die_out(
    arrival: &BallState,
    rebound: &BallState,
    next: &BallState,
    model: &dyn AeroModel,
    min_hop_height: f64,
) -> bool {
    // The hop keeps its spin, so its net return is all translational.
    let returned = 0.5 * (next.vel.length_squared() - rebound.vel.length_squared());
    // A ballistic hop on level ground returns exactly nothing; let its
    // rounding error count as nothing too.
    if returned <= 1e-9 * rebound.vel.length_squared() {
        return true;
    }
    let radius = model.radius_m();
    let spin_inertia = model.inertia_factor() * radius * radius;
    let energy =
        |s: &BallState| 0.5 * s.vel.length_squared() + 0.5 * spin_inertia * s.spin.length_squared();
    let removed = energy(arrival) - energy(rebound);
    removed - returned >= model.gravity() * min_hop_height
}

/// One impact of the ball with the ground or any other surface: the
/// post-impact state (same time and position, new velocity and spin).
///
/// `at.normal` is the surface's outward normal and may point any way: the
/// ground (`normal.y > 0`), a wall or tree trunk (horizontal), or the
/// underside of an overhang (`normal.y < 0`); `at.height` is not used. The
/// spin acts on every such bounce through the friction impulse below:
/// backspin checks a ball off a cart path and sidespin kicks it sideways
/// off a wall. The `normal.y > 0` restriction belongs to [`Terrain`] and
/// [`GroundPhysics`], which need a heightfield. For a hard surface (a cart
/// path, a wall, bark) use `compliance: 0.0`, a rigid flat contact, and a
/// [`Restitution::Constant`]; dimple has no published values for them.
///
/// A rigid-sphere point-contact impulse against Penner's (2002) effective
/// contact plane, in vector form. With `n` the surface normal, `t̂` the
/// direction of the incoming tangential velocity, `R` the radius and `k` the
/// inertia factor:
///
/// 1. Tilt: `tan β₀ = (v_t / 50 m/s) · (1 − exp(−(θ / 20°)³))` and
///    β = min(`compliance` · β₀, incidence angle), with `v_n = −v·n` the
///    normal impact speed, `v_t` the tangential speed and
///    θ = atan(`v_n` / `v_t`) the landing angle from the surface; effective
///    normal `n' = n cos β − t̂ sin β`, leaning back toward the incoming
///    ball. The law is fitted to 693 measured bounces on natural tee turf
///    at 2–39 m/s and 16–90° (Biber et al. 2023). The tilt follows the
///    tangential speed on steep landings, falls off for landings flatter
///    than about 20° and vanishes as `v_n` goes to zero, so the low,
///    near-grazing hops after a landing die out in a few impacts and the
///    ball goes on to slide and roll. No impact is measured below 16°.
/// 2. Normal impulse (per unit mass) `j_n = (1 + e)·|v·n'|`, with `e` from
///    [`Surface::restitution`] at the normal speed `|v·n'|`.
/// 3. Contact-point slip `u = (v − (v·n')n') + ω × (−R n')`. Friction
///    `j_t = min(|u|·k/(1+k), μ·j_n)` acts along `−û`: the ball either grips
///    (leaves with zero slip) or slides throughout (Coulomb limit).
/// 4. `v += j_n n' − j_t û`, `ω += (−R n') × (−j_t û) / (k R²)`.
/// 5. Pivot friction on the spin about the true normal `n` (not `n'`, so it
///    never removes spin along the heading):
///    `Δω_n = −sign(ω_n) · min(|ω_n|, pivot · j_n / (k R²))`.
/// 6. Absorption: with θ the landing angle from the surface
///    (`sin θ = v_n / |v|`), the ball leaves with
///    `kept = (1 − absorption)^(sin θ)` of the velocity and spin from steps
///    1–5. [`Surface::absorption`] is `0.0` on every preset, so `kept` is
///    exactly 1 there; a value above 1 counts as 1 (the ball stops).
///
/// No impulse is applied if the ball is not moving into the surface
/// (`v·n ≥ 0`, grazing; a normal speed within rounding of zero,
/// `|v·n| ≤ 1e-12·|v|`, counts as grazing, so a ball moving in a sloped
/// ground plane is not struck by its own rounding error). Whether the ball then hops or stays on the ground
/// is up to the caller ([`GroundPhysics`] uses the 5 mm hop rule; see
/// [`GroundPhysics::min_hop_height`]).
#[must_use]
pub fn impact(state: &BallState, at: &GroundSample, model: &dyn AeroModel) -> BallState {
    let n = at.normal.normalize();
    let v = state.vel;
    let v_dot_n = v.dot(n);
    if v_dot_n >= -GRAZING * v.length() {
        return *state;
    }
    let surface = &at.surface;
    let radius = model.radius_m();
    let k = model.inertia_factor();

    // Penner's effective contact plane.
    let v_t_vec = v - v_dot_n * n;
    let v_t = v_t_vec.length();
    let beta = surface.contact_angle(v_t, -v_dot_n);
    let t_hat = if v_t > 0.0 { v_t_vec / v_t } else { Vec3::ZERO };
    let n_eff = n * beta.cos() - t_hat * beta.sin();
    let v_n_eff = v.dot(n_eff);
    if v_n_eff >= 0.0 {
        return *state;
    }

    let e = surface.restitution.coefficient_at(v_n_eff);
    let j_n = (1.0 + e) * -v_n_eff;

    let lever = -radius * n_eff;
    let slip = (v - v_n_eff * n_eff) + state.spin.cross(lever);
    let slip_speed = slip.length();
    let j_t = (slip_speed * k / (1.0 + k)).min(surface.mu_slide * j_n);
    let j_t_vec = if slip_speed > 0.0 {
        slip * (-j_t / slip_speed)
    } else {
        Vec3::ZERO
    };

    let inertia = k * radius * radius; // per unit mass
    let vel = v + j_n * n_eff + j_t_vec;
    let mut spin = state.spin + lever.cross(j_t_vec) / inertia;

    let w_n = spin.dot(n);
    let dw = w_n.abs().min(surface.pivot * j_n / inertia);
    spin -= n * (w_n.signum() * dw);

    // The ground keeps `1 − absorption` of a straight-down rebound, and more
    // of a shallower one.
    let kept = (1.0 - surface.absorption)
        .max(0.0)
        .powf(-v_dot_n / v.length());
    BallState {
        vel: vel * kept,
        spin: spin * kept,
        ..*state
    }
}

/// Test fixture for the contact, slide and roll mechanics: every field
/// spelled out, so a change to a preset never moves the tests built on it.
/// Sliding friction 0.40, a constant rolling resistance of 0.055 (stimp
/// 14.3 ft) and no plowing, so the closed forms of a slide and a roll hold
/// exactly.
#[cfg(test)]
pub(crate) const TEST_TURF: Surface = Surface {
    restitution: Restitution::Penner2002,
    mu_slide: 0.40,
    mu_roll: 0.055,
    roll_law: RollLaw::Constant,
    compliance: 1.0,
    pivot: 0.001,
    absorption: 0.0,
    plow: 0.0,
};

#[cfg(test)]
#[allow(clippy::float_cmp, reason = "the tests assert exact results")]
mod absorption_tests;
#[cfg(test)]
mod regression_data;
#[cfg(test)]
#[allow(clippy::float_cmp, reason = "the tests assert exact results")]
mod terrain_tests;
#[cfg(test)]
#[allow(clippy::float_cmp, reason = "the tests assert exact results")]
mod tests;
#[cfg(test)]
mod validation_tests;
