use crate::Vec3;

use super::{BallState, GroundSample, Phase, find_root, holds_at_rest};
use crate::aero::AeroModel;
use crate::rk4::{self, Integrable};

/// Contact-point slip speed, or rolling speed, (m/s) treated as zero.
/// Slides end exactly at zero slip and rolls exactly at zero speed; this
/// only absorbs rounding.
const SLIP_EPS: f64 = 1e-9;

/// Smallest slip speed (m/s) at which a slide on a slope switches to
/// finishing along a fixed slip direction, whatever
/// [`GroundPhysics::roll_threshold`] is.
const MIN_SWITCH_SLIP: f64 = 1e-6;

/// Position, velocity and spin of a ball on the ground, for RK4.
#[derive(Debug, Clone, Copy)]
struct Motion {
    pos: Vec3,
    vel: Vec3,
    spin: Vec3,
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
pub(super) fn slip(vel: Vec3, spin: Vec3, n: Vec3, radius: f64) -> Vec3 {
    (vel - vel.dot(n) * n) + spin.cross(-radius * n)
}

/// Remove the contact slip with the friction impulse that brings the ball
/// to rolling (angular momentum about the contact point is conserved).
pub(super) fn lock_rolling(vel: Vec3, spin: Vec3, n: Vec3, radius: f64, k: f64) -> (Vec3, Vec3) {
    let u = slip(vel, spin, n, radius);
    (
        vel - u * (k / (1.0 + k)),
        spin + n.cross(u) / ((1.0 + k) * radius),
    )
}

/// Decay the spin component about `n` from `start` to `end` without letting
/// it cross zero (pivot friction saturates).
fn saturate_pivot(start: Vec3, end: Vec3, n: Vec3) -> Vec3 {
    let w0 = start.dot(n);
    let w1 = end.dot(n);
    if w0 == 0.0 || w0.signum() != w1.signum() {
        end - w1 * n
    } else {
        end
    }
}

/// Tangential gravity `g − (g·n)n` on ground with unit normal `n`, for
/// gravity of magnitude `g` along −Y: the downhill pull (m/s²), zero on
/// level ground.
pub(super) fn tangential_gravity(n: Vec3, g: f64) -> Vec3 {
    let gravity = Vec3::new(0.0, -g, 0.0);
    gravity - gravity.dot(n) * n
}

/// Put a ball that has just stepped across the ground back onto it.
///
/// `at` is the terrain under the ball's new (x, z) and `prev_normal` the
/// normal the step was taken on. The center goes to the terrain height and
/// the velocity into the tangent plane at `at`; a ball that ended the step
/// `rolling` (a roll, or a slide that reached zero slip) also gets its
/// tangential spin re-locked to that velocity, keeping its spin about the
/// normal. Where the normal and height already match (level ground, a
/// plane) the state is returned unchanged.
pub(super) fn onto_surface(
    state: &BallState,
    at: &GroundSample,
    prev_normal: Vec3,
    rolling: bool,
    radius: f64,
) -> BallState {
    #[allow(
        clippy::float_cmp,
        reason = "only an exact match may skip the projection"
    )]
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
    n: Vec3,
    /// Tangential gravity (m/s²).
    g_t: Vec3,
    /// Sliding friction magnitude `μ_slide·g·n_y` (m/s²).
    friction: f64,
    /// Ball radius (m).
    radius: f64,
    /// Inertia factor k.
    k: f64,
    /// Pivot-friction angular acceleration about `n` for this step.
    pivot_accel: Vec3,
    /// Plowing rate (1/s): a resistance of `plow` times the velocity.
    plow: f64,
}

impl Patch {
    /// Rate of change of a sliding ball with slip direction `u_hat`.
    fn slide_rate(&self, s: &Motion, u_hat: Vec3) -> Motion {
        Motion {
            pos: s.vel,
            vel: -self.friction * u_hat + self.g_t - self.plow * s.vel,
            spin: self.friction * self.n.cross(u_hat) / (self.k * self.radius) + self.pivot_accel,
        }
    }

    /// One RK4 slide step of length `h` along the fixed slip direction
    /// `u_hat`.
    fn slide_fixed(&self, s: Motion, t: f64, h: f64, u_hat: Vec3) -> Motion {
        let accel = -self.friction * u_hat + self.g_t;
        let spin_accel =
            self.friction * self.n.cross(u_hat) / (self.k * self.radius) + self.pivot_accel;
        rk4::step_state(s, t, h, |s, _| Motion {
            pos: s.vel,
            vel: accel - self.plow * s.vel,
            spin: spin_accel,
        })
    }

    /// One RK4 slide step of length `h`, with the slip direction evaluated
    /// at every stage (`fallback` where the slip vanishes).
    fn slide_turning(&self, s: Motion, t: f64, h: f64, fallback: Vec3) -> Motion {
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

/// Angular acceleration of pivot friction, `pivot_rate` (rad/s²) against
/// the spin `w_n` about the normal `n`. Pivot friction opposes that spin;
/// with none (of either sign of zero) there is none.
fn pivot_friction(w_n: f64, pivot_rate: f64, n: Vec3) -> Vec3 {
    if w_n == 0.0 {
        Vec3::ZERO
    } else {
        -w_n.signum() * pivot_rate * n
    }
}

/// One ground step of at most `dt` on the ground sample `at`; see
/// [`GroundPhysics::step`]. Returns the new state, the phase over the step
/// and whether the ball ends it rolling. The ball is left in the tangent
/// plane at `at`; putting it back on the terrain is up to the caller.
pub(super) fn contact_step(
    state: &BallState,
    at: &GroundSample,
    model: &dyn AeroModel,
    dt: f64,
    stop_speed: f64,
    roll_threshold: f64,
) -> (BallState, Phase, bool) {
    let n = at.normal.normalize();
    let surface = &at.surface;
    let radius = model.radius_m();
    let k = model.inertia_factor();
    let inertia = k * radius * radius; // per unit mass

    // Normal load and tangential gravity per unit mass.
    let g_n = model.gravity() * n.y;
    let g_t = tangential_gravity(n, model.gravity());

    let vel = state.vel - state.vel.dot(n) * n;
    let spin = state.spin;
    let w_n = spin.dot(n);
    let pivot_accel = pivot_friction(w_n, surface.pivot * g_n / inertia, n);

    if slip(vel, spin, n, radius).length() > SLIP_EPS {
        let patch = Patch {
            n,
            g_t,
            friction: surface.mu_slide * g_n,
            radius,
            k,
            pivot_accel,
            plow: surface.plow,
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
    let decel = |v: f64| (surface.rolling_friction(v) * g_n + surface.plow * v) / (1.0 + k);
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
        holds_at_rest(n, surface, model.gravity())
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
            vel: Vec3::ZERO,
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
    // d|u|/dt = −c + g_t·û − plow·v·û: the slip shrinks at `c` less the
    // slope's pull along it, and plowing slows the ball along it. The switch
    // to finishing along a fixed slip direction is kept above zero slip,
    // where û would flip within the RK4 stages.
    let roll_threshold = roll_threshold.max(MIN_SWITCH_SLIP);
    let c = (1.0 + patch.k) / patch.k * patch.friction;
    let finish = |s: Motion, t: f64, left: f64| -> (Motion, f64, bool) {
        if patch.plow > 0.0 {
            return finish_plowing(patch, s, t, left, c);
        }
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

    if patch.g_t == Vec3::ZERO && patch.plow == 0.0 {
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
    // with the slip direction tracked through the stages. Plowing changes
    // the slip by at most `plow` times the fastest speed within the step.
    let g_t_len = patch.g_t.length();
    let top_speed = start.vel.length() + (patch.friction + g_t_len) * dt;
    let fastest = c + g_t_len + patch.plow * top_speed;
    if s0 - roll_threshold > fastest * dt {
        return (patch.slide_turning(start, t, dt, u_hat0), dt, false);
    }
    // With plowing, within `s0 / fastest` no RK4 stage sees the slip pass
    // zero and flip.
    let dt = if patch.plow > 0.0 {
        dt.min(s0 / fastest)
    } else {
        dt
    };
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

/// Slide along the fixed slip direction of `s` for at most `left` with
/// plowing: the state reached, the time taken, and whether the slip reached
/// zero. Plowing makes the slip fall at `c − g_t·û + plow·v·û`, which
/// changes with the velocity, so the zero-slip time is solved for.
fn finish_plowing(patch: &Patch, s: Motion, t: f64, left: f64, c: f64) -> (Motion, f64, bool) {
    let u = slip(s.vel, s.spin, patch.n, patch.radius);
    let len = u.length();
    if len <= SLIP_EPS {
        return (s, 0.0, true);
    }
    let u_hat = u / len;
    let along = |m: &Motion| slip(m.vel, m.spin, patch.n, patch.radius).dot(u_hat);
    let end = patch.slide_fixed(s, t, left, u_hat);
    let at_end = along(&end);
    if at_end > SLIP_EPS {
        return (end, left, false);
    }
    if at_end >= 0.0 {
        return (end, left, true);
    }
    let rate = c - patch.g_t.dot(u_hat) + patch.plow * s.vel.dot(u_hat);
    let guess = if rate > 0.0 { len / rate } else { left };
    let tau = find_root(
        |tau| along(&patch.slide_fixed(s, t, tau, u_hat)),
        (0.0, len),
        (left, at_end),
        guess,
    );
    (patch.slide_fixed(s, t, tau, u_hat), tau, true)
}
