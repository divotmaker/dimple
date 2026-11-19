use crate::aero::AeroModel;
use crate::ground::{GroundPoint, Phase};
use crate::vec3::Vec3;

/// The hole's radius under the Rules of Golf (m): 4¼ in (108 mm) across.
pub const REGULATION_RADIUS_M: f64 = 0.053_975;

/// A hole on level ground: its center in dimple's frame and its radius.
///
/// Non-exhaustive so it can gain fields; build one with [`Cup::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct Cup {
    /// Center, downrange (m).
    pub x: f64,
    /// Center, right of the target line (m).
    pub z: f64,
    /// Radius (m).
    pub radius_m: f64,
}

/// What a ball reaching the hole's edge does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Outcome {
    /// It drops in.
    Holed,
    /// It runs onto the rim, rolls around it and leaves.
    LipOut,
    /// It flies across the hole and out, or its path misses the hole.
    Over,
}

/// A ball dropping into the hole.
///
/// Non-exhaustive so it can gain fields.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct Capture {
    /// Absolute shot time (s) the ball's center crosses the hole's edge.
    pub time: f64,
    /// Ball center there (m).
    pub pos: Vec3,
    /// Horizontal speed there (m/s).
    pub speed_mps: f64,
    /// Distance of the ball's path from the hole's center (m).
    pub offset_m: f64,
}

impl Cup {
    /// A regulation hole ([`REGULATION_RADIUS_M`]) centered at (x, z).
    #[must_use]
    pub const fn new(x: f64, z: f64) -> Self {
        Self {
            x,
            z,
            radius_m: REGULATION_RADIUS_M,
        }
    }

    /// This hole with another radius (m).
    #[must_use]
    pub const fn with_radius(self, radius_m: f64) -> Self {
        Self { radius_m, ..self }
    }

    /// What a rolling ball does when its center reaches the hole's edge at
    /// horizontal speed `speed_mps`, on a straight path passing `offset_m`
    /// from the center (either side).
    ///
    /// The ball's radius, its inertia factor and gravity come from `model`.
    /// With `r` the ball's radius, `R` the hole's, `δ` the offset, `k` the
    /// inertia factor and `g` gravity:
    ///
    /// - Above `c1 = √(g·r / (1 − δ²/R²))` the ball leaves the edge at once
    ///   and flies across the hole. It meets the far edge after travelling
    ///   `x` along the chord `2√(R² − δ²)`, at a contact angle `θ` above the
    ///   ball's equator. A strike with no bounce and no friction keeps only
    ///   the velocity along the rim's tangent there, so the ball drops in if
    ///   its flight path is at least as steep as `θ`, i.e. it is heading at or
    ///   below the edge, and flies [`Over`](Outcome::Over) otherwise. A ball
    ///   that reaches the hole's wall with its center below the edge drops in.
    ///   Through the center the limit is 1.61 m/s.
    /// - At or below `c1` it rolls onto the rim. It drops in below
    ///   `c2 = √( (2g·r/(1+k)) / (R·δ²/(R−r)³ − 1) )` (it falls off the
    ///   rim into the hole) or below the separatrix `c3` (it pitches into the
    ///   hole), and lips out between those and `c1`. `c3` is given in the
    ///   rim's saddle angle `β ∈ (0, π/2)` by `v² = g·D / ((1+k)·sin β)` and
    ///   `δ² = cos β·(R − r cos β)³ / D`, with
    ///   `D = cos β·(R − r cos β) − 2r·sin β·(1 − sin β)`; it is unbounded
    ///   for `δ ≤ R − r`, and `c2` for `R·δ² ≤ (R − r)³`.
    ///
    /// The green is level and the ball has no spin about the vertical. An
    /// offset of `R` or more misses the hole.
    #[must_use]
    pub fn outcome(&self, speed_mps: f64, offset_m: f64, model: &dyn AeroModel) -> Outcome {
        let rim = Rim {
            cup: self.radius_m,
            ball: model.radius_m(),
            inertia: model.inertia_factor(),
            gravity: model.gravity(),
        };
        rim.outcome(speed_mps, offset_m.abs())
    }

    /// The first time a ball sliding or rolling along `points` drops into
    /// this hole, or `None`.
    ///
    /// `points` are dimple's ground points (e.g.
    /// [`GroundResult::points`](crate::ground::GroundResult::points)). Each
    /// segment the ball slides or rolls along (never a hop: an airborne ball
    /// is not captured) that carries its center from outside the hole's
    /// edge to inside is judged by [`outcome`](Self::outcome) at the
    /// crossing, with the speed interpolated there and the offset of the
    /// segment's line from the center. A ball that lips out or flies over
    /// runs on along `points`; a path that starts inside the edge, or lands
    /// inside it from a hop, has no crossing.
    #[must_use]
    pub fn capture(&self, points: &[GroundPoint], model: &dyn AeroModel) -> Option<Capture> {
        let radius = self.radius_m;
        points.windows(2).find_map(|pair| {
            let (a, b) = (&pair[0], &pair[1]);
            if a.phase == Phase::Hop {
                return None;
            }
            let (ax, az) = (a.pos.x - self.x, a.pos.z - self.z);
            let (sx, sz) = (b.pos.x - a.pos.x, b.pos.z - a.pos.z);
            let len2 = sx * sx + sz * sz;
            let c = ax * ax + az * az - radius * radius;
            if len2 == 0.0 || c <= 0.0 {
                return None;
            }
            // |a + t·s|² = R²: the smaller root is the crossing inward.
            let half_b = ax * sx + az * sz;
            let disc = half_b * half_b - len2 * c;
            if disc < 0.0 {
                return None;
            }
            let t = (-half_b - disc.sqrt()) / len2;
            if !(0.0..=1.0).contains(&t) {
                return None;
            }
            let speed = |p: &GroundPoint| p.vel.x.hypot(p.vel.z);
            let speed_mps = speed(a) + (speed(b) - speed(a)) * t;
            let offset_m = (ax * sz - az * sx).abs() / len2.sqrt();
            (self.outcome(speed_mps, offset_m, model) == Outcome::Holed).then(|| Capture {
                time: a.time + (b.time - a.time) * t,
                pos: a.pos.lerp(b.pos, t),
                speed_mps,
                offset_m,
            })
        })
    }
}

/// The ball and the hole: radii (m), inertia factor and gravity (m/s²).
#[derive(Debug, Clone, Copy)]
struct Rim {
    cup: f64,
    ball: f64,
    inertia: f64,
    gravity: f64,
}

impl Rim {
    /// [`Cup::outcome`] for an offset `d` ≥ 0.
    fn outcome(&self, v: f64, d: f64) -> Outcome {
        if d >= self.cup {
            return Outcome::Over;
        }
        if v > self.c1(d) {
            if self.far_rim_holds(v, d) {
                Outcome::Holed
            } else {
                Outcome::Over
            }
        } else if v < self.c2(d) || v < self.c3(d) {
            Outcome::Holed
        } else {
            Outcome::LipOut
        }
    }

    /// Speed above which the ball leaves the edge at once.
    fn c1(&self, d: f64) -> f64 {
        let q = d / self.cup;
        (self.gravity * self.ball / (1.0 - q * q)).sqrt()
    }

    /// Speed below which the ball falls off the rim into the hole.
    fn c2(&self, d: f64) -> f64 {
        let (big, r) = (self.cup, self.ball);
        let k = big * d * d / (big - r).powi(3) - 1.0;
        if k <= 0.0 {
            return f64::INFINITY;
        }
        (2.0 * self.gravity * r / (1.0 + self.inertia) / k).sqrt()
    }

    /// The rim's separatrix: below it a ball on the rim pitches into the
    /// hole.
    fn c3(&self, d: f64) -> f64 {
        let (big, r) = (self.cup, self.ball);
        if d <= big - r {
            return f64::INFINITY;
        }
        // δ rises from R − r to R as β goes from 0 to π/2.
        let saddle = |beta: f64| {
            let (s, c) = beta.sin_cos();
            let rc = big - r * c;
            let den = c * rc - 2.0 * r * s * (1.0 - s);
            (den, (c * rc.powi(3) / den).sqrt(), s)
        };
        let (mut lo, mut hi) = (0.0, std::f64::consts::FRAC_PI_2);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if saddle(mid).1 < d {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let (den, _, s) = saddle(0.5 * (lo + hi));
        (self.gravity * den / ((1.0 + self.inertia) * s)).sqrt()
    }

    /// Whether a ball leaving the near edge at `v`, `d` off center, drops in
    /// when it meets the far side.
    fn far_rim_holds(&self, v: f64, d: f64) -> bool {
        let (big, r, g) = (self.cup, self.ball, self.gravity);
        let chord = 2.0 * (big * big - d * d).sqrt();
        // Center height above the edge after travelling x along the chord.
        let height = |x: f64| r - 0.5 * g * x * x / (v * v);
        let wall = (chord - r).max(0.0);
        if height(wall) <= 0.0 {
            return true; // reaches the wall below the edge
        }
        let level = v * (2.0 * r / g).sqrt();
        let end = chord.min(level);
        let gap = |x: f64| {
            let rise = height(x);
            (chord - x).powi(2) + rise * rise - r * r
        };
        // gap > 0 at `wall`, ≤ 0 at `end`, decreasing between.
        let (mut lo, mut hi) = (wall, end);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if gap(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let reach = hi;
        // Flight slope g·x/v² against the slope to the edge, h/(chord − x).
        g * reach / (v * v) * (chord - reach) >= height(reach)
    }
}

#[cfg(test)]
mod tests;
