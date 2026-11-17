use crate::units::{self, G};

/// Speed (m/s) at which a ball leaves the USGA stimpmeter onto the green,
/// already rolling: the 1.83 m/s (6.0 ft/s) of Penner (2002, "The physics of
/// putting", *Can. J. Phys.* 80, 83), from Holmes's (1986) analysis of the
/// device (*The Physics Teacher* 24, 401).
///
/// The ball is released from a notch 30 in up the ramp when the ramp reaches
/// about 20°. A solid ball rolling 30 in down a 20° plane without loss would
/// leave at about 1.91 m/s (6.3 ft/s); in the 145° V-groove the rails touch
/// the ball off its lowest point, so more of the drop goes into spin and the
/// lossless speed is about 1.88 m/s. Holmes's fuller analysis gives
/// 6.0 ft/s. dimple cites no direct measurement of it.
/// [`Surface::with_stimp`] and [`Surface::stimp_feet`] use this speed.
pub const STIMPMETER_SPEED: f64 = 1.83;

/// Inertia factor k = I/(m·R²) of the ball the stimpmeter conversions
/// assume: a uniform solid sphere, as in [`Tour`](crate::Tour).
const STIMP_INERTIA_FACTOR: f64 = 0.4;

/// Tangential speed (m/s) below which an impact has no tilt direction, so
/// the contact angle is zero.
const MIN_TILT_SPEED: f64 = 0.01;

/// Tangential impact speed scale V0 (m/s) of the contact angle: on a steep
/// landing `tan β = v_t / V0`. Fitted to 693 measured bounces on natural tee
/// turf (Biber et al. 2023; 95% interval 49.2–50.9 m/s).
const TILT_SPEED_SCALE: f64 = 50.0;

/// Landing angle scale θ0 (rad) of the contact angle, 20° from the surface:
/// the tilt falls off for landings flatter than about this. Fitted to the
/// same bounces (19.8–21.8°) and checked against the USGA's 2019 roll-out
/// line at 20–36°.
const TILT_LANDING_SCALE: f64 = 20.0 * (std::f64::consts::PI / 180.0);

/// Power p of the landing-angle factor `1 − exp(−(θ / θ0)^p)`. The bounces
/// allow about 2–4 and the 2019 roll-out line about 2–3.5; 3 is the lowest
/// round value in that band at which no grazing landing leaves faster along
/// the normal than it arrived.
const TILT_LANDING_POWER: f64 = 3.0;

/// Normal coefficient of restitution between ball and turf.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
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

/// Speed dependence of the rolling resistance:
/// `μ_eff(v)` = [`mu_roll`](Surface::mu_roll) · f(v), with f(0) = 1 for every
/// law, so `mu_roll` is always the zero-speed value.
///
/// Every preset uses [`Constant`](Self::Constant) and gets its rise with
/// speed from [`Surface::plow`], whose size is fitted to measured roll-outs.
/// The published speed-dependent fits of the resistance itself (Burritt &
/// Henrikson 2024) cover putts, 0–3 m/s. On slopes the published checks
/// favor `Constant`: Brede's (1991) slope correction of stimpmeter readings
/// is exact under it and Cross's (2024) speed-independent break angle holds
/// under it alone. See `docs/AERO.md`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum RollLaw {
    /// f = 1: constant rolling resistance, Penner's (2002) run model (his
    /// `ρ_g` = 0.131 reads about 6.0 ft on the stimpmeter).
    Constant,
    /// f = exp(v / `v_ref`): an exponential law of the kind Burritt &
    /// Henrikson (2024, "Modeling golf ball roll with a variable coefficient
    /// of friction", *Sports Eng.* 27:39) fit to camera-tracked putts,
    /// 0–3 m/s, one fit per surface. Positive `v_ref` makes a fast ball
    /// resist more (three of their four surfaces, `v_ref` roughly 1.6 to
    /// 4.9 m/s read from their figure); negative makes it resist less (one
    /// of their mats). `v_ref` (m/s) must be nonzero; ±∞ is
    /// [`Constant`](Self::Constant).
    ///
    /// The law grows quickly past the fitted range: `v_ref` = 1.6 m/s gives
    /// about 150 × `mu_roll` at 8 m/s, so a fast roll stops almost at once.
    /// With a negative `v_ref` a ball has no terminal speed on a slope
    /// steeper than its hold angle, and a fast roll on level ground can run
    /// to [`GroundPhysics::max_time`](super::GroundPhysics::max_time).
    Exponential {
        /// Speed scale (m/s); negative for resistance that falls with speed.
        v_ref: f64,
    },
}

impl RollLaw {
    /// The factor f(`speed`) above (`speed` in m/s, ≥ 0).
    ///
    /// A `v_ref` of zero or NaN, outside the documented range, is treated as
    /// [`Constant`](Self::Constant) (f = 1) rather than giving NaN.
    #[must_use]
    pub fn speed_factor(&self, speed: f64) -> f64 {
        match *self {
            Self::Constant => 1.0,
            Self::Exponential { v_ref } => {
                let s = speed / v_ref;
                if s.is_finite() { s.exp() } else { 1.0 }
            }
        }
    }

    /// I(v₀) = 2 ∫₀^v₀ v / f(v) dv (m²/s²): a ball rolling from `v0` on level
    /// ground travels (1+k)·I/(2·`mu_roll`·g). With s₀ = v₀/`v_ref`:
    ///
    /// - `Constant`: I = v₀²
    /// - `Exponential`: I = 2·`v_ref`²·(1 − e^(−s₀)·(1 + s₀))
    ///
    /// The exponential form tends to v₀² as `v_ref` → ±∞, which s₀ = 0
    /// gives exactly. NaN when `v_ref` is zero or NaN.
    fn stimp_integral(&self, v0: f64) -> f64 {
        let v0_sq = v0 * v0;
        match *self {
            Self::Constant => v0_sq,
            Self::Exponential { v_ref } => {
                if v_ref == 0.0 || v_ref.is_nan() {
                    return f64::NAN;
                }
                let s0 = v0 / v_ref;
                if s0 == 0.0 {
                    v0_sq
                } else if s0.abs() < 0.01 {
                    // (1 − e^(−s)(1 + s))/s² by its series, which the direct
                    // form loses to cancellation at small s.
                    let tail = 0.5
                        + s0 * (-1.0 / 3.0
                            + s0 * (1.0 / 8.0
                                + s0 * (-1.0 / 30.0 + s0 * (1.0 / 144.0 - s0 / 840.0))));
                    2.0 * v0_sq * tail
                } else {
                    2.0 * v_ref * v_ref * (1.0 - (-s0).exp() * (1.0 + s0))
                }
            }
        }
    }
}

/// Material properties of a patch of ground.
///
/// Plain public fields with `const` presets; adjust one with struct update
/// syntax, e.g. `Surface { mu_roll: 0.09, ..Surface::FAIRWAY }`. A `Surface`
/// is also a flat [`Terrain`](super::Terrain) (the plane y = 0 made of this surface) and a
/// complete [`GroundModel`](super::GroundModel) with the default [`GroundPhysics`](super::GroundPhysics).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Surface {
    /// Normal coefficient of restitution at impact.
    pub restitution: Restitution,
    /// Coulomb friction coefficient between the ball and the turf, used at
    /// impact and while the ball slides.
    pub mu_slide: f64,
    /// Rolling-resistance coefficient at zero speed: the effective friction
    /// of a slowly rolling ball, and what holds a ball at rest on a slope.
    pub mu_roll: f64,
    /// How the rolling resistance grows (or falls) with speed:
    /// `μ_eff` = `mu_roll` · f(v). Every preset uses `RollLaw::Constant`
    /// and gets its rise with speed from [`plow`](Self::plow).
    pub roll_law: RollLaw,
    /// Turf compliance: dimensionless multiplier on the effective contact
    /// angle β of [`impact`](super::impact),
    /// `tan β = (v_t / 50 m/s) · (1 − exp(−(θ / 20°)³))`, with `v_t` the
    /// tangential impact speed and θ the landing angle from the surface.
    /// `1.0` is firm, short-mown turf (the natural tee turf the law is fitted
    /// on, and the firm green of [`FIRM_GREEN`](Self::FIRM_GREEN)); higher
    /// is softer turf that the ball presses deeper into; `0.0` is a rigid
    /// flat surface. Measured single sites range 0.86–1.55.
    ///
    /// Compliance is how firmness enters the model. On greens it follows
    /// the penetration depth read by the USGA firmness meter: compliance ≈
    /// 0.090 × depth (mm) + 0.25 over six greens at 8.8–14.3 mm (R&A 2008;
    /// r = 0.96). On fairways a similar relation is indicative only (0.074 ×
    /// depth + 0.23 over five fairways at 8.9–15.4 mm, r = 0.70).
    pub compliance: f64,
    /// Torsional friction lever (m): the friction moment about the surface
    /// normal is `pivot` × normal force. It opposes spin about the normal,
    /// both at impact and on the ground, saturating at zero.
    ///
    /// It never changes the velocity of the impact it acts in, nor a slide
    /// or roll on level ground. The path depends on it only through spin
    /// about the normal carried into a later impact, where the tilted
    /// contact plane turns part of that spin into slip, and through the
    /// Magnus force that spin gives the hops in between. With the presets'
    /// 1 mm a full-shot landing removes all of it at once: a driver landing
    /// on `FAIRWAY` can remove up to about 160 rad/s (~1500 rpm), against
    /// the 47 rad/s a 600 rpm-sidespin drive lands with. Any lever that
    /// clears the spin at every impact gives the same path.
    ///
    /// The presets use 1 mm, the full-slip torsional couple of a Hertzian
    /// contact, (3π/16)·μ·a (Johnson 1985), with μ = 0.4 and a contact
    /// radius a ≈ 4 mm (a = 3–5 mm gives 0.7–1.2 mm). Uncalibrated.
    pub pivot: f64,
    /// Fraction of a rebound the ground absorbs, 0 to 1, stated for a ball
    /// that lands straight down. A ball that lands at angle θ from the
    /// surface leaves with `(1 − absorption)^(sin θ)` of the velocity and
    /// of the spin that the contact model above gives it: a steep landing
    /// loses nearly the whole fraction, a grazing touch almost nothing.
    /// `0.0` (every preset) is turf, which returns the whole rebound; `1.0`
    /// stops every landing on its mark with no spin.
    ///
    /// It stands for a bed the ball drives into, such as loose sand. The
    /// bed resists the ball's motion through it with a drag proportional to
    /// v² (Katsuragi & Durian 2007, 2013), which takes the same fraction of
    /// the velocity in every direction, and the ball drives in as deep as
    /// the normal part of its velocity takes it (Wright et al. 2020), which
    /// is the `sin θ`. That the spin loses the same fraction is a modeling
    /// choice; no published measurement covers it. It acts at impacts only;
    /// a ball rolling on the surface is slowed by `mu_roll` and `plow`. See
    /// `docs/AERO.md` for sand.
    pub absorption: f64,
    /// Plowing rate (1/s, ≥ 0): while the ball is on the ground, sliding or
    /// rolling, the turf resists it with a force per unit mass of `plow`
    /// times its velocity, on top of the sliding friction and the rolling
    /// resistance. A ball slowed by it alone travels `v / plow`. `0.0` is
    /// no plowing. A Stimpmeter reading fixes the resistance at 1.8 m/s, and
    /// measured fairway roll-outs that start at 2–6 m/s need about twice
    /// that, so the resistance has to rise with speed; `plow` gives it that
    /// rise. The form (resistance proportional to speed) is assumed and only
    /// its size is fitted: [`FAIRWAY`](Self::FAIRWAY)'s 0.36 to roll-outs on
    /// five fairways (R&A 2008) and the USGA's 2019 roll-out line (0.2–0.6
    /// fit about as well), [`GREEN`](Self::GREEN)'s 0.2 to the
    /// bounce-and-roll totals of 57 landings on one firm green (USGA and R&A
    /// 2007; 0.1–0.3 fit about as well). It does not act during an impact or
    /// a hop.
    ///
    /// Friction can hold a plowing ball rolling only below
    /// `(1+k)/k · mu_slide · g / plow` (about 17 m/s on `FAIRWAY` and 38 m/s
    /// on `GREEN`); a ball keeps rolling from its first zero slip even above
    /// that speed.
    pub plow: f64,
}

impl Surface {
    /// Standard fairway. `compliance` 1.05 and `mu_slide` 0.18 are the
    /// medians of eight fairways (R&A 2008: a 64 mph, 36° drive landing on
    /// five; USGA 2005: 252 landings at 23–48° on three at Pinehurst No. 2);
    /// single fairways range 0.86–1.55 and 0.12–0.33. Rolling resistance: a
    /// Stimpmeter reading of 5 ft 8 in on a creeping bentgrass fairway
    /// (USGA/R&A 2018 tests), with the constant law. Plowing 0.36 /s:
    /// fitted to the R&A roll-outs and the USGA 2019 roll-out line; 0.2–0.6
    /// fit about as well. See `docs/AERO.md`.
    pub const FAIRWAY: Self = Self {
        restitution: Restitution::Penner2002,
        mu_slide: 0.18,
        mu_roll: 0.0955,
        roll_law: RollLaw::Constant,
        compliance: 1.05,
        pivot: 0.001,
        absorption: 0.0,
        plow: 0.36,
    };

    /// Standard green. `compliance` 1.10 and `mu_slide` 0.22 are the medians
    /// of nine greens (R&A 2008: 5-iron and sand wedge landings on six; USGA
    /// 2005: 4-iron to pitching wedge on three); single greens range
    /// 1.05–1.52 and 0.17–0.26. Fitted to Cross's (2025) printed fits to 184
    /// bounces on four greens at 17–27 m/s and 35–55°, the two constants
    /// alone come out at 1.12 and 0.23.
    ///
    /// The friction is fitted to landings at 17–27 m/s. One green measured
    /// at 9–10 m/s (Cross 2025, Table 1) slides at 0.29; a chip with
    /// 3,300 rpm of backspin leaves that green 1.3 m/s slower than this
    /// preset gives.
    ///
    /// Rolling resistance: stimp 10.5 ft ([`stimp_feet`](Self::stimp_feet)),
    /// the reading of the one green with both a Stimpmeter reading and
    /// measured roll-outs (USGA and R&A 2007); use
    /// [`with_stimp`](Self::with_stimp) for another speed. Plowing 0.2 /s, fitted to that green's
    /// roll-outs.
    pub const GREEN: Self = Self {
        restitution: Restitution::Penner2002,
        mu_slide: 0.22,
        mu_roll: 0.0509,
        roll_law: RollLaw::Constant,
        compliance: 1.10,
        pivot: 0.001,
        absorption: 0.0,
        plow: 0.2,
    };

    /// Firm green: 57 impacts at 20–26 m/s and 36–59° on one green dried to
    /// championship firmness, Stimpmeter 10–11 ft (USGA and R&A 2007), fit
    /// at `compliance` 1.0 and `mu_slide` 0.17. Rolling resistance and
    /// plowing are [`GREEN`](Self::GREEN)'s, set on the same green.
    pub const FIRM_GREEN: Self = Self {
        mu_slide: 0.17,
        compliance: 1.0,
        ..Self::GREEN
    };

    /// This surface with `mu_roll` set so that it reads `feet` on the
    /// stimpmeter: a ball leaving the ramp rolling at [`STIMPMETER_SPEED`]
    /// rolls `feet` on level ground before it stops. Every other field,
    /// including [`roll_law`](Self::roll_law), is kept.
    ///
    /// Closed form: a rolling ball decelerates at `μ_eff·g/(1+k)` with
    /// `μ_eff = mu_roll·f(v)` from [`rolling_friction`](Self::rolling_friction),
    /// so from `v₀` it rolls `d = (1+k)·I / (2·mu_roll·g)` with
    /// `I = 2∫₀^v₀ v/f(v) dv`, and `mu_roll = (1+k)·I / (2·g·d)`. With
    /// s₀ = v₀/`v_ref`, per [`RollLaw`]:
    ///
    /// ```text
    /// Constant:     I = v₀²
    /// Exponential:  I = 2·v_ref²·(1 − e^(−s₀)·(1 + s₀))
    /// ```
    ///
    /// k = 2/5 (a uniform solid sphere, the inertia factor of
    /// [`Tour`](crate::Tour)). A ball with another inertia factor rolls
    /// `(1+k)/1.4` times this distance.
    ///
    /// With [`plow`](Self::plow) `c` > 0 the deceleration is
    /// `(mu_roll·f(v)·g + c·v)/(1+k)`, so `d = (1+k)·∫₀^v₀ v/(mu_roll·f(v)·g + c·v) dv`:
    /// for `Constant`, `d = (1+k)·[v₀/c − (μg/c²)·ln(1 + c·v₀/(μg))]`; for
    /// the other laws the integral is taken numerically. `d` falls as
    /// `mu_roll` grows, and `mu_roll` is found by bisection. Plowing alone
    /// stops the ball within `(1+k)·v₀/c`, so a reading at or beyond that
    /// (about 23 ft for [`FAIRWAY`](Self::FAIRWAY)'s 0.36 /s) has no
    /// `mu_roll`.
    ///
    /// Returns the surface unchanged if `feet` is not finite and positive,
    /// if the law's `v_ref` is zero or NaN, or if no `mu_roll` ≥ 0 reads
    /// `feet`.
    ///
    /// ```
    /// use dimple::ground::Surface;
    ///
    /// let fast = Surface::GREEN.with_stimp(13.0);
    /// assert!((fast.stimp_feet() - 13.0).abs() < 1e-9);
    /// assert!(fast.mu_roll < Surface::GREEN.mu_roll);
    /// ```
    #[must_use]
    pub fn with_stimp(self, feet: f64) -> Self {
        let distance = units::feet_to_meters(feet);
        if self.plow != 0.0 {
            return self.with_stimp_plowing(distance);
        }
        let integral = self.stimp_integral();
        if !(distance.is_finite() && distance > 0.0 && integral.is_finite()) {
            return self;
        }
        Self {
            mu_roll: (1.0 + STIMP_INERTIA_FACTOR) * integral / (2.0 * G * distance),
            ..self
        }
    }

    /// Stimpmeter reading (feet) of this surface: how far a ball leaving the
    /// ramp rolling at [`STIMPMETER_SPEED`] rolls on level ground before it
    /// stops, by the closed form in [`with_stimp`](Self::with_stimp)
    /// (k = 2/5), including [`plow`](Self::plow). 10.5 ft for
    /// [`GREEN`](Self::GREEN) and [`FIRM_GREEN`](Self::FIRM_GREEN), 5 ft
    /// 8 in for [`FAIRWAY`](Self::FAIRWAY).
    ///
    /// Infinite when `mu_roll` and `plow` are both zero (the ball never
    /// stops); NaN when the law's `v_ref` is zero or NaN.
    #[must_use]
    pub fn stimp_feet(&self) -> f64 {
        if self.plow != 0.0 {
            return units::meters_to_feet(self.stimp_distance_plowing(self.mu_roll));
        }
        let distance =
            (1.0 + STIMP_INERTIA_FACTOR) * self.stimp_integral() / (2.0 * G * self.mu_roll);
        units::meters_to_feet(distance)
    }

    /// The stimp distance times `2·mu_roll·g/(1+k)` (m²/s²), without
    /// plowing.
    fn stimp_integral(&self) -> f64 {
        self.roll_law.stimp_integral(STIMPMETER_SPEED)
    }

    /// Level-ground roll (m) from [`STIMPMETER_SPEED`] with rolling
    /// resistance `mu_roll` and this surface's law and plowing (k = 2/5):
    /// `(1+k)·∫₀^v₀ v/(mu_roll·f(v)·g + plow·v) dv`. NaN when the law's
    /// `v_ref` is zero or NaN.
    fn stimp_distance_plowing(&self, mu_roll: f64) -> f64 {
        let v0 = STIMPMETER_SPEED;
        if self.stimp_integral().is_nan() {
            return f64::NAN;
        }
        let friction = mu_roll * G;
        let plow = self.plow;
        let integral = if friction == 0.0 {
            v0 / plow
        } else if let RollLaw::Constant = self.roll_law {
            let x = plow * v0 / friction;
            friction / (plow * plow) * x_minus_ln_1p(x)
        } else {
            let law = self.roll_law;
            simpson(|v| v / (friction * law.speed_factor(v) + plow * v), v0)
        };
        (1.0 + STIMP_INERTIA_FACTOR) * integral
    }

    /// [`with_stimp`](Self::with_stimp) for a surface with plowing:
    /// bisection on `mu_roll` for a stimp roll of `distance` (m).
    fn with_stimp_plowing(self, distance: f64) -> Self {
        // Plowing alone (mu_roll 0) gives the longest roll.
        if !(distance.is_finite() && distance > 0.0 && distance < self.stimp_distance_plowing(0.0))
        {
            return self;
        }
        // The mu_roll that reads `distance` without plowing rolls shorter
        // with it, so it bounds the answer from above.
        let mut lo = 0.0;
        let mut hi = (1.0 + STIMP_INERTIA_FACTOR) * self.stimp_integral() / (2.0 * G * distance);
        if !hi.is_finite() {
            return self;
        }
        for _ in 0..STIMP_BISECTIONS {
            let mid = 0.5 * (lo + hi);
            if !(lo < mid && mid < hi) {
                break;
            }
            if self.stimp_distance_plowing(mid) > distance {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Self {
            mu_roll: hi,
            ..self
        }
    }

    /// Effective rolling-resistance coefficient at `speed` (m/s):
    /// `mu_roll · f(speed)` with f from [`roll_law`](Self::roll_law).
    /// [`plow`](Self::plow) is separate and not included.
    #[must_use]
    pub fn rolling_friction(&self, speed: f64) -> f64 {
        self.mu_roll * self.roll_law.speed_factor(speed)
    }

    /// Effective contact angle β (rad) for an impact with tangential speed
    /// `v_t` and normal approach speed `v_n` (both ≥ 0):
    ///
    /// ```text
    /// tan β₀ = (v_t / V0) · (1 − exp(−(θ / θ0)^p)),   θ = atan(v_n / v_t)
    /// β = min(compliance · β₀, incidence)
    /// ```
    ///
    /// with V0 = 50 m/s, θ0 = 20° and p = 3, fitted to 693 measured bounces
    /// on natural tee turf at 2–39 m/s and landing angles of 16–90° (Biber
    /// et al. 2023). On steep landings the factor in θ is 1 and the tilt
    /// follows the tangential speed alone; it falls off below about 20°
    /// and vanishes as `v_n` goes to zero, so the low hops after a landing
    /// meet nearly flat ground. The result never passes the incidence angle
    /// from the normal, where the ball would strike the effective plane
    /// head-on. Zero when there is (almost) no tangential motion, since the
    /// tilt direction is undefined.
    pub(super) fn contact_angle(&self, v_t: f64, v_n: f64) -> f64 {
        if v_t <= MIN_TILT_SPEED {
            return 0.0;
        }
        let incidence = v_t.atan2(v_n);
        let landing = v_n.atan2(v_t) / TILT_LANDING_SCALE;
        let factor = -(-landing.powf(TILT_LANDING_POWER)).exp_m1();
        let tilt = (v_t / TILT_SPEED_SCALE * factor).atan();
        (self.compliance * tilt).min(incidence)
    }
}

/// Bisection steps [`Surface::with_stimp`] allows when solving for
/// `mu_roll` with plowing; it stops earlier, once the bracket cannot shrink.
const STIMP_BISECTIONS: usize = 200;

/// Intervals of the composite Simpson rule for the stimp integral with
/// plowing (even).
const SIMPSON_INTERVALS: u32 = 256;

/// `x − ln(1 + x)` for x ≥ 0, by its series near zero, where the direct
/// form loses to cancellation.
fn x_minus_ln_1p(x: f64) -> f64 {
    if x < 0.01 {
        // x²/2 − x³/3 + … − x⁹/9
        let mut term = x;
        let mut sum = 0.0;
        for n in 2..10_i32 {
            term *= -x;
            sum -= term / f64::from(n);
        }
        sum
    } else {
        x - x.ln_1p()
    }
}

/// ∫₀^`upper` `f` by the composite Simpson rule.
fn simpson(f: impl Fn(f64) -> f64, upper: f64) -> f64 {
    let n = SIMPSON_INTERVALS;
    let h = upper / f64::from(n);
    let mut sum = f(0.0) + f(upper);
    for i in 1..n {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * f(h * f64::from(i));
    }
    sum * h / 3.0
}
