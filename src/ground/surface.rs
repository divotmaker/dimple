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

/// Penner (2002) eq. 8 reference contact angle θc (rad): 15.4°, his fit to
/// Haake's green impact.
pub(super) const PENNER_THETA_C: f64 = 15.4 * std::f64::consts::PI / 180.0;
/// Impact speed (m/s) of Haake's measured green impact.
pub(super) const PENNER_REF_SPEED: f64 = 18.6;
/// Incidence angle from the normal (rad) of Haake's measured green impact:
/// 44.4°.
pub(super) const PENNER_REF_ANGLE: f64 = 44.4 * std::f64::consts::PI / 180.0;

/// Tangential speed (m/s) below which an impact has no tilt direction, so
/// Penner's contact angle is zero.
const MIN_TILT_SPEED: f64 = 0.01;

/// Normal impact speed (m/s) at and above which Penner's contact angle
/// applies in full. Below it the angle fades in linearly with the normal
/// speed, `θc · v_n / 2 m/s`: the tilt stands for the crater the ball
/// presses into the turf, and a contact this soft (the fall from a bounce
/// under about 20 cm) barely presses one. Lofted landings, from about 4 m/s
/// for a low chip to 10–20 m/s for full shots, get the full angle; the fade
/// acts on the low hops that follow them and on landings launched under
/// about 12°.
const TILT_ONSET_SPEED: f64 = 2.0;

/// Landing angle from the surface (rad) of the shallowest landing Penner
/// (2002) applies eq. 8 to: his low drive, 68.0° from the vertical (22.0°).
const SHALLOW_LANDING: f64 = 22.0 * (std::f64::consts::PI / 180.0);

/// Eq. 8 at that landing (26.8 m/s): 34.0°.
const SHALLOW_TILT: f64 = 34.0 * (std::f64::consts::PI / 180.0);

/// Power of the shallow-landing cap `SHALLOW_TILT · (θ / SHALLOW_LANDING)^p`
/// on the contact angle. Fitted to a simulator's fairway rollouts for
/// landings under 12°, not measured.
const SHALLOW_POWER: f64 = 3.5;

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
/// No law is validated above 3 m/s: the only published speed-dependent fits
/// (Burritt & Henrikson 2024) cover putts, 0–3 m/s. On slopes the published
/// checks favor [`Constant`](Self::Constant): Brede's (1991) slope
/// correction of stimpmeter readings is exact under it and Cross's (2024)
/// speed-independent break angle holds under it alone; at the same stimp
/// the quadratic law also holds a ball on a gentler grade. See
/// `docs/AERO.md`.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum RollLaw {
    /// f = 1: constant rolling resistance, Penner's (2002) run model (his
    /// `ρ_g` = 0.131 reads about 6.0 ft on the stimpmeter).
    Constant,
    /// f = 1 + (v / `v_ref`)²: dimple's own empirical law, not a fit to
    /// published data. [`Surface::GREEN`] and [`Surface::FIRM_GREEN`] use
    /// `v_ref` = 2.5 m/s. `v_ref` (m/s) must be nonzero; ±∞ is
    /// [`Constant`](Self::Constant).
    Quadratic {
        /// Speed scale (m/s).
        v_ref: f64,
    },
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
            Self::Quadratic { v_ref } => {
                let s = speed / v_ref;
                if s.is_finite() { 1.0 + s * s } else { 1.0 }
            }
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
    /// - `Quadratic`: I = `v_ref`²·ln(1 + s₀²)
    /// - `Exponential`: I = 2·`v_ref`²·(1 − e^(−s₀)·(1 + s₀))
    ///
    /// Both parameterized forms tend to v₀² as `v_ref` → ±∞, which s₀ = 0
    /// gives exactly. NaN when `v_ref` is zero or NaN.
    fn stimp_integral(&self, v0: f64) -> f64 {
        let v0_sq = v0 * v0;
        match *self {
            Self::Constant => v0_sq,
            Self::Quadratic { v_ref } => {
                let x = v0_sq / (v_ref * v_ref);
                if x == 0.0 {
                    v0_sq
                } else {
                    v0_sq * x.ln_1p() / x
                }
            }
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
    /// `μ_eff` = `mu_roll` · f(v). [`GREEN`](Self::GREEN) and
    /// [`FIRM_GREEN`](Self::FIRM_GREEN) use
    /// `RollLaw::Quadratic { v_ref: 2.5 }`; [`FAIRWAY`](Self::FAIRWAY) uses
    /// `RollLaw::Constant` with [`plow`](Self::plow).
    pub roll_law: RollLaw,
    /// Turf compliance: dimensionless multiplier on Penner's (2002) eq. 8
    /// effective contact angle, θc = 15.4° · (v / 18.6 m/s) · (φ / 44.4°),
    /// with v the impact speed and φ the incidence angle from the normal,
    /// faded in below 2 m/s of normal impact speed and capped for shallow
    /// landings: below 22° the cap falls to zero at grazing, and above 22°
    /// it grows quickly and binds only on landings faster than Penner's
    /// (see [`impact`](super::impact)). `1.0` is Penner's green; `0.0` is a
    /// rigid flat surface.
    pub compliance: f64,
    /// Torsional friction lever (m): the friction moment about the surface
    /// normal is `pivot` × normal force. It opposes spin about the normal,
    /// both at impact and on the ground, saturating at zero.
    ///
    /// It never changes the velocity of the impact it acts in, nor a slide
    /// or roll on level ground. The path depends on it only through spin
    /// about the normal carried into a later impact, where Penner's tilted
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
    /// no plowing ([`GREEN`](Self::GREEN), [`FIRM_GREEN`](Self::FIRM_GREEN)).
    /// [`FAIRWAY`](Self::FAIRWAY)'s 0.36 is fitted to a simulator's fairway
    /// rollouts for landings at 2–50°, not measured: published work finds
    /// that rolling resistance changes with speed (Burritt & Henrikson 2024)
    /// and gives no law above 3 m/s. It does not act during an impact or a
    /// hop.
    ///
    /// Friction can hold a plowing ball rolling only below
    /// `(1+k)/k · mu_slide · g / plow` (about 38 m/s on `FAIRWAY`, which no
    /// shot reaches); a ball keeps rolling from its first zero slip even
    /// above that speed.
    pub plow: f64,
}

impl Surface {
    /// Standard fairway: Penner's (2002) green impact model (restitution,
    /// friction and compliance 1.0) with a constant rolling resistance and
    /// a plowing rate of 0.36 /s ([`plow`](Self::plow)), stimp 4.8 ft. It
    /// bounces like [`GREEN`](Self::GREEN), and a full shot, which reaches
    /// the ground at 2–4 m/s, rolls shorter on it; a roll faster than about
    /// 11 m/s, such as a thin shot's, runs farther, since `GREEN` keeps its
    /// quadratic law with no plowing fit of its own. The
    /// compliance is chosen, not measured; the run is sensitive to it (see
    /// `docs/AERO.md`).
    pub const FAIRWAY: Self = Self {
        restitution: Restitution::Penner2002,
        mu_slide: 0.40,
        mu_roll: 0.12,
        roll_law: RollLaw::Constant,
        compliance: 1.0,
        pivot: 0.001,
        absorption: 0.0,
        plow: 0.36,
    };

    /// Standard green: Penner's (2002) green impact model with a rolling
    /// resistance of stimp 11.4 ft ([`stimp_feet`](Self::stimp_feet)).
    pub const GREEN: Self = Self {
        restitution: Restitution::Penner2002,
        mu_slide: 0.40,
        mu_roll: 0.055,
        roll_law: RollLaw::Quadratic { v_ref: 2.5 },
        compliance: 1.0,
        pivot: 0.001,
        absorption: 0.0,
        plow: 0.0,
    };

    /// Firm green: [`GREEN`](Self::GREEN) with Penner's (2002, Fig. 8)
    /// sliding friction μ = 0.25, his model of a firm tournament green.
    /// Penner changes only μ; restitution, the contact-angle law and the
    /// rolling resistance stay those of his green, so this has `GREEN`'s
    /// stimp. Penner's high-backspin 9-iron landing here bounces forward
    /// still carrying backspin, then runs back.
    pub const FIRM_GREEN: Self = Self {
        mu_slide: 0.25,
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
    /// Quadratic:    I = v_ref²·ln(1 + s₀²)
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
    /// (k = 2/5), including [`plow`](Self::plow). 11.4 ft for
    /// [`GREEN`](Self::GREEN) and [`FIRM_GREEN`](Self::FIRM_GREEN), 4.8 ft
    /// for [`FAIRWAY`](Self::FAIRWAY).
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
    /// `mu_roll · f(speed)` with f from [`roll_law`](Self::roll_law). Under
    /// the greens' quadratic law a fast-rolling ball decelerates much
    /// harder than a stimpmeter roll. [`plow`](Self::plow) is separate and
    /// not included.
    #[must_use]
    pub fn rolling_friction(&self, speed: f64) -> f64 {
        self.mu_roll * self.roll_law.speed_factor(speed)
    }

    /// Effective contact angle β (rad) for an impact with tangential speed
    /// `v_t` and normal approach speed `v_n` (both ≥ 0).
    ///
    /// `compliance` times the smaller of Penner's eq. 8 θc times the fade
    /// `min(1, v_n / TILT_ONSET_SPEED)` and the shallow-landing cap
    /// `SHALLOW_TILT · (θ / SHALLOW_LANDING)^SHALLOW_POWER`, with θ the
    /// landing angle from the surface; the result is capped at the
    /// incidence angle from the normal, where the ball strikes the
    /// effective plane head-on. At 22° the cap equals eq. 8 at Penner's
    /// shallowest published landing and above it the cap grows out of the
    /// way; below it the tilt falls to zero at grazing, so a fast shallow
    /// landing meets nearly flat ground instead of a steep wall. Zero when
    /// there is (almost) no tangential motion, since the tilt direction is
    /// undefined.
    pub(super) fn contact_angle(&self, v_t: f64, v_n: f64) -> f64 {
        if v_t <= MIN_TILT_SPEED {
            return 0.0;
        }
        let speed = v_t.hypot(v_n);
        let incidence = v_t.atan2(v_n);
        let theta_c = PENNER_THETA_C * (speed / PENNER_REF_SPEED) * (incidence / PENNER_REF_ANGLE);
        let fade = (v_n / TILT_ONSET_SPEED).min(1.0);
        let shallow = SHALLOW_TILT * (v_n.atan2(v_t) / SHALLOW_LANDING).powf(SHALLOW_POWER);
        (self.compliance * (theta_c * fade).min(shallow)).min(incidence)
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
