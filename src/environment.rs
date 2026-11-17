use crate::Vec3;

/// The air a shot flies through: its density, its viscosity and its wind.
///
/// Built from a preset with `..`, e.g.
/// `Environment { temperature_c: 30.0, ..Environment::SEA_LEVEL }`. Every
/// preset is still air with the density from altitude and temperature and
/// the calibrated viscosity. The fields [`measured_air_density`](Self::measured_air_density),
/// [`viscosity`](Self::viscosity) and [`wind`](Self::wind) are **beta**.
#[derive(Debug, Clone, Copy)]
pub struct Environment {
    /// Altitude above sea level (meters). Sets the density unless
    /// [`measured_air_density`](Self::measured_air_density) is given.
    pub altitude_m: f64,
    /// Ambient temperature (Celsius). Sets the density unless
    /// [`measured_air_density`](Self::measured_air_density) is given.
    pub temperature_c: f64,
    /// A measured air density (kg/m³), for callers that know the air better
    /// than an altitude and temperature do (a weather station, or pressure,
    /// temperature and humidity). Used as given: it is not converted to an
    /// altitude or clamped to sea level, so cold dense air flies as dense.
    /// It must be finite and non-negative; `Some(0.0)` gives ballistic
    /// flight, under gravity alone. `None` (every preset): the density comes
    /// from [`altitude_m`](Self::altitude_m) and
    /// [`temperature_c`](Self::temperature_c).
    pub measured_air_density: Option<f64>,
    /// The air's dynamic viscosity (Pa·s), which sets the Reynolds number
    /// and so the drag and lift coefficients. Used as given; it must be
    /// finite and positive. `None` (every preset): dimple's model default,
    /// currently the constant 1.81e-5 Pa·s (air at about 19 °C) that the
    /// ball models were fitted with, whatever the temperature. A later
    /// calibration may derive the default from
    /// [`temperature_c`](Self::temperature_c).
    ///
    /// To follow the temperature, pass
    /// `Some(sutherland_viscosity(env.temperature_c))`
    /// ([`sutherland_viscosity`]). That is **beta, and not calibrated:** the
    /// ball models were fitted at the constant viscosity, and it moves the
    /// one temperature benchmark in `docs/AERO.md` away from the data. See
    /// there before using it.
    pub viscosity: Option<f64>,
    /// The wind, in flight and in the hops between bounces. Every preset is
    /// [`Wind::None`].
    pub wind: Wind,
}

impl Environment {
    /// Sea level, 15°C (ISA standard).
    pub const SEA_LEVEL: Self = Self {
        altitude_m: 0.0,
        temperature_c: 15.0,
        measured_air_density: None,
        viscosity: None,
        wind: Wind::None,
    };

    /// Denver, CO — 5,280 ft (1,609 m), ~17°C summer average.
    pub const DENVER: Self = Self {
        altitude_m: 1609.0,
        temperature_c: 17.0,
        ..Self::SEA_LEVEL
    };

    /// Mexico City — 7,350 ft (2,240 m).
    pub const MEXICO_CITY: Self = Self {
        altitude_m: 2240.0,
        temperature_c: 16.0,
        ..Self::SEA_LEVEL
    };

    /// Air density (kg/m³): [`measured_air_density`](Self::measured_air_density)
    /// if given, otherwise the barometric formula for the troposphere at
    /// [`altitude_m`](Self::altitude_m) and [`temperature_c`](Self::temperature_c):
    ///
    /// ρ = ρ₀ × (1 - L·h / T₀)^(g·M / (R·L))
    ///
    /// where:
    /// - ρ₀ = 1.225 kg/m³ (ISA sea level density at 15°C)
    /// - L = 0.0065 K/m (temperature lapse rate)
    /// - T₀ = 288.15 K (ISA sea level temperature)
    /// - g = 9.80665 m/s²
    /// - M = 0.0289644 kg/mol (molar mass of dry air)
    /// - R = 8.31447 J/(mol·K) (universal gas constant)
    ///
    /// Temperature correction: scales by `T_standard` / `T_actual` at altitude.
    #[must_use]
    pub fn air_density(&self) -> f64 {
        const RHO_0: f64 = 1.225;
        const L: f64 = 0.0065; // K/m
        const T0: f64 = 288.15; // K
        const EXPONENT: f64 = 5.2559; // g·M / (R·L)

        if let Some(rho) = self.measured_air_density {
            return rho;
        }

        let h = self.altitude_m.max(0.0);

        // Standard temperature at this altitude
        let t_standard = T0 - L * h;
        // Actual temperature
        let t_actual = self.temperature_c + 273.15;

        // Barometric formula (pressure-based density)
        let rho_standard = RHO_0 * (1.0 - L * h / T0).powf(EXPONENT);

        // Temperature correction: ideal gas law, ρ ∝ 1/T at constant pressure
        rho_standard * (t_standard / t_actual)
    }

    /// Dynamic viscosity of the air (Pa·s): [`viscosity`](Self::viscosity)
    /// if given, otherwise dimple's model default (currently the calibrated
    /// constant 1.81e-5 Pa·s).
    #[must_use]
    pub fn dynamic_viscosity(&self) -> f64 {
        self.viscosity.unwrap_or(AIR_DYNAMIC_VISCOSITY)
    }
}

impl Default for Environment {
    fn default() -> Self {
        Self::SEA_LEVEL
    }
}

/// Dynamic viscosity of air (Pa·s) when [`Environment::viscosity`] is
/// `None` (dimple's model default, currently this constant): Sutherland's law
/// gives this value at about 19 °C. The `TOUR` and `PATENT_1997` coefficients
/// were fitted with it.
pub(crate) const AIR_DYNAMIC_VISCOSITY: f64 = 1.81e-5;

/// Dynamic viscosity (Pa·s) of dry air at `temperature_c` (°C) by
/// Sutherland's law: μ = b T^(3/2) / (T + S), with T in kelvin,
/// b = 1.458e-6 kg/(m·s·K^½) and S = 110.4 K. That gives 1.72e-5 Pa·s at
/// 0 °C, 1.81e-5 at about 19 °C and 1.88e-5 at 35 °C.
///
/// Equation and constants from *U.S. Standard Atmosphere, 1976*
/// (NOAA, NASA, USAF; NASA TM-X-74335), p. 19:
/// <https://ntrs.nasa.gov/citations/19770009539>
///
/// For [`Environment::viscosity`], which is **beta**:
///
/// ```
/// use dimple::Environment;
/// use dimple::environment::sutherland_viscosity;
///
/// let cold = Environment { temperature_c: 2.0, ..Environment::SEA_LEVEL };
/// let env = Environment {
///     viscosity: Some(sutherland_viscosity(cold.temperature_c)),
///     ..cold
/// };
/// assert!(env.dynamic_viscosity() < 1.81e-5);
/// ```
#[must_use]
pub fn sutherland_viscosity(temperature_c: f64) -> f64 {
    let t = temperature_c + 273.15;
    1.458e-6 * t.powf(1.5) / (t + 110.4)
}

/// Horizontal wind, in the public frame: `+x` blows downrange (a
/// tailwind), `+z` blows to the right. **Beta.** The vertical component of any wind
/// vector is ignored. Drag and Magnus lift act on the ball's velocity
/// relative to the air, in flight and in the hops between bounces; the
/// ground contacts themselves are unchanged, because the turf does not move.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
#[non_exhaustive]
pub enum Wind {
    /// Still air.
    #[default]
    None,
    /// The same wind (m/s) at every height.
    Uniform(Vec3),
    /// A wind measured 10 m above the ground, scaled with height by the
    /// power law `wind(h) = at_10m (h / 10)^α`, where `h` is the ball's
    /// height above the ground directly beneath it. The shear exponent
    /// follows Justus & Mikhail (1976), `α = 0.37 − 0.0881 ln(speed)`, with
    /// the 10 m speed in m/s. At or below 0.1 m/s the wind is applied
    /// uniformly. The exponent is floored at zero, so above about 66.7 m/s
    /// (where it would turn negative) there is no shear and the wind is
    /// applied uniformly as well.
    ///
    /// This is a neutral-stability modeling assumption, not a measured wind
    /// field.
    ///
    /// Justus, C. G. & Mikhail, A. (1976), "Height variation of wind speed
    /// and wind distributions statistics", *Geophysical Research Letters*
    /// 3(5), 261–264. <https://doi.org/10.1029/GL003i005p00261>
    PowerLaw {
        /// Horizontal wind (m/s) at 10 m.
        at_10m: Vec3,
    },
}

impl Wind {
    /// Horizontal wind velocity (m/s) at `height_m` above the ground beneath
    /// the ball. Heights at or below the ground give zero for
    /// [`PowerLaw`](Self::PowerLaw) with a positive exponent.
    #[must_use]
    pub fn at_height(self, height_m: f64) -> Vec3 {
        let (wind, factor) = match self {
            Self::None => return Vec3::ZERO,
            Self::Uniform(wind) => (wind, 1.0),
            Self::PowerLaw { at_10m } => {
                let speed = at_10m.x.hypot(at_10m.z);
                let factor = if speed > 0.1 {
                    (height_m.max(0.0) / 10.0).powf((0.37 - 0.0881 * speed.ln()).max(0.0))
                } else {
                    1.0
                };
                (at_10m, factor)
            }
        };
        Vec3::new(wind.x * factor, 0.0, wind.z * factor)
    }

    /// The ball's velocity relative to the air: `vel` less the wind at the
    /// ball, evaluating the height above ground only when the wind depends
    /// on it. Still air returns `vel` untouched.
    pub(crate) fn airspeed(self, vel: Vec3, height_above_ground: impl FnOnce() -> f64) -> Vec3 {
        match self {
            Self::None => vel,
            Self::Uniform(_) => vel - self.at_height(0.0),
            Self::PowerLaw { .. } => vel - self.at_height(height_above_ground()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sea_level_density() {
        let env = Environment::SEA_LEVEL;
        let rho = env.air_density();
        assert!(
            (rho - 1.225).abs() < 0.001,
            "sea level density should be ~1.225, got {rho:.4}"
        );
    }

    #[test]
    fn denver_less_dense() {
        let sea = Environment::SEA_LEVEL.air_density();
        let denver = Environment::DENVER.air_density();
        let reduction = 1.0 - denver / sea;
        assert!(
            reduction > 0.12 && reduction < 0.25,
            "Denver should be 12-25% less dense, got {:.1}%",
            reduction * 100.0
        );
    }

    #[test]
    fn mexico_city_less_than_denver() {
        let denver = Environment::DENVER.air_density();
        let mexico = Environment::MEXICO_CITY.air_density();
        assert!(
            mexico < denver,
            "Mexico City ({mexico:.4}) should be less dense than Denver ({denver:.4})"
        );
    }

    #[test]
    fn hot_day_less_dense() {
        let cool = Environment {
            temperature_c: 10.0,
            ..Environment::SEA_LEVEL
        };
        let hot = Environment {
            temperature_c: 35.0,
            ..Environment::SEA_LEVEL
        };
        assert!(
            hot.air_density() < cool.air_density(),
            "hot air should be less dense"
        );
    }

    #[test]
    #[allow(clippy::float_cmp, reason = "the density is used as given")]
    fn measured_density_is_used_as_given() {
        for rho in [0.0, 0.8, 1.225, 1.5] {
            let env = Environment {
                altitude_m: 3000.0,
                temperature_c: -20.0,
                measured_air_density: Some(rho),
                ..Environment::SEA_LEVEL
            };
            assert_eq!(env.air_density(), rho);
        }
    }

    #[test]
    #[allow(clippy::float_cmp, reason = "the viscosity is returned as is")]
    fn viscosity_given_or_calibrated_and_sutherland_reference_values() {
        assert!((sutherland_viscosity(15.0) - 1.7894e-5).abs() < 1e-9);
        assert!((sutherland_viscosity(0.0) - 1.716e-5).abs() < 1e-9);
        assert!(sutherland_viscosity(35.0) > sutherland_viscosity(0.0));
        // `None` is the calibrated constant, whatever the temperature; a
        // given viscosity is used as is.
        for t in [0.0, 15.0, 35.0] {
            let env = Environment {
                temperature_c: t,
                ..Environment::SEA_LEVEL
            };
            assert_eq!(env.dynamic_viscosity(), AIR_DYNAMIC_VISCOSITY);
            let given = Environment {
                viscosity: Some(2.0e-5),
                ..env
            };
            assert_eq!(given.dynamic_viscosity(), 2.0e-5);
        }
    }

    #[test]
    fn power_law_profile_and_threshold() {
        let w = Wind::PowerLaw {
            at_10m: Vec3::new(3.0, 99.0, 4.0),
        };
        assert_eq!(w.at_height(10.0), Vec3::new(3.0, 0.0, 4.0));
        assert_eq!(w.at_height(0.0), Vec3::ZERO);
        assert_eq!(w.at_height(-1.0), Vec3::ZERO);
        let expected = 0.2_f64.powf(0.37 - 0.0881 * 5.0_f64.ln());
        assert!((w.at_height(2.0).x - 3.0 * expected).abs() < 1e-12);
        for speed in [0.0, 0.01, 0.1] {
            let weak = Wind::PowerLaw {
                at_10m: Vec3::new(speed, 9.0, 0.0),
            };
            assert_eq!(weak.at_height(0.0), Vec3::new(speed, 0.0, 0.0));
        }
        let above = Wind::PowerLaw {
            at_10m: Vec3::new(0.100_001, 0.0, 0.0),
        };
        assert_eq!(above.at_height(0.0), Vec3::ZERO);
        assert_eq!(Wind::None.at_height(10.0), Vec3::ZERO);
        assert_eq!(
            Wind::Uniform(Vec3::new(3.0, 9.0, -4.0)).at_height(0.0),
            Vec3::new(3.0, 0.0, -4.0)
        );
    }

    #[test]
    fn power_law_exponent_never_goes_negative() {
        // Above e^(0.37 / 0.0881) = 66.7 m/s the fitted exponent is negative;
        // the wind is then uniform, never infinite at the ground.
        for speed in [67.0, 100.0, 1000.0] {
            let w = Wind::PowerLaw {
                at_10m: Vec3::new(speed, 0.0, 0.0),
            };
            for h in [0.0, 0.5, 10.0, 80.0] {
                assert_eq!(w.at_height(h), Vec3::new(speed, 0.0, 0.0));
            }
        }
        let below = Wind::PowerLaw {
            at_10m: Vec3::new(66.0, 0.0, 0.0),
        };
        assert!(below.at_height(0.0).x.is_finite());
    }
}
