use crate::Vec3;

/// Atmospheric environment for ball flight simulation.
#[derive(Debug, Clone, Copy)]
pub struct Environment {
    /// Altitude above sea level (meters).
    pub altitude_m: f64,
    /// Ambient temperature (Celsius).
    pub temperature_c: f64,
}

impl Environment {
    /// Sea level, 15°C (ISA standard).
    pub const SEA_LEVEL: Self = Self {
        altitude_m: 0.0,
        temperature_c: 15.0,
    };

    /// Denver, CO — 5,280 ft (1,609 m), ~17°C summer average.
    pub const DENVER: Self = Self {
        altitude_m: 1609.0,
        temperature_c: 17.0,
    };

    /// Mexico City — 7,350 ft (2,240 m).
    pub const MEXICO_CITY: Self = Self {
        altitude_m: 2240.0,
        temperature_c: 16.0,
    };

    /// Air density (kg/m³) using the barometric formula for the troposphere.
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
}

impl Default for Environment {
    fn default() -> Self {
        Self::SEA_LEVEL
    }
}

/// Dynamic viscosity of air (Pa·s) used wherever no temperature is given:
/// Sutherland's law ([`dynamic_viscosity`]) gives this value at about
/// 19 °C. The `TOUR` and `PATENT_1997` coefficients were fitted with it.
pub(crate) const AIR_DYNAMIC_VISCOSITY: f64 = 1.81e-5;

/// Air properties given directly, for callers that already know the density
/// (from a weather station, or from pressure, temperature and humidity)
/// instead of an altitude and temperature as in [`Environment`].
#[derive(Debug, Clone, Copy)]
pub struct Atmosphere {
    /// Air density (kg/m³). [`from_density`](Self::from_density) validates it.
    pub air_density: f64,
    /// Dynamic viscosity of the air (Pa·s), which sets the Reynolds number
    /// and so the drag and lift coefficients.
    pub dynamic_viscosity: f64,
}

impl Atmosphere {
    /// `air_density` (kg/m³) with the constant viscosity used by
    /// [`Environment`]-based entry points, so that results through those
    /// entry points stay unchanged. Not validated: a zero density gives
    /// ballistic flight, under gravity alone.
    #[must_use]
    pub const fn legacy(air_density: f64) -> Self {
        Self {
            air_density,
            dynamic_viscosity: AIR_DYNAMIC_VISCOSITY,
        }
    }

    /// A measured `air_density` (kg/m³), with the constant viscosity of
    /// [`legacy`](Self::legacy). `None` unless the density is finite and
    /// positive. The density is used as given: it is not converted to an
    /// altitude or clamped to sea level.
    #[must_use]
    pub fn from_density(air_density: f64) -> Option<Self> {
        (air_density.is_finite() && air_density > 0.0).then_some(Self::legacy(air_density))
    }

    /// A measured `air_density` (kg/m³) with the viscosity of dry air at
    /// `temperature_c` (°C) from [`dynamic_viscosity`]. `None` if either is
    /// invalid.
    #[must_use]
    pub fn from_density_and_temperature(air_density: f64, temperature_c: f64) -> Option<Self> {
        let mut atmosphere = Self::from_density(air_density)?;
        atmosphere.dynamic_viscosity = dynamic_viscosity(temperature_c)?;
        Some(atmosphere)
    }
}

/// Dynamic viscosity of dry air (Pa·s) at `temperature_c` (°C), by
/// Sutherland's law: μ = b T^(3/2) / (T + S), with T in kelvin,
/// b = 1.458e-6 kg/(m·s·K^½) and S = 110.4 K. `None` at or below absolute
/// zero, or for a non-finite temperature.
///
/// Constants from Cambridge University Engineering Department, *Properties
/// of Fluids*, "Formulas for viscosity":
/// <https://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/fprops/propsoffluids/node5.html>
#[must_use]
pub fn dynamic_viscosity(temperature_c: f64) -> Option<f64> {
    let t = temperature_c + 273.15;
    if !t.is_finite() || t <= 0.0 {
        return None;
    }
    let mu = 1.458e-6 * t.powf(1.5) / (t + 110.4);
    (mu.is_finite() && mu > 0.0).then_some(mu)
}

/// Horizontal wind sampled at each RK4 stage, including airborne ground hops.
#[derive(Debug, Clone, Copy)]
pub enum WindProfile {
    /// The same wind at every height. Its vertical component is ignored.
    Uniform(Vec3),
    /// A wind measured 10 m above the ground, scaled with height by the
    /// power law `wind(h) = wind10 (h / 10)^α`, where `h` is the ball's
    /// height above the ground directly beneath it. The shear exponent
    /// follows Justus & Mikhail (1976), `α = 0.37 − 0.0881 ln(speed)`, with
    /// the 10 m speed in m/s. At or below 0.1 m/s the logarithm is
    /// meaningless and the wind is applied uniformly instead.
    ///
    /// This is a neutral-stability modeling assumption, not a measured wind
    /// field. Its vertical component is ignored.
    ///
    /// Justus, C. G. & Mikhail, A. (1976), "Height variation of wind speed
    /// and wind distributions statistics", *Geophysical Research Letters*
    /// 3(5), 261–264. <https://doi.org/10.1029/GL003i005p00261>
    PowerLaw {
        /// Horizontal wind (m/s) at 10 m.
        wind10: Vec3,
    },
}

impl WindProfile {
    /// Horizontal wind velocity (m/s) at `height_m` above the ground beneath
    /// the ball. Heights at or below the ground give zero for
    /// [`PowerLaw`](Self::PowerLaw).
    #[must_use]
    pub fn at_height(self, height_m: f64) -> Vec3 {
        let (wind, factor) = match self {
            Self::Uniform(wind) => (wind, 1.0),
            Self::PowerLaw { wind10 } => {
                let speed = wind10.x.hypot(wind10.z);
                let factor = if speed > 0.1 {
                    (height_m.max(0.0) / 10.0).powf(0.37 - 0.0881 * speed.ln())
                } else {
                    1.0
                };
                (wind10, factor)
            }
        };
        Vec3::new(wind.x * factor, 0.0, wind.z * factor)
    }

    /// [`at_height`](Self::at_height), evaluating the height above ground
    /// only when the profile depends on it, so uniform wind costs no terrain
    /// lookups.
    pub(crate) fn sample(self, height_above_ground: impl FnOnce() -> f64) -> Vec3 {
        match self {
            Self::Uniform(_) => self.at_height(0.0),
            Self::PowerLaw { .. } => self.at_height(height_above_ground()),
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
            altitude_m: 0.0,
            temperature_c: 10.0,
        };
        let hot = Environment {
            altitude_m: 0.0,
            temperature_c: 35.0,
        };
        assert!(
            hot.air_density() < cool.air_density(),
            "hot air should be less dense"
        );
    }

    #[test]
    #[allow(clippy::float_cmp, reason = "the density is stored as given")]
    fn density_is_used_as_given() {
        for rho in [0.8, 1.225, 1.5] {
            assert_eq!(Atmosphere::from_density(rho).unwrap().air_density, rho);
        }
        for rho in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(Atmosphere::from_density(rho).is_none());
        }
    }

    #[test]
    fn sutherland_reference_values() {
        assert!((dynamic_viscosity(15.0).unwrap() - 1.7894e-5).abs() < 1e-9);
        assert!((dynamic_viscosity(0.0).unwrap() - 1.716e-5).abs() < 1e-9);
        assert!(dynamic_viscosity(35.0).unwrap() > dynamic_viscosity(0.0).unwrap());
        for t in [-273.15, -300.0, f64::NAN, f64::INFINITY] {
            assert!(dynamic_viscosity(t).is_none());
        }
        assert!(Atmosphere::from_density_and_temperature(1.2, f64::NAN).is_none());
        assert!(Atmosphere::from_density_and_temperature(0.0, 20.0).is_none());
    }

    #[test]
    fn power_law_profile_and_threshold() {
        use Vec3;
        let w = WindProfile::PowerLaw {
            wind10: Vec3::new(3.0, 99.0, 4.0),
        };
        assert_eq!(w.at_height(10.0), Vec3::new(3.0, 0.0, 4.0));
        assert_eq!(w.at_height(0.0), Vec3::ZERO);
        assert_eq!(w.at_height(-1.0), Vec3::ZERO);
        let expected = 0.2_f64.powf(0.37 - 0.0881 * 5.0_f64.ln());
        assert!((w.at_height(2.0).x - 3.0 * expected).abs() < 1e-12);
        for speed in [0.0, 0.01, 0.1] {
            let weak = WindProfile::PowerLaw {
                wind10: Vec3::new(speed, 9.0, 0.0),
            };
            assert_eq!(weak.at_height(0.0), Vec3::new(speed, 0.0, 0.0));
        }
        let above = WindProfile::PowerLaw {
            wind10: Vec3::new(0.100_001, 0.0, 0.0),
        };
        assert_eq!(above.at_height(0.0), Vec3::ZERO);
    }
}
