use std::fmt;
use std::sync::Arc;

use crate::Vec3;
use crate::ground::BallState;

/// Air density (kg/m³) on a day of `temperature_c` (°C) at `altitude_m`
/// above sea level: the standard pressure at that altitude over the actual
/// temperature, by the ideal gas law.
///
/// ρ = ρ₀ × (1 − L·h / T₀)^(g·M / (R·L)) × T₀ / T
///
/// where:
/// - ρ₀ = 1.225 kg/m³, T₀ = 288.15 K (ISA sea level)
/// - L = 0.0065 K/m (standard lapse rate), h = `altitude_m` (floored at 0)
/// - g·M / (R·L) = 5.2559 (g = 9.80665 m/s², M = 0.0289644 kg/mol, R =
///   8.31447 J/(mol·K))
/// - T = `temperature_c` + 273.15 K, the temperature at the course, used as
///   given (it is not lowered with altitude)
///
/// The first two factors are ρ₀ times the standard pressure ratio at `h`.
/// On a standard day (15 °C at sea level, 5.1 °C at 1,524 m) this gives the
/// ISA density: 1.225 and 1.0555 kg/m³.
///
/// With a weather station's pressure `p` (Pa) and temperature `T` (K), use
/// the ideal gas law directly: ρ = p / (287.05 · T). Humidity lowers the
/// density by under 1 % in normal weather.
#[must_use]
pub fn standard_air_density(altitude_m: f64, temperature_c: f64) -> f64 {
    const RHO_0: f64 = 1.225;
    const L: f64 = 0.0065; // K/m
    const T0: f64 = 288.15; // K
    const EXPONENT: f64 = 5.2559; // g·M / (R·L)

    let h = altitude_m.max(0.0);
    let pressure_ratio = (1.0 - L * h / T0).powf(EXPONENT);
    RHO_0 * pressure_ratio * (T0 / (temperature_c + 273.15))
}

/// Horizontal wind, in the public frame: `+x` blows downrange (a
/// tailwind), `+z` blows to the right. **Beta.** The vertical component of
/// any wind vector is ignored. [`Tour`](crate::Tour) applies drag and
/// Magnus lift to the ball's velocity relative to the air, in flight and in
/// the hops between bounces; the ground contacts themselves are unchanged,
/// because the turf does not move.
///
/// Every variant is a pure function of the shot time, the ball's position
/// and its height above the ground beneath it.
#[derive(Clone, Default)]
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
    /// A wind field of your own; build one with [`Wind::custom`].
    Custom(Arc<dyn WindField>),
}

/// A wind field: the wind (m/s, public frame) at shot time `t` (s), at the
/// ball's position `pos` (m), `height` (m) above the ground beneath it.
///
/// It must be a pure function of its arguments, with no interior
/// mutability: the integrator evaluates it several times per step and must
/// get the same value for the same arguments. Closures
/// `Fn(f64, Vec3, f64) -> Vec3` implement it. Used through
/// [`Wind::custom`]; the vertical component is ignored.
pub trait WindField: Send + Sync {
    /// The wind at time `t`, position `pos`, `height` above the ground.
    #[must_use]
    fn at(&self, t: f64, pos: Vec3, height: f64) -> Vec3;
}

impl<F: Fn(f64, Vec3, f64) -> Vec3 + Send + Sync> WindField for F {
    fn at(&self, t: f64, pos: Vec3, height: f64) -> Vec3 {
        self(t, pos, height)
    }
}

impl Wind {
    /// A wind of your own, any [`WindField`]. A headwind that starts two
    /// seconds into the shot:
    ///
    /// ```
    /// use dimple::air::Wind;
    /// use dimple::{ShotInput, Tour, Vec3, simulate_flight};
    ///
    /// let gust = Wind::custom(|t: f64, _pos: Vec3, _height: f64| {
    ///     if t > 2.0 { Vec3::new(-8.0, 0.0, 0.0) } else { Vec3::ZERO }
    /// });
    /// let calm = simulate_flight(&ShotInput::driver(), &Tour::SEA_LEVEL);
    /// let gusted = simulate_flight(&ShotInput::driver(), &Tour::new(1.225, gust));
    /// assert!(gusted.carry_m < calm.carry_m);
    /// ```
    #[must_use]
    pub fn custom(field: impl WindField + 'static) -> Self {
        Self::Custom(Arc::new(field))
    }

    /// Horizontal wind velocity (m/s) for a ball at time `t`, position
    /// `pos`, `height` above the ground beneath it. Heights at or below the
    /// ground give zero for [`PowerLaw`](Self::PowerLaw) with a positive
    /// exponent.
    pub(crate) fn at(&self, t: f64, pos: Vec3, height: f64) -> Vec3 {
        let (wind, factor) = match self {
            Self::None => return Vec3::ZERO,
            Self::Uniform(wind) => (*wind, 1.0),
            Self::PowerLaw { at_10m } => {
                let speed = at_10m.x.hypot(at_10m.z);
                let factor = if speed > 0.1 {
                    (height.max(0.0) / 10.0).powf((0.37 - 0.0881 * speed.ln()).max(0.0))
                } else {
                    1.0
                };
                (*at_10m, factor)
            }
            Self::Custom(field) => (field.at(t, pos, height), 1.0),
        };
        Vec3::new(wind.x * factor, 0.0, wind.z * factor)
    }

    /// The ball's velocity relative to the air: its velocity less the wind
    /// at the ball. Still air returns the velocity untouched.
    pub(crate) fn airspeed(&self, ball: &BallState, height: f64) -> Vec3 {
        match self {
            Self::None => ball.vel,
            _ => ball.vel - self.at(ball.time, ball.pos, height),
        }
    }
}

/// Equal variants with equal vectors; two [`Custom`](Wind::Custom) winds
/// are equal only if they share one field (the same `Arc`).
impl PartialEq for Wind {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::None, Self::None) => true,
            (Self::Uniform(a), Self::Uniform(b))
            | (Self::PowerLaw { at_10m: a }, Self::PowerLaw { at_10m: b }) => a == b,
            (Self::Custom(a), Self::Custom(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }
}

impl fmt::Debug for Wind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::None => f.write_str("None"),
            Self::Uniform(w) => f.debug_tuple("Uniform").field(w).finish(),
            Self::PowerLaw { at_10m } => {
                f.debug_struct("PowerLaw").field("at_10m", at_10m).finish()
            }
            Self::Custom(_) => f.debug_tuple("Custom").finish_non_exhaustive(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Wind at `height` for a ball at the origin at t = 0.
    fn at_height(w: &Wind, height: f64) -> Vec3 {
        w.at(0.0, Vec3::ZERO, height)
    }

    #[test]
    fn sea_level_density() {
        let rho = standard_air_density(0.0, 15.0);
        assert!((rho - 1.225).abs() < 1e-12, "{rho}");
    }

    #[test]
    fn isa_density_at_1524_m_on_a_standard_day() {
        // ISA at 5,000 ft (1,524 m): 5.094 °C, 1.0555 kg/m³.
        let rho = standard_air_density(1524.0, 15.0 - 0.0065 * 1524.0);
        assert!((rho - 1.0555).abs() < 0.001, "{rho}");
    }

    #[test]
    fn altitude_at_a_fixed_temperature_scales_with_the_pressure_only() {
        // At one temperature, the density falls exactly as the standard
        // pressure: 0.8320 of sea level at 1,524 m (ISA), so 1.0192 kg/m³ at
        // 15 °C.
        let ratio = standard_air_density(1524.0, 15.0) / standard_air_density(0.0, 15.0);
        assert!((ratio - 0.8320).abs() < 0.0005, "{ratio}");
        let denver = standard_air_density(1609.0, 17.0);
        let reduction = 1.0 - denver / 1.225;
        assert!((0.17..0.19).contains(&reduction), "{reduction}");
    }

    #[test]
    fn hot_day_less_dense() {
        assert!(standard_air_density(0.0, 35.0) < standard_air_density(0.0, 10.0));
    }

    #[test]
    fn below_sea_level_is_floored() {
        assert!((standard_air_density(-100.0, 15.0) - 1.225).abs() < 1e-12);
    }

    #[test]
    fn power_law_profile_and_threshold() {
        let w = Wind::PowerLaw {
            at_10m: Vec3::new(3.0, 99.0, 4.0),
        };
        assert_eq!(at_height(&w, 10.0), Vec3::new(3.0, 0.0, 4.0));
        assert_eq!(at_height(&w, 0.0), Vec3::ZERO);
        assert_eq!(at_height(&w, -1.0), Vec3::ZERO);
        let expected = 0.2_f64.powf(0.37 - 0.0881 * 5.0_f64.ln());
        assert!((at_height(&w, 2.0).x - 3.0 * expected).abs() < 1e-12);
        for speed in [0.0, 0.01, 0.1] {
            let weak = Wind::PowerLaw {
                at_10m: Vec3::new(speed, 9.0, 0.0),
            };
            assert_eq!(at_height(&weak, 0.0), Vec3::new(speed, 0.0, 0.0));
        }
        let above = Wind::PowerLaw {
            at_10m: Vec3::new(0.100_001, 0.0, 0.0),
        };
        assert_eq!(at_height(&above, 0.0), Vec3::ZERO);
        assert_eq!(at_height(&Wind::None, 10.0), Vec3::ZERO);
        assert_eq!(
            at_height(&Wind::Uniform(Vec3::new(3.0, 9.0, -4.0)), 0.0),
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
                assert_eq!(at_height(&w, h), Vec3::new(speed, 0.0, 0.0));
            }
        }
        let below = Wind::PowerLaw {
            at_10m: Vec3::new(66.0, 0.0, 0.0),
        };
        assert!(at_height(&below, 0.0).x.is_finite());
    }

    #[test]
    fn custom_wind_reads_time_position_and_height_and_drops_the_vertical() {
        let w = Wind::custom(|t: f64, pos: Vec3, h: f64| Vec3::new(t, 7.0, pos.x + h));
        assert_eq!(
            w.at(2.0, Vec3::new(5.0, 1.0, 0.0), 0.5),
            Vec3::new(2.0, 0.0, 5.5)
        );
    }

    #[test]
    fn equality_and_debug() {
        let v = Vec3::new(1.0, 0.0, 2.0);
        assert_eq!(Wind::Uniform(v), Wind::Uniform(v));
        assert_ne!(Wind::Uniform(v), Wind::PowerLaw { at_10m: v });
        let c = Wind::custom(move |_: f64, _: Vec3, _: f64| v);
        assert_eq!(c, c.clone());
        assert_ne!(c, Wind::custom(move |_: f64, _: Vec3, _: f64| v));
        assert_ne!(c, Wind::Uniform(v));
        assert_eq!(format!("{c:?}"), "Custom(..)");
        assert_eq!(format!("{:?}", Wind::None), "None");
    }
}
