use crate::Vec3;

use crate::air::Wind;
use crate::ground::BallState;
use crate::units;

/// The air's effect on one ball: the right-hand side of the flight
/// equations, less gravity.
///
/// dimple integrates the flight and the hops between bounces with RK4 and
/// asks the model for the acceleration at every stage; it adds gravity
/// itself, along −Y with magnitude [`gravity`](Self::gravity). The ground
/// reads [`radius_m`](Self::radius_m), [`inertia_factor`](Self::inertia_factor)
/// and [`gravity`](Self::gravity) for the contacts. [`Tour`] is dimple's
/// model; implement the trait for your own ball, air or wind.
///
/// Every method is a pure function of its arguments: the integrator calls
/// [`acceleration`](Self::acceleration) several times per step and must
/// get the same value for the same arguments, so a model has no interior
/// mutability. A model that changes with time (a gust) reads
/// [`BallState::time`], the absolute shot time.
///
/// A ball with drag only, `a = −k·|v|·v`:
///
/// ```
/// use dimple::{AeroModel, BallState, ShotInput, Surface, Vec3, simulate_shot};
///
/// struct DragOnly {
///     k: f64,
/// }
///
/// impl AeroModel for DragOnly {
///     fn acceleration(&self, ball: &BallState, _height: f64) -> Vec3 {
///         ball.vel * (-self.k * ball.vel.length())
///     }
///     fn spin_after(&self, ball: &BallState, _dt: f64) -> Vec3 {
///         ball.spin
///     }
///     fn radius_m(&self) -> f64 {
///         0.021_335
///     }
///     fn inertia_factor(&self) -> f64 {
///         0.4
///     }
/// }
///
/// let shot = simulate_shot(&ShotInput::driver(), &DragOnly { k: 0.004 }, &Surface::FAIRWAY);
/// assert!(shot.flight.carry_m > 100.0);
/// ```
///
/// Another gravity, consistent in flight and on the ground: wrap [`Tour`]
/// and override [`gravity`](Self::gravity).
///
/// ```
/// use dimple::{AeroModel, BallState, ShotInput, Surface, Tour, Vec3, simulate_shot};
///
/// /// The tour ball under a dome on the Moon, in sea-level air.
/// struct MoonDome(Tour);
///
/// impl AeroModel for MoonDome {
///     fn acceleration(&self, ball: &BallState, height: f64) -> Vec3 {
///         self.0.acceleration(ball, height)
///     }
///     fn spin_after(&self, ball: &BallState, dt: f64) -> Vec3 {
///         self.0.spin_after(ball, dt)
///     }
///     fn radius_m(&self) -> f64 {
///         self.0.radius_m()
///     }
///     fn inertia_factor(&self) -> f64 {
///         self.0.inertia_factor()
///     }
///     fn gravity(&self) -> f64 {
///         1.62
///     }
/// }
///
/// let input = ShotInput::seven_iron();
/// let earth = simulate_shot(&input, &Tour::SEA_LEVEL, &Surface::FAIRWAY);
/// let moon = simulate_shot(&input, &MoonDome(Tour::SEA_LEVEL), &Surface::FAIRWAY);
/// assert!(moon.flight.carry_m > earth.flight.carry_m);
/// ```
pub trait AeroModel {
    /// Acceleration (m/s²) the air gives the ball: drag, lift and anything
    /// else the model has. Gravity is not part of it. `height` is the
    /// ball's height (m) above the ground beneath it.
    #[must_use]
    fn acceleration(&self, ball: &BallState, height: f64) -> Vec3;

    /// The ball's spin (rad/s) after `dt` seconds of flight from `ball`.
    /// dimple calls it once per flight step; the hops between bounces hold
    /// the spin fixed and do not call it.
    #[must_use]
    fn spin_after(&self, ball: &BallState, dt: f64) -> Vec3;

    /// Ball radius (m).
    #[must_use]
    fn radius_m(&self) -> f64;

    /// Moment-of-inertia factor k, with I = k·m·R². `0.4` (2/5) is a
    /// uniform solid sphere. Used on the ground, where it sets how friction
    /// splits between slowing the ball and changing its spin: a ball sliding
    /// on flat ground without spin leaves rolling at 1/(1+k) of its speed.
    #[must_use]
    fn inertia_factor(&self) -> f64;

    /// Magnitude (m/s²) of gravity, which acts along −Y, in flight and on
    /// the ground. Standard gravity, [`units::G`], unless overridden.
    #[must_use]
    fn gravity(&self) -> f64 {
        units::G
    }
}

/// dimple's model of a modern tour ball in a given air: its density and its
/// wind.
///
/// Drag and Magnus lift act on the ball's velocity relative to the air, and
/// spin decays at 4 % per second. The drag and lift coefficients are fitted
/// constants (see `docs/AERO.md` in the repository for the formulas, the
/// values, the data they were fitted to and the fit quality); they are
/// **beta** and may change between releases.
///
/// ```
/// use dimple::air::Wind;
/// use dimple::{ShotInput, Tour, Vec3, simulate_flight, standard_air_density};
///
/// // 1,500 m up on a 25 °C day, into a 4 m/s headwind.
/// let tour = Tour::new(
///     standard_air_density(1500.0, 25.0),
///     Wind::Uniform(Vec3::new(-4.0, 0.0, 0.0)),
/// );
/// let flight = simulate_flight(&ShotInput::driver(), &tour);
/// assert!(flight.carry_m > 200.0);
/// ```
#[derive(Debug, Clone)]
pub struct Tour {
    air_density: f64,
    wind: Wind,
}

/// Compile-time check: models and winds can be shared across threads.
const _: () = {
    const fn send_sync<T: Send + Sync>() {}
    send_sync::<Tour>();
    send_sync::<Wind>();
};

/// Ball mass (kg), the USGA maximum (1.62 oz).
const MASS_KG: f64 = 0.04593;
/// Ball radius (m): half the USGA minimum diameter (1.68 in).
const RADIUS_M: f64 = 0.04267 / 2.0;
/// Moment-of-inertia factor of a uniform solid sphere.
const INERTIA_FACTOR: f64 = 0.4;
/// Drag at zero spin.
const CD_BASE: f64 = 0.10;
/// Slope of the spin drag.
const CD_SPIN: f64 = 1.69;
/// Spin ratio at which the spin drag reaches half its cap.
const SR_SAT: f64 = 0.24;
/// Lift ceiling.
const CL_MAX: f64 = 0.359;
/// Lift saturation rate in spin ratio.
const SR_SCALE: f64 = 0.149;
/// Height of the low-spin lift term.
const CL_BUMP: f64 = 0.0185;
/// Spin ratio at which the low-spin lift term peaks.
const SR_BUMP: f64 = 0.027;
/// Spin decay rate (s⁻¹).
const SPIN_DECAY: f64 = 0.04;

impl Tour {
    /// Sea level on a standard day (1.225 kg/m³, 15 °C), still air.
    pub const SEA_LEVEL: Self = Self::new(1.225, Wind::None);

    /// The tour ball in air of density `air_density` (kg/m³) with `wind`.
    ///
    /// The density is used as given, so cold dense air flies as dense. It
    /// must be finite and non-negative; `0.0` is no air at all, under
    /// gravity alone. [`standard_air_density`](crate::air::standard_air_density)
    /// gives it from an altitude and a temperature.
    #[must_use]
    pub const fn new(air_density: f64, wind: Wind) -> Self {
        Self { air_density, wind }
    }
}

impl AeroModel for Tour {
    fn acceleration(&self, ball: &BallState, height: f64) -> Vec3 {
        let airspeed = self.wind.airspeed(ball, height);
        air_acceleration(airspeed, ball.spin, self.air_density)
    }

    fn spin_after(&self, ball: &BallState, dt: f64) -> Vec3 {
        ball.spin * (-SPIN_DECAY * dt).exp()
    }

    fn radius_m(&self) -> f64 {
        RADIUS_M
    }

    fn inertia_factor(&self) -> f64 {
        INERTIA_FACTOR
    }
}

/// Drag and lift coefficients of the tour ball at angular speed
/// `spin_rate` (rad/s) and airspeed `speed` (m/s), from the spin ratio
/// SR = ω·R / V:
///
/// - Cd = `CD_BASE` + `CD_SPIN`·SR / (1 + SR / `SR_SAT`)
/// - Cl = `CL_MAX`·(1 − exp(−SR / `SR_SCALE`)) + `CL_BUMP`·u·exp(1 − u),
///   u = (SR / `SR_BUMP`)⁴
fn cd_cl(spin_rate: f64, speed: f64) -> (f64, f64) {
    if speed < 1e-10 {
        return (CD_BASE, 0.0);
    }
    let sr = (spin_rate * RADIUS_M / speed).abs();
    let cd = CD_BASE + CD_SPIN * sr / (1.0 + sr / SR_SAT);
    let u = (sr / SR_BUMP).powi(4);
    let cl = CL_MAX * (1.0 - (-sr / SR_SCALE).exp()) + CL_BUMP * u * (1.0 - u).exp();
    (cd, cl)
}

/// Drag and Magnus lift (m/s²) on the tour ball moving at `vel` relative to
/// air of density `air_density`, with angular velocity `spin`. Drag opposes
/// `vel`; lift is along ω̂ × v̂.
fn air_acceleration(vel: Vec3, spin: Vec3, air_density: f64) -> Vec3 {
    let speed = vel.length();
    if speed < 1e-10 {
        return Vec3::ZERO;
    }
    let vel_hat = vel / speed;
    let rate = spin.length();
    let (cd, cl) = cd_cl(rate, speed);
    // Dynamic pressure times area over mass.
    let area = std::f64::consts::PI * RADIUS_M * RADIUS_M;
    let scale = 0.5 * air_density * speed * speed * area / MASS_KG;
    let drag = vel_hat * (-cd * scale);
    let magnus = if rate > 1e-10 {
        let lift_dir = (spin / rate).cross(vel_hat);
        let lift_mag = lift_dir.length();
        if lift_mag > 1e-10 {
            (lift_dir / lift_mag) * (cl * scale)
        } else {
            Vec3::ZERO
        }
    } else {
        Vec3::ZERO
    };
    drag + magnus
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::air::standard_air_density;
    use crate::units::rpm_to_rads;

    fn moving(vel: Vec3, spin: Vec3) -> BallState {
        BallState {
            time: 0.0,
            pos: Vec3::ZERO,
            vel,
            spin,
        }
    }

    /// Backspin and sidespin (rpm) for a straight launch.
    fn spin(backspin_rpm: f64, sidespin_rpm: f64) -> Vec3 {
        Vec3::new(0.0, -rpm_to_rads(sidespin_rpm), rpm_to_rads(backspin_rpm))
    }

    fn accel(vel: Vec3, spin: Vec3) -> Vec3 {
        Tour::SEA_LEVEL.acceleration(&moving(vel, spin), 0.0)
    }

    /// Spin rate (rad/s) for spin ratio `sr` at `speed`.
    fn sr_to_spin(sr: f64, speed: f64) -> f64 {
        sr * speed / RADIUS_M
    }

    // ══════════════════════════════════════════════════════════════════
    //  FORCE MODEL — structural invariants
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn no_spin_gives_drag_only() {
        let a = accel(Vec3::new(50.0, 0.0, 0.0), Vec3::ZERO);
        assert!(a.x < 0.0, "drag should oppose +X velocity");
        assert!(a.y.abs() < 1e-12, "no lift without spin");
        assert!(a.z.abs() < 1e-12, "no lateral force without spin");
    }

    #[test]
    fn positive_sidespin_pushes_right() {
        let a = accel(Vec3::new(50.0, 0.0, 0.0), spin(3000.0, 1000.0));
        assert!(a.z > 0.0, "positive sidespin should push toward +Z: {a:?}");
    }

    #[test]
    fn negative_sidespin_pushes_left() {
        let a = accel(Vec3::new(50.0, 0.0, 0.0), spin(3000.0, -1000.0));
        assert!(a.z < 0.0, "negative sidespin should push toward -Z: {a:?}");
    }

    #[test]
    fn backspin_lifts() {
        let a = accel(Vec3::new(50.0, 0.0, 0.0), spin(3000.0, 0.0));
        assert!(a.y > 0.0, "backspin should add upward force: {a:?}");
    }

    #[test]
    #[allow(clippy::float_cmp, reason = "no air gives exactly nothing")]
    fn no_air_no_force() {
        let vacuum = Tour::new(0.0, Wind::None);
        let a = vacuum.acceleration(&moving(Vec3::new(50.0, 3.0, 1.0), spin(3000.0, 500.0)), 1.0);
        assert_eq!(a, Vec3::ZERO);
    }

    #[test]
    fn less_drag_in_thinner_air() {
        let ball = moving(Vec3::new(70.0, 10.0, 0.0), spin(2700.0, 0.0));
        let sea = Tour::SEA_LEVEL.acceleration(&ball, 0.0);
        let high =
            Tour::new(standard_air_density(1609.0, 17.0), Wind::None).acceleration(&ball, 0.0);
        assert!(high.x > sea.x, "less drag at altitude: {sea:?} {high:?}");
    }

    #[test]
    #[allow(clippy::float_cmp, reason = "the default is the constant itself")]
    fn tour_uses_standard_gravity() {
        assert_eq!(Tour::SEA_LEVEL.gravity(), units::G);
    }

    // ══════════════════════════════════════════════════════════════════
    //  Cd / Cl
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn cl_zero_without_spin() {
        let (_, cl) = cd_cl(0.0, 70.0);
        assert!(cl.abs() < 1e-10, "Cl should be 0 without spin, got {cl:.6}");
    }

    #[test]
    fn cl_increases_with_spin() {
        let (_, cl_high) = cd_cl(400.0, 70.0);
        let (_, cl_low) = cd_cl(200.0, 70.0);
        assert!(cl_high > cl_low, "{cl_high:.3} vs {cl_low:.3}");
    }

    #[test]
    fn cl_rises_with_spin_ratio_through_the_low_spin_term() {
        // The low-spin term must not make lift fall with spin anywhere
        // below SR 0.06.
        let speed = 70.0;
        let mut prev = 0.0;
        for i in 1..=600 {
            let sr = f64::from(i) * 1e-4;
            let (_, cl) = cd_cl(sr_to_spin(sr, speed), speed);
            assert!(cl > prev, "Cl falls at SR {sr}: {cl} after {prev}");
            prev = cl;
        }
    }

    /// Cl without the low-spin term.
    fn cl_without_bump(sr: f64) -> f64 {
        CL_MAX * (1.0 - (-sr / SR_SCALE).exp())
    }

    #[test]
    fn low_spin_term_is_zero_at_zero_spin_and_smooth() {
        let speed = 60.0;
        let (_, at_zero) = cd_cl(0.0, speed);
        assert!(at_zero.abs() < 1e-15);
        // The term is flat at zero (u⁴ starts with zero slope) and has no
        // step anywhere: its slope stays under 3 per unit of spin ratio.
        let bump = |sr: f64| cd_cl(sr_to_spin(sr, speed), speed).1 - cl_without_bump(sr);
        assert!(bump(1e-5) < 1e-14, "{}", bump(1e-5));
        let h = 1e-6;
        let mut sr = h;
        while sr < 0.06 {
            let slope = (bump(sr + h) - bump(sr)).abs() / h;
            assert!(slope < 3.0, "slope {slope} at SR {sr}");
            sr += 1e-4;
        }
        // It peaks at SR_BUMP with height CL_BUMP.
        assert!((bump(SR_BUMP) - CL_BUMP).abs() < 1e-12);
    }

    #[test]
    fn low_spin_term_fades_out_above_sr_0_045() {
        // At SR 0.045 the term adds 1.7e-4 (0.2 % of Cl); from SR 0.06 on,
        // under 1e-9.
        let speed = 60.0;
        for i in 0..=500 {
            let sr = 0.045 + f64::from(i) * 1e-3;
            let (_, cl) = cd_cl(sr_to_spin(sr, speed), speed);
            let diff = (cl - cl_without_bump(sr)).abs();
            let bound = if sr >= 0.06 { 1e-9 } else { 2e-4 };
            assert!(diff < bound, "SR {sr}: differs by {diff}");
        }
    }

    #[test]
    fn spin_drag_saturates() {
        // Doubling SR from 0.25 to 0.5 adds less drag than the first 0.25.
        let speed = 70.0;
        let (cd0, _) = cd_cl(0.0, speed);
        let (cd1, _) = cd_cl(sr_to_spin(0.25, speed), speed);
        let (cd2, _) = cd_cl(sr_to_spin(0.50, speed), speed);
        let (first, second) = (cd1 - cd0, cd2 - cd1);
        assert!(
            second < first * 0.75,
            "first ΔCd={first:.3}, second ΔCd={second:.3}"
        );
    }

    #[test]
    fn cd_increases_with_spin() {
        let (cd_high, _) = cd_cl(400.0, 70.0);
        let (cd_low, _) = cd_cl(200.0, 70.0);
        assert!(cd_high > cd_low, "{cd_high:.3} vs {cd_low:.3}");
    }

    #[test]
    fn coefficients_do_not_depend_on_speed_at_a_given_spin_ratio() {
        for sr in [0.0, 0.03, 0.08, 0.3] {
            let slow = cd_cl(sr_to_spin(sr, 20.0), 20.0);
            let fast = cd_cl(sr_to_spin(sr, 80.0), 80.0);
            assert!((slow.0 - fast.0).abs() < 1e-12 && (slow.1 - fast.1).abs() < 1e-12);
        }
    }

    #[test]
    fn driver_and_iron_coefficients() {
        // Driver, SR 0.074: Cd 0.196, Cl 0.141. 7-iron, SR 0.29: Cd 0.32,
        // Cl 0.31 (docs/AERO.md).
        let (cd, cl) = cd_cl(sr_to_spin(0.074, 70.0), 70.0);
        assert!(
            (cd - 0.196).abs() < 0.001 && (cl - 0.141).abs() < 0.001,
            "{cd} {cl}"
        );
        let (cd, cl) = cd_cl(sr_to_spin(0.29, 50.0), 50.0);
        assert!(
            (cd - 0.32).abs() < 0.005 && (cl - 0.31).abs() < 0.005,
            "{cd} {cl}"
        );
    }

    #[test]
    fn tour_cl_near_li_2017() {
        // Li et al. (2017): Cl = 0.135 at Re 1.1e5, SR 0.1 (LES on one
        // production ball). The fitted model gives about 0.18 there.
        let (_, cl) = cd_cl(sr_to_spin(0.1, 50.0), 50.0);
        assert!((cl - 0.135).abs() < 0.08, "Cl at SR 0.1: {cl:.3}");
    }

    // ══════════════════════════════════════════════════════════════════
    //  SPIN DECAY
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn spin_decays_exponentially_about_a_fixed_axis() {
        let start = moving(Vec3::new(60.0, 0.0, 0.0), spin(5000.0, 800.0));
        let after = Tour::SEA_LEVEL.spin_after(&start, 1.0);
        let expected = start.spin * (-0.04_f64).exp();
        assert!((after - expected).length() < 1e-9);
        assert!(after.normalize().dot(start.spin.normalize()) > 1.0 - 1e-12);
    }

    #[test]
    fn derived_geometry() {
        assert!((Tour::SEA_LEVEL.radius_m() - 0.04267 / 2.0).abs() < 1e-15);
        assert!((Tour::SEA_LEVEL.inertia_factor() - 0.4).abs() < 1e-15);
    }
}
