//! A custom `AeroModel` whose drag follows the Reynolds number, so the air's
//! temperature changes the flight through its viscosity as well as its
//! density. `Tour`'s coefficients depend on the spin ratio alone; this model
//! adds a smooth drag crisis: high drag at low Reynolds number, low drag
//! above it.
//!
//! The same drive and 7-iron fly in air of one density at 0 °C and at 35 °C.
//! The density is held fixed so that only the viscosity differs; in real air
//! pass `standard_air_density(altitude_m, temperature_c)` for each day.
//!
//! `cargo run --release --example viscosity`

use dimple::air::Wind;
use dimple::units::meters_to_yards;
use dimple::{AeroModel, BallState, ShotInput, Tour, Vec3, simulate_flight, standard_air_density};

/// Ball mass (kg), the USGA maximum.
const MASS_KG: f64 = 0.04593;

// Illustrative coefficients in the range of published dimpled-ball wind
// tunnel data, not fitted to anything.
/// Drag at zero spin below the drag crisis.
const CD_SUBCRITICAL: f64 = 0.45;
/// Drag at zero spin above the drag crisis.
const CD_SUPERCRITICAL: f64 = 0.21;
/// Reynolds number at the middle of the drag crisis.
const RE_CRITICAL: f64 = 7.0e4;
/// Width of the drag crisis in Reynolds number.
const RE_WIDTH: f64 = 1.0e4;
/// Drag added per unit of spin ratio.
const CD_SPIN: f64 = 0.35;
/// Lift ceiling.
const CL_MAX: f64 = 0.33;
/// Lift saturation rate in spin ratio.
const SR_SCALE: f64 = 0.15;

/// Dynamic viscosity of air (Pa·s) at `temperature_k`, by Sutherland's law
/// with the constants of the U.S. Standard Atmosphere (1976):
/// μ = μ₀ (T / T₀)^(3/2) (T₀ + S) / (T + S),
/// μ₀ = 1.716e-5 Pa·s at T₀ = 273.15 K, S = 110.4 K.
fn sutherland_viscosity(temperature_k: f64) -> f64 {
    const MU_REF: f64 = 1.716e-5;
    const T_REF: f64 = 273.15;
    const S: f64 = 110.4;
    MU_REF * (temperature_k / T_REF).powf(1.5) * (T_REF + S) / (temperature_k + S)
}

/// A ball whose drag depends on the Reynolds number, in still air of a given
/// density and temperature. Spin decay, radius and inertia come from `Tour`.
struct ReynoldsBall {
    tour: Tour,
    air_density: f64,
    viscosity: f64,
}

impl ReynoldsBall {
    fn new(air_density: f64, temperature_c: f64) -> Self {
        Self {
            tour: Tour::new(air_density, Wind::None),
            air_density,
            viscosity: sutherland_viscosity(temperature_c + 273.15),
        }
    }

    /// Drag and lift coefficients at airspeed `speed` (m/s) and spin rate
    /// `spin_rate` (rad/s).
    fn cd_cl(&self, speed: f64, spin_rate: f64) -> (f64, f64) {
        let radius = self.tour.radius_m();
        let reynolds = self.air_density * speed * 2.0 * radius / self.viscosity;
        let sr = spin_rate * radius / speed;
        // Logistic blend from sub- to supercritical drag across the crisis.
        let crisis = 1.0 / (1.0 + ((reynolds - RE_CRITICAL) / RE_WIDTH).exp());
        let cd = CD_SUPERCRITICAL + (CD_SUBCRITICAL - CD_SUPERCRITICAL) * crisis + CD_SPIN * sr;
        let cl = CL_MAX * (1.0 - (-sr / SR_SCALE).exp());
        (cd, cl)
    }
}

impl AeroModel for ReynoldsBall {
    fn acceleration(&self, ball: &BallState, _height: f64) -> Vec3 {
        let speed = ball.vel.length();
        if speed < 1e-10 {
            return Vec3::ZERO;
        }
        let vel_hat = ball.vel / speed;
        let spin_rate = ball.spin.length();
        let (cd, cl) = self.cd_cl(speed, spin_rate);
        let radius = self.tour.radius_m();
        let area = std::f64::consts::PI * radius * radius;
        let scale = 0.5 * self.air_density * speed * speed * area / MASS_KG;
        let drag = vel_hat * (-cd * scale);
        // Magnus lift along ω̂ × v̂.
        let lift_dir = ball.spin.cross(vel_hat);
        let lift = if lift_dir.length() > 1e-10 {
            lift_dir / lift_dir.length() * (cl * scale)
        } else {
            Vec3::ZERO
        };
        drag + lift
    }

    fn spin_after(&self, ball: &BallState, dt: f64) -> Vec3 {
        self.tour.spin_after(ball, dt)
    }

    fn radius_m(&self) -> f64 {
        self.tour.radius_m()
    }

    fn inertia_factor(&self) -> f64 {
        self.tour.inertia_factor()
    }
}

fn main() {
    let density = standard_air_density(0.0, 15.0);
    println!("air density {density:.3} kg/m³ for every case\n");
    println!(
        "{:<8} {:>14} {:>14} {:>14}",
        "shot", "Reynolds 0 °C", "Reynolds 35 °C", "Tour"
    );
    for (label, input) in [
        ("driver", ShotInput::driver()),
        ("7-iron", ShotInput::seven_iron()),
    ] {
        let carry = |model: &dyn AeroModel| meters_to_yards(simulate_flight(&input, model).carry_m);
        println!(
            "{label:<8} {:>11.1} yd {:>11.1} yd {:>11.1} yd",
            carry(&ReynoldsBall::new(density, 0.0)),
            carry(&ReynoldsBall::new(density, 35.0)),
            carry(&Tour::new(density, Wind::None)),
        );
    }
}
