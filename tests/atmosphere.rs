//! Measured air, viscosity and wind profile through the public API.

use dimple::environment::{Viscosity, Wind};
use dimple::ground::{BallState, GroundModel};
use dimple::{BallModel, Environment, ShotInput, Surface, Vec3, simulate_shot};

/// Measured air of `air_density` (kg/m³) at `temperature_c`, with the
/// viscosity from Sutherland's law.
fn measured(air_density: f64, temperature_c: f64) -> Environment {
    Environment {
        temperature_c,
        measured_air_density: Some(air_density),
        viscosity: Viscosity::Sutherland,
        ..Environment::SEA_LEVEL
    }
}

#[test]
#[allow(
    clippy::float_cmp,
    reason = "a measured density equal to the computed one is the same computation"
)]
fn measured_density_equal_to_the_computed_one_gives_the_same_shot() {
    for env in [Environment::SEA_LEVEL, Environment::DENVER] {
        for wind in [Wind::None, Wind::Uniform(Vec3::new(-5.0, 0.0, 3.0))] {
            let env = Environment { wind, ..env };
            let given = Environment {
                altitude_m: 0.0,
                temperature_c: 40.0,
                measured_air_density: Some(env.air_density()),
                ..env
            };
            let input = ShotInput::driver();
            let a = simulate_shot(&input, &BallModel::TOUR, &env, &Surface::FAIRWAY);
            let b = simulate_shot(&input, &BallModel::TOUR, &given, &Surface::FAIRWAY);
            assert_eq!(a.flight.carry_m, b.flight.carry_m);
            assert_eq!(a.final_pos, b.final_pos);
        }
    }
}

#[test]
fn warmer_air_is_more_viscous_and_shortens_carry_at_constant_density() {
    let shot = ShotInput::seven_iron();
    let carry = |t| {
        simulate_shot(
            &shot,
            &BallModel::TOUR,
            &measured(1.2, t),
            &Surface::FAIRWAY,
        )
        .flight
        .carry_m
    };
    let (cold, mild, hot) = (carry(0.0), carry(20.0), carry(35.0));
    // A higher viscosity lowers the Reynolds number toward the drag crisis,
    // so drag rises: about 2.5 m from 0 °C to 35 °C for this 7-iron.
    assert!(cold > mild && mild > hot, "{cold} {mild} {hot}");
    assert!((1.0..5.0).contains(&(cold - hot)), "{}", cold - hot);
}

#[test]
fn bounce_hops_sample_the_same_profile() {
    let still = measured(1.2, 20.0);
    let shot = simulate_shot(
        &ShotInput::driver(),
        &BallModel::TOUR,
        &still,
        &Surface::FAIRWAY,
    );
    let p = shot.flight.points.last().unwrap();
    let landing = BallState {
        time: p.time,
        pos: p.pos,
        vel: p.vel,
        spin: shot.flight.landing_spin,
    };
    let at_10m = Vec3::new(0.0, 0.0, 5.0);
    let settle =
        |wind| Surface::FAIRWAY.settle(&landing, &BallModel::TOUR, &Environment { wind, ..still });
    let uniform = settle(Wind::Uniform(at_10m));
    let profile = settle(Wind::PowerLaw { at_10m });
    assert!(profile.contacts.len() > 1);
    assert!((uniform.rest.pos.z - profile.rest.pos.z).abs() > 0.001);
}
