use dimple::aero::BallModel;
use dimple::environment::{Atmosphere, Environment, WindProfile};
use dimple::ground::{BallState, GroundModel, Surface};
use dimple::trajectory::{ShotInput, simulate_shot_with_atmosphere, simulate_shot_with_wind};
use glam::DVec3;

#[test]
fn legacy_entry_point_keeps_identical_results() {
    let env = Environment::SEA_LEVEL;
    let atmosphere = Atmosphere::legacy(env.air_density());
    for wind in [DVec3::ZERO, DVec3::new(-5.0, 0.0, 3.0)] {
        let input = ShotInput::driver();
        let old = simulate_shot_with_wind(&input, &BallModel::TOUR, &env, &Surface::FAIRWAY, wind);
        let new = simulate_shot_with_atmosphere(
            &input,
            &BallModel::TOUR,
            &atmosphere,
            &Surface::FAIRWAY,
            WindProfile::Uniform(wind),
        );
        assert_eq!(old.flight.carry_m, new.flight.carry_m);
        assert_eq!(old.final_pos, new.final_pos);
    }
}

#[test]
fn warmer_air_is_more_viscous_and_shortens_carry_at_constant_density() {
    let shot = ShotInput::seven_iron();
    let carry = |t| {
        simulate_shot_with_atmosphere(
            &shot,
            &BallModel::TOUR,
            &Atmosphere::from_density_and_temperature(1.2, t).unwrap(),
            &Surface::FAIRWAY,
            WindProfile::Uniform(DVec3::ZERO),
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
    let atmosphere = Atmosphere::from_density_and_temperature(1.2, 20.0).unwrap();
    let shot = simulate_shot_with_atmosphere(
        &ShotInput::driver(),
        &BallModel::TOUR,
        &atmosphere,
        &Surface::FAIRWAY,
        WindProfile::Uniform(DVec3::ZERO),
    );
    let p = shot.flight.points.last().unwrap();
    let landing = BallState {
        time: p.time,
        pos: p.pos,
        vel: p.vel,
        spin: shot.flight.landing_spin,
    };
    let wind10 = DVec3::new(0.0, 0.0, 5.0);
    let uniform = Surface::FAIRWAY.settle_with_atmosphere(
        &landing,
        &BallModel::TOUR,
        &atmosphere,
        WindProfile::Uniform(wind10),
    );
    let profile = Surface::FAIRWAY.settle_with_atmosphere(
        &landing,
        &BallModel::TOUR,
        &atmosphere,
        WindProfile::PowerLaw { wind10 },
    );
    assert!(profile.contacts.len() > 1);
    assert!((uniform.rest.pos.z - profile.rest.pos.z).abs() > 0.001);
}
