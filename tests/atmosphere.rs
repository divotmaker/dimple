//! Air density and wind profile through the public API.

use dimple::air::Wind;
use dimple::ground::GroundModel;
use dimple::{BallState, ShotInput, Surface, Tour, Vec3, simulate_shot, standard_air_density};

#[test]
#[allow(
    clippy::float_cmp,
    reason = "standard air at sea level is exactly the preset's density"
)]
fn sea_level_preset_is_standard_air_at_sea_level() {
    assert_eq!(standard_air_density(0.0, 15.0), 1.225);
    for wind in [Wind::None, Wind::Uniform(Vec3::new(-5.0, 0.0, 3.0))] {
        let preset = Tour::new(1.225, wind.clone());
        let computed = Tour::new(standard_air_density(0.0, 15.0), wind);
        let input = ShotInput::driver();
        let a = simulate_shot(&input, &preset, &Surface::FAIRWAY);
        let b = simulate_shot(&input, &computed, &Surface::FAIRWAY);
        assert_eq!(a.flight.carry_m, b.flight.carry_m);
        assert_eq!(a.ground.rest.pos, b.ground.rest.pos);
    }
    let still = simulate_shot(&ShotInput::driver(), &Tour::SEA_LEVEL, &Surface::FAIRWAY);
    let built = simulate_shot(
        &ShotInput::driver(),
        &Tour::new(1.225, Wind::None),
        &Surface::FAIRWAY,
    );
    assert_eq!(still.ground.rest.pos, built.ground.rest.pos);
}

#[test]
fn altitude_gain_for_the_preset_driver() {
    // 0 → 1,524 m (5,000 ft) on a standard day: ISA density 1.0555 kg/m³.
    let input = ShotInput::driver();
    let sea = simulate_shot(&input, &Tour::SEA_LEVEL, &Surface::FAIRWAY)
        .flight
        .carry_m;
    let high_air = Tour::new(
        standard_air_density(1524.0, 15.0 - 0.0065 * 1524.0),
        Wind::None,
    );
    let high = simulate_shot(&input, &high_air, &Surface::FAIRWAY)
        .flight
        .carry_m;
    let gain = high / sea - 1.0;
    assert!((0.04..0.08).contains(&gain), "gain {:.2} %", gain * 100.0);
}

#[test]
fn bounce_hops_sample_the_same_profile() {
    let still = Tour::new(1.2, Wind::None);
    let shot = simulate_shot(&ShotInput::driver(), &still, &Surface::FAIRWAY);
    let p = shot.flight.points.last().unwrap();
    let landing = BallState {
        time: p.time,
        pos: p.pos,
        vel: p.vel,
        spin: shot.flight.landing_spin,
    };
    let at_10m = Vec3::new(0.0, 0.0, 5.0);
    let settle = |wind| Surface::FAIRWAY.settle(&landing, &Tour::new(1.2, wind));
    let uniform = settle(Wind::Uniform(at_10m));
    let profile = settle(Wind::PowerLaw { at_10m });
    assert!(profile.contacts.len() > 1);
    assert!((uniform.rest.pos.z - profile.rest.pos.z).abs() > 0.001);
}
