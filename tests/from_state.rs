//! Flights and shots that start from any ball state.

use dimple::air::Wind;
use dimple::ground::{FnTerrain, GroundPhysics, GroundSample};
use dimple::{
    AeroModel, BallState, ShotInput, Surface, Tour, Vec3, simulate_flight, simulate_flight_from,
    simulate_shot, simulate_shot_from, units,
};

/// The launch state of `input` at `(0, tee_height, 0)`, built from the
/// public API: the launch velocity along the launch angle and azimuth, and
/// [`ShotInput::launch_spin_rads`].
fn launch(input: &ShotInput, tee_height: f64) -> BallState {
    let speed = units::mph_to_mps(input.ball_speed_mph);
    let launch = units::deg_to_rad(input.launch_angle_deg);
    let azimuth = units::deg_to_rad(input.launch_azimuth_deg);
    BallState {
        time: 0.0,
        pos: Vec3::new(0.0, tee_height, 0.0),
        vel: Vec3::new(
            speed * launch.cos() * azimuth.cos(),
            speed * launch.sin(),
            speed * launch.cos() * azimuth.sin(),
        ),
        spin: input.launch_spin_rads(),
    }
}

fn inputs() -> [ShotInput; 4] {
    [
        ShotInput::driver(),
        ShotInput::seven_iron(),
        ShotInput::wedge(),
        ShotInput {
            ball_speed_mph: 120.0,
            launch_angle_deg: 16.3,
            launch_azimuth_deg: 1.5,
            backspin_rpm: 7097.0,
            sidespin_rpm: 900.0,
        },
    ]
}

/// A sloped green: down 3% downrange, up 2% to the right.
fn slope() -> GroundPhysics<FnTerrain<impl Fn(f64, f64) -> GroundSample>> {
    GroundPhysics::new(FnTerrain(|x: f64, z: f64| GroundSample {
        height: 0.5 - 0.03 * x + 0.02 * z,
        normal: Vec3::new(0.03, 1.0, -0.02).normalize(),
        surface: Surface::GREEN,
    }))
}

#[test]
#[allow(clippy::float_cmp, reason = "the wrappers must be bit-identical")]
fn the_wrappers_equal_the_from_calls_bit_for_bit() {
    let wind = Tour::new(1.2, Wind::Uniform(Vec3::new(-4.0, 0.0, 2.0)));
    for model in [&Tour::SEA_LEVEL, &wind] {
        for input in inputs() {
            let a = simulate_flight(&input, model);
            let b = simulate_flight_from(&launch(&input, 0.0), model);
            assert_eq!(format!("{a:?}"), format!("{b:?}"));
            assert_eq!(a.carry_m, b.carry_m);

            let a = simulate_shot(&input, model, &Surface::FAIRWAY);
            let b = simulate_shot_from(&launch(&input, 0.0), model, &Surface::FAIRWAY);
            assert_eq!(format!("{a:?}"), format!("{b:?}"));
            assert_eq!(a.ground, b.ground);

            let ground = slope();
            let a = simulate_shot(&input, model, &ground);
            let b = simulate_shot_from(&launch(&input, 0.5), model, &ground);
            assert_eq!(format!("{a:?}"), format!("{b:?}"));
            assert_eq!(a.total_m, b.total_m);
        }
    }
}

/// The state of `flight` at its sampled point `index` (taken every 10 flight
/// steps, before the landing), with the spin the model gives after as many
/// steps.
fn state_at(
    flight: &dimple::trajectory::FlightResult,
    index: usize,
    input: &ShotInput,
) -> BallState {
    let p = flight.points[index];
    let mut ball = launch(input, 0.0);
    for _ in 0..index * 10 {
        ball.spin = Tour::SEA_LEVEL.spin_after(&ball, 0.001);
    }
    BallState {
        time: p.time,
        pos: p.pos,
        vel: p.vel,
        spin: ball.spin,
    }
}

#[test]
#[allow(
    clippy::float_cmp,
    reason = "the resumed flight starts at the given time"
)]
fn a_flight_split_and_resumed_lands_where_the_whole_flight_does() {
    // Still air, and a crosswind that strengthens with time: the resumed
    // flight must carry on with the shot's own clock.
    let gusting = Tour::new(
        1.225,
        Wind::custom(|t: f64, _pos: Vec3, _h: f64| Vec3::new(-2.0 * t, 0.0, 1.5 * t)),
    );
    for model in [Tour::SEA_LEVEL, gusting] {
        for input in inputs() {
            let whole = simulate_flight(&input, &model);
            let land = whole.points.last().unwrap();
            for index in [1, 50, 100, whole.points.len() - 3] {
                let start = state_at(&whole, index, &input);
                let rest = simulate_flight_from(&start, &model);
                let end = rest.points.last().unwrap();
                assert!(
                    (end.pos - land.pos).length() < 1e-6,
                    "split at {}: {:?} vs {:?}",
                    start.time,
                    end.pos,
                    land.pos
                );
                assert!((end.time - land.time).abs() < 1e-9);
                assert_eq!(rest.points[0].time, start.time);
                assert!((rest.flight_time - (land.time - start.time)).abs() < 1e-9);
                let from_start = (end.pos.x - start.pos.x).hypot(end.pos.z - start.pos.z);
                assert!((rest.carry_m - from_start).abs() < 1e-12);
                assert!((rest.lateral_m - (end.pos.z - start.pos.z)).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn a_shot_split_mid_flight_comes_to_rest_where_the_whole_shot_does() {
    let input = ShotInput::seven_iron();
    let whole = simulate_shot(&input, &Tour::SEA_LEVEL, &Surface::FAIRWAY);
    let start = state_at(&whole.flight, 200, &input);
    let rest = simulate_shot_from(&start, &Tour::SEA_LEVEL, &Surface::FAIRWAY);
    assert!((rest.ground.rest.pos - whole.ground.rest.pos).length() < 1e-6);
    assert_eq!(rest.ground.contacts.len(), whole.ground.contacts.len());
    let to_rest = (whole.ground.rest.pos.x - start.pos.x).hypot(whole.ground.rest.pos.z);
    assert!((rest.total_m - to_rest).abs() < 1e-6);
}

#[test]
#[allow(clippy::float_cmp, reason = "an immediate landing is exact")]
fn a_start_below_the_ground_lands_at_once_on_it() {
    let start = BallState {
        time: 2.0,
        pos: Vec3::new(10.0, -0.3, 1.0),
        vel: Vec3::new(5.0, 3.0, 0.0),
        spin: Vec3::new(0.0, 0.0, 100.0),
    };
    let flight = simulate_flight_from(&start, &Tour::SEA_LEVEL);
    assert_eq!(flight.points.len(), 1);
    let p = flight.points[0];
    assert_eq!(
        (p.time, p.pos, p.vel),
        (2.0, Vec3::new(10.0, 0.0, 1.0), start.vel)
    );
    assert_eq!(flight.flight_time, 0.0);
    assert_eq!(flight.carry_m, 0.0);
    assert_eq!(flight.apex_m, 0.0);
    assert_eq!(flight.landing_spin, start.spin);

    // A full shot from under the ground, flat and on terrain: no panic, a
    // finite rest on the surface.
    let ground = slope();
    let flat = simulate_shot_from(&start, &Tour::SEA_LEVEL, &Surface::FAIRWAY);
    let sloped = simulate_shot_from(&start, &Tour::SEA_LEVEL, &ground);
    for shot in [flat, sloped] {
        let first = shot.ground.contacts[0];
        assert_eq!((first.time, first.pos.x, first.pos.z), (2.0, 10.0, 1.0));
        assert!(shot.total_m.is_finite());
        assert!(shot.ground.rest.pos.length().is_finite());
    }
    let first = simulate_shot_from(&start, &Tour::SEA_LEVEL, &ground)
        .ground
        .contacts[0];
    assert_eq!(first.pos.y, 0.5 - 0.03 * 10.0 + 0.02 * 1.0);
}

#[test]
fn a_start_on_the_ground_moving_into_it_lands_within_the_first_step() {
    let start = BallState {
        time: 1.0,
        pos: Vec3::new(3.0, 0.0, 0.0),
        vel: Vec3::new(10.0, -4.0, 0.0),
        spin: Vec3::ZERO,
    };
    let flight = simulate_flight_from(&start, &Tour::SEA_LEVEL);
    let end = flight.points.last().unwrap();
    assert!(flight.flight_time < 1e-9, "{}", flight.flight_time);
    assert!((end.pos - start.pos).length() < 1e-9);
}
