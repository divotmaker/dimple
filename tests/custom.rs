//! Models and winds of the caller's own, through the public API.

use dimple::air::Wind;
use dimple::ground::GroundModel;
use dimple::units::G;
use dimple::{
    AeroModel, BallState, ShotInput, Surface, Tour, Vec3, simulate_flight, simulate_shot,
};

/// A ball with drag only: a = −k |v| v.
struct DragOnly {
    k: f64,
}

impl AeroModel for DragOnly {
    fn acceleration(&self, ball: &BallState, _height: f64) -> Vec3 {
        ball.vel * (-self.k * ball.vel.length())
    }
    fn spin_after(&self, ball: &BallState, _dt: f64) -> Vec3 {
        ball.spin
    }
    fn radius_m(&self) -> f64 {
        0.021_335
    }
    fn inertia_factor(&self) -> f64 {
        0.4
    }
}

/// [`Tour`] under gravity `g`, or under the trait's default when `None`.
struct Gravity {
    tour: Tour,
    g: Option<f64>,
}

impl AeroModel for Gravity {
    fn acceleration(&self, ball: &BallState, height: f64) -> Vec3 {
        self.tour.acceleration(ball, height)
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
    fn gravity(&self) -> f64 {
        self.g.unwrap_or_else(|| self.tour.gravity())
    }
}

fn with_gravity(g: f64) -> Gravity {
    Gravity {
        tour: Tour::SEA_LEVEL,
        g: Some(g),
    }
}

#[test]
fn a_custom_model_flies_and_lands_through_dyn() {
    let model: Box<dyn AeroModel + Send + Sync> = Box::new(DragOnly { k: 0.004 });
    let shot = simulate_shot(&ShotInput::driver(), &*model, &Surface::FAIRWAY);
    assert!(shot.flight.carry_m > 50.0);
    assert!(shot.total_m.is_finite() && shot.total_m > 50.0);
}

#[test]
#[allow(clippy::float_cmp, reason = "the default gravity is the same constant")]
fn default_gravity_is_bit_identical_to_standard_gravity() {
    let input = ShotInput::seven_iron();
    let tour = simulate_shot(&input, &Tour::SEA_LEVEL, &Surface::FAIRWAY);
    for model in [
        Gravity {
            tour: Tour::SEA_LEVEL,
            g: None,
        },
        with_gravity(G),
    ] {
        let shot = simulate_shot(&input, &model, &Surface::FAIRWAY);
        assert_eq!(shot.flight.carry_m, tour.flight.carry_m);
        assert_eq!(shot.flight.apex_m, tour.flight.apex_m);
        assert_eq!(shot.ground, tour.ground);
    }
}

#[test]
fn lower_gravity_flies_and_rolls_farther() {
    let input = ShotInput::seven_iron();
    let earth = simulate_shot(&input, &Tour::SEA_LEVEL, &Surface::FAIRWAY);
    let moon = simulate_shot(&input, &with_gravity(1.62), &Surface::FAIRWAY);
    assert!(
        moon.flight.carry_m > earth.flight.carry_m,
        "{} vs {}",
        moon.flight.carry_m,
        earth.flight.carry_m
    );
    assert!(moon.flight.apex_m > earth.flight.apex_m);

    // The driver stays up for more than 30 s under 1.62 m/s²; it lands
    // rather than running out of flight time.
    let driver = simulate_shot(&ShotInput::driver(), &with_gravity(1.62), &Surface::FAIRWAY);
    let last = driver.flight.points.last().unwrap();
    assert!(driver.flight.flight_time > 30.0 && driver.flight.flight_time < 120.0);
    assert!(
        last.pos.y.abs() < 1e-9,
        "last flight point at y = {}",
        last.pos.y
    );

    // The same ball rolling away at 5 m/s on GREEN: rolling resistance
    // scales with gravity, so it rolls about 9.81 / 1.62 times as far.
    // (FAIRWAY's plowing is a drag that does not scale with gravity.)
    let rolling = BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel: Vec3::new(5.0, 0.0, 0.0),
        spin: Vec3::new(0.0, 0.0, -5.0 / Tour::SEA_LEVEL.radius_m()),
    };
    let on_earth = Surface::GREEN.settle(&rolling, &Tour::SEA_LEVEL).rest.pos.x;
    let on_moon = Surface::GREEN
        .settle(&rolling, &with_gravity(1.62))
        .rest
        .pos
        .x;
    assert!(
        on_moon > on_earth * 3.0,
        "roll {on_moon:.2} m vs {on_earth:.2} m"
    );
}

#[test]
#[allow(clippy::float_cmp, reason = "the same subtraction either way")]
fn a_constant_custom_wind_equals_the_uniform_wind() {
    let w = Vec3::new(-5.0, 0.0, 3.0);
    let uniform = Tour::new(1.2, Wind::Uniform(w));
    let custom = Tour::new(1.2, Wind::custom(move |_t: f64, _pos: Vec3, _h: f64| w));
    for input in [ShotInput::seven_iron(), ShotInput::driver()] {
        let a = simulate_shot(&input, &uniform, &Surface::FAIRWAY);
        let b = simulate_shot(&input, &custom, &Surface::FAIRWAY);
        assert_eq!(a.flight.carry_m, b.flight.carry_m);
        assert_eq!(a.flight.lateral_m, b.flight.lateral_m);
        assert_eq!(a.ground, b.ground);
    }
}

#[test]
fn a_gust_that_starts_mid_flight_needs_no_wrapper() {
    let head = Vec3::new(-8.0, 0.0, 0.0);
    let gust = Wind::custom(
        move |t: f64, _pos: Vec3, _h: f64| {
            if t > 2.0 { head } else { Vec3::ZERO }
        },
    );
    let input = ShotInput::driver();
    let calm = simulate_flight(&input, &Tour::SEA_LEVEL).carry_m;
    let gusted = simulate_flight(&input, &Tour::new(1.225, gust)).carry_m;
    let steady = simulate_flight(&input, &Tour::new(1.225, Wind::Uniform(head))).carry_m;
    assert!(
        steady < gusted && gusted < calm,
        "steady {steady}, gusted {gusted}, calm {calm}"
    );
}

#[test]
fn an_updraft_raises_the_apex_and_a_downdraft_lowers_it() {
    let draft = |vy: f64| {
        Tour::new(
            1.225,
            Wind::custom(move |_t: f64, _p: Vec3, _h: f64| Vec3::new(0.0, vy, 0.0)),
        )
    };
    for input in [ShotInput::driver(), ShotInput::seven_iron()] {
        let still = simulate_flight(&input, &Tour::SEA_LEVEL);
        let up = simulate_flight(&input, &draft(3.0));
        let down = simulate_flight(&input, &draft(-3.0));
        assert!(
            down.apex_m < still.apex_m && still.apex_m < up.apex_m,
            "apex {} / {} / {}",
            down.apex_m,
            still.apex_m,
            up.apex_m
        );
        assert!(up.flight_time > still.flight_time);
    }
}
