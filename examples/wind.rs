//! Wind on a full shot: the same drive in still air, into a 5 m/s headwind and
//! with it, plus a 5 m/s crosswind.
//!
//! `cargo run --release --example wind`

use dimple::aero::BallModel;
use dimple::environment::Environment;
use dimple::ground::Surface;
use dimple::trajectory::{ShotInput, simulate_shot_with_wind};
use glam::DVec3;

fn main() {
    let input = ShotInput {
        ball_speed_mph: 150.0,
        launch_angle_deg: 13.0,
        launch_azimuth_deg: 0.0,
        backspin_rpm: 2700.0,
        sidespin_rpm: 0.0,
    };
    let env = Environment::SEA_LEVEL;

    // +x is downrange, +z is right of the target line.
    let cases = [
        ("still", DVec3::ZERO),
        ("5 m/s headwind", DVec3::new(-5.0, 0.0, 0.0)),
        ("5 m/s tailwind", DVec3::new(5.0, 0.0, 0.0)),
        ("5 m/s left to right", DVec3::new(0.0, 0.0, 5.0)),
    ];

    println!(
        "{:<22} {:>9} {:>9} {:>9}",
        "wind", "carry m", "total m", "lateral m"
    );
    for (label, wind) in cases {
        let shot = simulate_shot_with_wind(&input, &BallModel::TOUR, &env, &Surface::FAIRWAY, wind);
        println!(
            "{label:<22} {:>9.1} {:>9.1} {:>9.1}",
            shot.flight.carry_m, shot.total_m, shot.final_pos.z
        );
    }
}
