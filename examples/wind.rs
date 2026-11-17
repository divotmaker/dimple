//! Wind on a full shot: the same drive in still air, into a 5 m/s headwind and
//! with it, plus a 5 m/s crosswind. First as a uniform wind, then as a 10 m
//! reading scaled with height (`WindProfile::PowerLaw`).
//!
//! `cargo run --release --example wind`

use dimple::Vec3;
use dimple::aero::BallModel;
use dimple::environment::{Atmosphere, Environment, WindProfile};
use dimple::ground::Surface;
use dimple::trajectory::{ShotInput, simulate_shot_with_atmosphere, simulate_shot_with_wind};

fn main() {
    let input = ShotInput {
        ball_speed_mph: 150.0,
        launch_angle_deg: 13.0,
        launch_azimuth_deg: 0.0,
        backspin_rpm: 2700.0,
        sidespin_rpm: 0.0,
    };
    let env = Environment::SEA_LEVEL;
    let atmosphere = Atmosphere::legacy(env.air_density());

    // +x is downrange, +z is right of the target line.
    let cases = [
        ("still", Vec3::ZERO),
        ("5 m/s headwind", Vec3::new(-5.0, 0.0, 0.0)),
        ("5 m/s tailwind", Vec3::new(5.0, 0.0, 0.0)),
        ("5 m/s left to right", Vec3::new(0.0, 0.0, 5.0)),
    ];

    for profile in [false, true] {
        println!(
            "\n{}",
            if profile {
                "10 m wind scaled with height (WindProfile::PowerLaw)"
            } else {
                "uniform wind (simulate_shot_with_wind)"
            }
        );
        println!(
            "{:<22} {:>9} {:>9} {:>9}",
            "wind", "carry m", "total m", "lateral m"
        );
        for (label, wind) in cases {
            let shot = if profile {
                simulate_shot_with_atmosphere(
                    &input,
                    &BallModel::TOUR,
                    &atmosphere,
                    &Surface::FAIRWAY,
                    WindProfile::PowerLaw { wind10: wind },
                )
            } else {
                simulate_shot_with_wind(&input, &BallModel::TOUR, &env, &Surface::FAIRWAY, wind)
            };
            println!(
                "{label:<22} {:>9.1} {:>9.1} {:>9.1}",
                shot.flight.carry_m, shot.total_m, shot.final_pos.z
            );
        }
    }
}
