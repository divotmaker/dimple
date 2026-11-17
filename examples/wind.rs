//! Wind on a full shot: the same drive in still air, into a 5 m/s headwind and
//! with it, plus a 5 m/s crosswind. First as a uniform wind, then as a 10 m
//! reading scaled with height (`Wind::PowerLaw`).
//!
//! `cargo run --release --example wind`

use dimple::air::Wind;
use dimple::{ShotInput, Surface, Tour, Vec3, simulate_shot};

fn main() {
    let input = ShotInput {
        ball_speed_mph: 150.0,
        launch_angle_deg: 13.0,
        launch_azimuth_deg: 0.0,
        backspin_rpm: 2700.0,
        sidespin_rpm: 0.0,
    };

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
                "10 m wind scaled with height (Wind::PowerLaw)"
            } else {
                "uniform wind (Wind::Uniform)"
            }
        );
        println!(
            "{:<22} {:>9} {:>9} {:>9}",
            "wind", "carry m", "total m", "lateral m"
        );
        for (label, wind) in cases {
            let wind = if profile {
                Wind::PowerLaw { at_10m: wind }
            } else {
                Wind::Uniform(wind)
            };
            let shot = simulate_shot(&input, &Tour::new(1.225, wind), &Surface::FAIRWAY);
            println!(
                "{label:<22} {:>9.1} {:>9.1} {:>9.1}",
                shot.flight.carry_m, shot.total_m, shot.ground.rest.pos.z
            );
        }
    }
}
