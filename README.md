# dimple

[![CI](https://github.com/divotmaker/dimple/actions/workflows/ci.yml/badge.svg)](https://github.com/divotmaker/dimple/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/dimple.svg)](https://crates.io/crates/dimple)
[![docs.rs](https://docs.rs/dimple/badge.svg)](https://docs.rs/dimple)

## Disclaimer

dimple is an independent open-source project. It is not affiliated with,
endorsed by, or sponsored by TrackMan A/S, PGA TOUR, Inc., Acushnet Company,
Karsten Manufacturing Corporation, or the United States Golf Association.
TrackMan is a trademark of TrackMan A/S. PGA TOUR is a trademark of PGA TOUR,
Inc. Titleist is a trademark of Acushnet Company. PING is a trademark of
Karsten Manufacturing Corporation. USGA is a trademark of the United States
Golf Association. All other trademarks are the property of their respective
owners. These names appear only to identify the published sources the model
is compared against.

Validation figures are summary values quoted from publicly available
sources, cited in [docs/AERO.md](docs/AERO.md). That data belongs to its
publishers; dimple does not include or redistribute their datasets.

dimple produces estimates from a physics model. Its results are not
measurements and may differ from any launch monitor or from real ball
flight. The software is provided "as is", without warranty of any kind; see
the licenses below.

## Description

Golf ball-flight physics in pure Rust. Give it launch-monitor numbers (ball
speed, launch angle, azimuth, backspin, sidespin) and it returns the full
shot:

- **Flight**: fixed-step RK4 integration (1 ms) of gravity, Reynolds- and
  spin-dependent drag, Magnus lift, and exponential spin decay, with air
  density from altitude and temperature.
- **Bounce**: oblique impact on turf (Penner/Cross model) with restitution,
  friction, and turf compliance, including spin-back.
- **Rollout**: sliding then pure rolling with speed-dependent friction.
- **Stats**: per-shot carry, total, apex, descent angle, and lateral offset,
  plus session averages and dispersion.

The aerodynamic model is calibrated against published Trackman PGA Tour
averages. See [docs/AERO.md](docs/AERO.md) for the model, its sources, and
its validation.

`dimple` is pure math: no I/O, no async runtime, one dependency (`glam`).

## Quick Start

```toml
[dependencies]
dimple = "0.0"
```

```rust
use dimple::aero::BallModel;
use dimple::bounce::BounceSurface;
use dimple::environment::Environment;
use dimple::rollout::RolloutSurface;
use dimple::stats::ShotStats;
use dimple::trajectory::{ShotInput, simulate_flight, simulate_shot};

let input = ShotInput {
    ball_speed_mph: 167.0,
    launch_angle_deg: 10.5,
    launch_azimuth_deg: -1.0, // starts 1° left
    backspin_rpm: 2700.0,
    sidespin_rpm: -300.0,     // curves left: a draw for a right-hander
};

// Flight only: launch to first ground contact.
let flight = simulate_flight(&input, &BallModel::TOUR, &Environment::SEA_LEVEL);
println!("carry {:.1} yds, lateral {:+.1} m", flight.carry_yards, flight.lateral_m);

// Full shot: flight, bounces, and rollout.
let shot = simulate_shot(
    &input,
    &BallModel::TOUR,
    &Environment::DENVER,
    &BounceSurface::FAIRWAY,
    &RolloutSurface::FAIRWAY,
);
let stats = ShotStats::from_shot(&shot);
println!(
    "carry {:.1}, total {:.1} yds, apex {:.0} ft, lateral {:+.1} yds",
    stats.carry_yards, stats.total_yards, stats.apex_feet, stats.lateral_yards
);
```

`flight.points` and `shot.rollout.points` hold the sampled path for
animation; `trajectory::interpolate_trajectory` gives the position at any
time.

## Conventions

**Units.** Internally SI (m, s, kg, rad, rad/s). `ShotInput` takes
launch-monitor units (mph, degrees, rpm). `ShotStats` reports yards, feet,
and degrees. `dimple::units` has the conversions.

**Coordinate frame.** Right-handed, origin at the tee:

| Axis | Direction |
|---|---|
| +X | downrange, along the target line |
| +Y | up |
| +Z | right of the target line, from the golfer's view |

**Lateral signs: positive = right**, matching the
[Flight Relay Protocol](https://github.com/flightrelay/spec) and common
launch-monitor usage:

| Quantity | Positive means |
|---|---|
| `launch_azimuth_deg` | ball starts right of target |
| `sidespin_rpm` | ball curves right (fade/slice for a right-hander, draw/hook for a left-hander) |
| `lateral_m`, `lateral_yards` | ball lands right of target |
| `landing_sidespin_rads` | same as input sidespin |
| `z` of any position/velocity | right of / moving right of the target line |

**Spin axis.** Backspin and sidespin form one angular velocity vector
(right-hand rule): `ω = (0, -sidespin, backspin)`. Positive backspin points
along +Z and lifts the ball; positive sidespin (clockwise viewed from above)
points along -Y and curves it toward +Z.

## Ball Models

| Preset | Description |
|---|---|
| `BallModel::TOUR` | Modern (2020s) tour ball. Calibrated against Trackman 2024 PGA Tour averages. |
| `BallModel::PATENT_1997` | 1990s tour ball, fitted to published wind-tunnel data (US6186002B1). |

Every coefficient is a public field, so you can build your own `BallModel`.

## Validation

Trackman 2024 PGA Tour averages (25°C, sea level): carry RMSE 2.8 yds,
descent-angle RMSE 1.3°, 9 of 12 clubs within ±3 yds. Altitude and
temperature response track published data. Per-club results, sources, and
known limitations are in [docs/AERO.md](docs/AERO.md).

## License

Licensed under either of:

- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)

at your option.
