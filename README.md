# dimple

[![CI](https://github.com/divotmaker/dimple/actions/workflows/ci.yml/badge.svg)](https://github.com/divotmaker/dimple/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/dimple.svg)](https://crates.io/crates/dimple)
[![docs.rs](https://docs.rs/dimple/badge.svg)](https://docs.rs/dimple)

## Description

Golf ball-flight physics in pure Rust. Give it launch-monitor numbers (ball
speed, launch angle, azimuth, backspin, sidespin) and it returns the full
shot:

- **Flight**: fixed-step RK4 integration (1 ms) of gravity, Reynolds- and
  spin-dependent drag, Magnus lift, and exponential spin decay, with air
  density from altitude and temperature.
- **Ground**: 3D rigid-ball impacts on turf (Penner's compliant-turf model
  with speed-dependent restitution, friction, and pivot friction), hops,
  then sliding and rolling with speed-dependent rolling resistance. The
  full spin vector is carried through every contact, so spin-back, the
  small kick of a curving ball off the turf, and drift from spin along the
  direction of travel all fall out of the physics.
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
use dimple::environment::Environment;
use dimple::ground::Surface;
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

// Full shot: flight, then bounces, slide and roll on a fairway.
let shot = simulate_shot(&input, &BallModel::TOUR, &Environment::DENVER, &Surface::FAIRWAY);
let stats = ShotStats::from_shot(&shot);
println!(
    "carry {:.1}, total {:.1} yds, apex {:.0} ft, lateral {:+.1} yds",
    stats.carry_yards, stats.total_yards, stats.apex_feet, stats.lateral_yards
);
```

`flight.points` followed by `shot.ground.points` is the sampled path for
animation (ground times are absolute shot times);
`trajectory::interpolate_trajectory` gives the flight position at any time.

## Ground

The ground API has three layers:

```rust
use dimple::aero::BallModel;
use dimple::environment::Environment;
use dimple::ground::{FnTerrain, GroundPhysics, GroundSample, Surface};
use dimple::trajectory::{ShotInput, simulate_shot};
use glam::DVec3;

let (input, ball, env) = (ShotInput::driver(), BallModel::TOUR, Environment::SEA_LEVEL);

// 1. One flat surface: a preset (FAIRWAY, GREEN, FIRM_GREEN), a green
//    set to a stimpmeter reading, or a preset with any field changed.
let shot = simulate_shot(&input, &ball, &env, &Surface::FAIRWAY);
let fast = Surface::GREEN.with_stimp(12.0);
assert!((fast.stimp_feet() - 12.0).abs() < 1e-9);
let shot = simulate_shot(&input, &ball, &env, &fast);
let firm = Surface { compliance: 0.75, ..Surface::FAIRWAY }.with_stimp(7.0);
let shot = simulate_shot(&input, &ball, &env, &firm);

// 2. dimple's physics over ground you describe point by point: here a
//    green from 260 m to 290 m downrange, fairway elsewhere. The surface is
//    sampled under the ball at every step, so a ball that lands on the
//    fairway and rolls onto the green feels the green.
let range = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| GroundSample {
    height: 0.0,
    normal: DVec3::Y,
    surface: if (260.0..290.0).contains(&x) { Surface::GREEN } else { Surface::FAIRWAY },
}));
let shot = simulate_shot(&input, &ball, &env, &range);

// Slopes and height: a heightfield with its normal. Here the fairway
// rises 2° downrange; for height = a·x + b·z the normal is (−a, 1, −b)
// normalised. The ball lands where its flight first crosses the surface
// and slides and rolls with the slope's pull.
let a = 2.0_f64.to_radians().tan();
let uphill = GroundPhysics::new(FnTerrain(move |x: f64, _z: f64| GroundSample {
    height: a * x,
    normal: DVec3::new(-a, 1.0, 0.0).normalize(),
    surface: Surface::FAIRWAY,
}));
let shot = simulate_shot(&input, &ball, &env, &uphill);

// Tunables (step, hop transition, limits) are public fields too.
let quick = GroundPhysics { max_time: 10.0, ..GroundPhysics::new(Surface::GREEN) };
let shot = simulate_shot(&input, &ball, &env, &quick);
```

3. Implement `ground::GroundModel` to replace the physics entirely (it is
   object-safe, so `&dyn GroundModel` works); `ground::impact` and
   `GroundPhysics::step` are public building blocks, and `GroundResult`,
   `Contact` and `GroundPoint` have `new` constructors.

The terrain is a heightfield: one height, normal and surface per (x, z),
sampled under the ball's centre. Flights launch from `(0, height(0, 0), 0)`
and land where they first cross the surface; carry is the horizontal
distance to that point. Each 1 ms flight or hop step is checked for a
crossing at 8 sub-steps, so a ridge or wall narrower than a sub-step (about
1 cm at full speed) can be missed. On a slope a ball at rest stays put where `tan θ ≤ mu_roll`
and rolls downhill otherwise. There are no walls, lips or cup, and on
curved ground the ball follows the surface (it never launches off a
crest). See [docs/AERO.md](docs/AERO.md#ground-model) for the contact
model and its limitations.

`GREEN` and `FIRM_GREEN` impact like Penner's (2002) green and firm green;
`FAIRWAY` uses the same impact model with a chosen compliance.
`Surface::with_stimp(feet)` sets the rolling resistance from a stimpmeter
reading in closed form, and `Surface::stimp_feet()` reads it back (`GREEN`
is 11.4 ft). There are no rough or bunker presets, because no published
set of this model's parameters exists for them; build such surfaces with
struct update syntax and `with_stimp` (see
[docs/AERO.md](docs/AERO.md#other-surfaces)).

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
(right-hand rule). For launch azimuth `az`,
`ω = backspin · (−sin az, 0, cos az) − sidespin · ŷ`: backspin is about the
horizontal axis perpendicular to the start line (+Z for a straight shot) and
lifts the ball; positive sidespin (clockwise viewed from above) points along
-Y and curves it toward +Z. A straight launch gives
`ω = (0, -sidespin, backspin)`.

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
the licenses above.
