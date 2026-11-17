# dimple

[![CI](https://github.com/divotmaker/dimple/actions/workflows/ci.yml/badge.svg)](https://github.com/divotmaker/dimple/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/dimple.svg)](https://crates.io/crates/dimple)
[![docs.rs](https://docs.rs/dimple/badge.svg)](https://docs.rs/dimple)

## Description

Golf ball-flight physics in pure Rust. Give it launch-monitor numbers (ball
speed, launch angle, azimuth, backspin, sidespin) and it returns the full
shot:

- **Flight**: fixed-step RK4 integration (1 ms) of gravity, spin-dependent
  drag, Magnus lift, and exponential spin decay. One `AeroModel` trait
  describes the ball and the air; `Tour`, a modern tour ball, is the
  built-in model. Your own ball, air or gravity (the Moon's, say) is a type
  that implements the trait.
- **Air**: density from altitude and temperature (`standard_air_density`)
  or a measured value, and wind: uniform, power-law with height, or a
  custom field of time, position, and height. **Beta:** wind is not
  yet compared with measured shots in wind, and the wind and measured
  density may change between releases.
- **Any starting state**: `simulate_flight_from` and `simulate_shot_from`
  fly on from any ball state, so a ball that hits a tree, a wall, or a cart
  path can bounce off it with its spin acting in the contact and continue
  from there.
- **Ground**: 3D rigid-ball impacts on turf (Penner's compliant-turf model
  with speed-dependent restitution, friction, and pivot friction), hops,
  then sliding and rolling with speed-dependent rolling resistance. The
  full spin vector is carried through every contact, so spin-back, the
  small kick of a curving ball off the turf, and drift from spin along the
  direction of travel all fall out of the physics. Impacts are checked
  against a published fit to bounce measurements, and rolling on slopes
  against published slope studies. **Beta:** run and total distances of
  full shots are compared with measured data only for drives on a fairway,
  and may change between releases.
- **Results**: carry, total, apex, descent angle, lateral offset, flight
  time, every sampled point, every impact, and the rest position, all in SI.

The aerodynamic model is calibrated against per-shot launch-monitor and
reference data and checked against published tour averages. See
[docs/AERO.md](docs/AERO.md) for the model, its sources, and its validation.

`dimple` is pure math: no I/O, no async runtime, no dependencies.

## Quick Start

```toml
[dependencies]
dimple = "0.0"
```

```rust
use dimple::units::{meters_to_feet, meters_to_yards};
use dimple::{ShotInput, Surface, Tour, Wind, simulate_flight, simulate_shot, standard_air_density};

let input = ShotInput {
    ball_speed_mph: 167.0,
    launch_angle_deg: 10.5,
    launch_azimuth_deg: -1.0, // starts 1° left
    backspin_rpm: 2700.0,
    sidespin_rpm: -300.0,     // curves left: a draw for a right-hander
};

// Flight only: launch to first ground contact, at sea level in still air.
let flight = simulate_flight(&input, &Tour::SEA_LEVEL);
println!(
    "carry {:.1} yds, lateral {:+.1} m",
    meters_to_yards(flight.carry_m),
    flight.lateral_m
);

// Full shot: flight, then bounces, slide and roll on a fairway, in Denver
// (1,609 m) at 17 °C with a 3 m/s headwind.
let air = Tour::new(standard_air_density(1609.0, 17.0), Wind::Uniform([-3.0, 0.0, 0.0].into()));
let shot = simulate_shot(&input, &air, &Surface::FAIRWAY);
println!(
    "carry {:.1}, total {:.1} yds, apex {:.0} ft, lateral {:+.1} yds",
    meters_to_yards(shot.flight.carry_m),
    meters_to_yards(shot.total_m),
    meters_to_feet(shot.flight.apex_m),
    meters_to_yards(shot.ground.rest.pos.z),
);
```

`flight.points` followed by `shot.ground.points` is the sampled path for
animation (ground times are absolute shot times);
`trajectory::interpolate_trajectory` gives the flight position at any time.

## The air

Every shot is flown by one `AeroModel`: the ball and the air it flies
through. Functions take it as `&dyn AeroModel`.

| You know | Model |
|---|---|
| nothing (sea level, standard day, still air) | `Tour::SEA_LEVEL` |
| altitude and temperature | `Tour::new(standard_air_density(altitude_m, temperature_c), wind)` |
| a measured density (weather station: ρ = p / (287.05 · T)) | `Tour::new(air_density, wind)` |
| a wind | `Wind::Uniform`, `Wind::PowerLaw` (a 10 m reading scaled with height), or `Wind::custom` (any function of time, position, and height, vertical component included) |
| another ball, air, or gravity | implement `AeroModel` yourself, or wrap `Tour` |

A wind is a horizontal air velocity in the frame below: +X blows downrange
(a tailwind), +Z blows to the right.

```rust
use dimple::{ShotInput, Tour, Vec3, Wind, simulate_flight};

// A gust that arrives two seconds into the flight: 8 m/s straight at the golfer.
let gust = Wind::custom(|t: f64, _pos: Vec3, _height: f64| {
    if t > 2.0 { Vec3::new(-8.0, 0.0, 0.0) } else { Vec3::ZERO }
});
let calm = simulate_flight(&ShotInput::driver(), &Tour::SEA_LEVEL);
let gusted = simulate_flight(&ShotInput::driver(), &Tour::new(1.225, gust));
assert!(gusted.carry_m < calm.carry_m);
```

A ball model of your own needs four methods (the air's acceleration on the
ball, the spin after a step, the radius, and the inertia factor); gravity is
a fifth, optional one, and dimple then uses it in flight and in every ground
contact. See [docs/AERO.md](docs/AERO.md#your-own-model) for a drag-only
ball and the details of the air, the wind, and measured density.

## From any ball state

`simulate_flight_from` and `simulate_shot_from` take a `BallState` (time,
position, velocity, spin) instead of a launch. Time continues from the
state's, carry and lateral offset are measured from its position, and a
state below the ground lands at once. With `ground::impact`, which bounces a
ball off a surface with any outward normal, a ball can hit a wall or a tree
and carry on:

```rust
use dimple::ground::{GroundSample, impact};
use dimple::ground::Restitution;
use dimple::{BallState, Surface, Tour, Vec3, simulate_shot_from};

// A hard, rigid wall (illustrative values) facing the golfer, at x = 40 m.
let wall = GroundSample {
    height: 0.0, // unused for an impact
    normal: Vec3::new(-1.0, 0.0, 0.0),
    surface: Surface {
        restitution: Restitution::Constant(0.5),
        compliance: 0.0,
        ..Surface::FAIRWAY
    },
};
let model = Tour::SEA_LEVEL;

// A ball arrives at the wall with sidespin about the vertical axis.
let arriving = BallState {
    time: 1.2,
    pos: Vec3::new(40.0, 6.0, 0.0),
    vel: Vec3::new(45.0, 8.0, 0.0),
    spin: Vec3::new(0.0, -250.0, 0.0),
};
let off_the_wall = impact(&arriving, &wall, &model);
assert!(off_the_wall.vel.x < 0.0); // back toward the golfer
assert!(off_the_wall.vel.z != 0.0); // the spin kicks it sideways

// Fly on from the rebound until it stops on the fairway.
let shot = simulate_shot_from(&off_the_wall, &model, &Surface::FAIRWAY);
assert!(shot.ground.rest.pos.x < 40.0);
```

## Ground

The ground API has three layers:

```rust
use dimple::ground::{FnTerrain, GroundPhysics, GroundSample};
use dimple::{ShotInput, Surface, Tour, Vec3, simulate_shot};

let (input, air) = (ShotInput::driver(), Tour::SEA_LEVEL);

// 1. One flat surface: a preset (FAIRWAY, GREEN, FIRM_GREEN), a green
//    set to a stimpmeter reading, or a preset with any field changed.
let shot = simulate_shot(&input, &air, &Surface::FAIRWAY);
let fast = Surface::GREEN.with_stimp(12.0);
assert!((fast.stimp_feet() - 12.0).abs() < 1e-9);
let shot = simulate_shot(&input, &air, &fast);
let firm = Surface { compliance: 0.75, ..Surface::FAIRWAY }.with_stimp(7.0);
let shot = simulate_shot(&input, &air, &firm);

// 2. dimple's physics over ground you describe point by point: here a
//    green from 260 m to 290 m downrange, fairway elsewhere. The surface is
//    sampled under the ball at every step, so a ball that lands on the
//    fairway and rolls onto the green feels the green.
let range = GroundPhysics::new(FnTerrain(|x: f64, _z: f64| GroundSample {
    height: 0.0,
    normal: Vec3::Y,
    surface: if (260.0..290.0).contains(&x) { Surface::GREEN } else { Surface::FAIRWAY },
}));
let shot = simulate_shot(&input, &air, &range);

// Slopes and height: a heightfield with its normal. Here the fairway
// rises 2° downrange; for height = a·x + b·z the normal is (−a, 1, −b)
// normalized. The ball lands where its flight first crosses the surface
// and slides and rolls with the slope's pull.
let a = 2.0_f64.to_radians().tan();
let uphill = GroundPhysics::new(FnTerrain(move |x: f64, _z: f64| GroundSample {
    height: a * x,
    normal: Vec3::new(-a, 1.0, 0.0).normalize(),
    surface: Surface::FAIRWAY,
}));
let shot = simulate_shot(&input, &air, &uphill);

// Tunables (step, hop transition, limits) are public fields too.
let quick = GroundPhysics { max_time: 10.0, ..GroundPhysics::new(Surface::GREEN) };
let shot = simulate_shot(&input, &air, &quick);
```

3. Implement `ground::GroundModel` to replace the physics entirely (it is
   object-safe, so `&dyn GroundModel` works); `ground::impact` and
   `GroundPhysics::step` are public building blocks, and `GroundResult`,
   `Contact` and `GroundPoint` have `new` constructors. A model for an
   obstacle such as a tree trunk finds where the flight meets it, calls
   `impact` with the trunk's outward normal, and continues with
   `simulate_shot_from`.

The terrain is a heightfield: one height, normal and surface per (x, z),
sampled under the ball's center. Flights launch from `(0, height(0, 0), 0)`
and land where they first cross the surface; carry is the horizontal
distance to that point. Each 1 ms flight or hop step is checked for a
crossing at 8 sub-steps, so a ridge or wall narrower than a sub-step (about
1 cm at full speed) can be missed. On a slope a ball at rest stays put where `tan θ ≤ mu_roll`
and rolls downhill otherwise. The heightfield has no walls, lips or cup, and on
curved ground the ball follows the surface (it never launches off a
crest). Walls, trunks and overhangs are not terrain: `ground::impact` takes
any normal, and a custom `GroundModel` places them. See [docs/AERO.md](docs/AERO.md#ground-model) for the contact
model and its limitations.

`GREEN` and `FIRM_GREEN` impact like Penner's (2002) green and firm green;
`FAIRWAY` bounces like `GREEN` and rolls slower (stimp 5.2 ft).
`Surface::with_stimp(feet)` sets the rolling resistance from a stimpmeter
reading in closed form, and `Surface::stimp_feet()` reads it back (`GREEN`
is 11.4 ft). There are no rough or bunker presets, because no published
set of this model's parameters exists for them; build such surfaces with
struct update syntax and `with_stimp` (see
[docs/AERO.md](docs/AERO.md#other-surfaces)).

## Conventions

**Units.** Internally SI (m, s, kg, rad, rad/s). `ShotInput` takes
launch-monitor units (mph, degrees, rpm). Results are SI, with the descent
angle in degrees. `dimple::units` converts to yards and feet.

**Vectors.** Every position, velocity, and spin in the API is a
`dimple::Vec3` in the frame below; convert it to your own math library's
vector through `[f64; 3]`.

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
| `lateral_m` | ball lands right of target |
| `z` of any position/velocity | right of / moving right of the target line |

**Spin axis.** Backspin and sidespin form one angular velocity vector
(right-hand rule). For launch azimuth `az`,
`ω = backspin · (−sin az, 0, cos az) − sidespin · ŷ`: backspin is about the
horizontal axis perpendicular to the start line (+Z for a straight shot) and
lifts the ball; positive sidespin (clockwise viewed from above) points along
-Y and curves it toward +Z. A straight launch gives
`ω = (0, -sidespin, backspin)`.

## The Tour model

`Tour` is a modern (2020s) tour ball, USGA conforming (45.93 g, 42.67 mm),
in a given air. Its drag and lift coefficients depend on the spin ratio
alone and are fitted to per-shot launch-monitor and reference data; they are
**beta** and may change between releases, so the constants are not part of
the API. A different ball is your own `AeroModel`.

## Validation

Per-shot launch-monitor and reference data (held out from the fit):
1,178 shots from 26 to 153 mph and 0.8k to 14k rpm come out at +0.1 yd mean
carry error and 1.6 yd RMSE, and the reference carries from sea level to
5,000 ft gain 5.7% on average against the model's 5.8%. Altitude and
temperature response track published data.

Trackman 2024 PGA Tour averages (25°C, sea level): carry bias +4.2 yd, RMSE
7.6 yd, with 5 of 12 clubs within ±3 yd. The model carries the woods 11 to
15 yd longer than that table and lands the irons 3° to 7° shallower. The
per-shot data agree with each other to about 1% and disagree with the tour
table the same way, so the fit follows the per-shot data. The preset
driver (167 mph, 10.5°, 2700 rpm) carries 282.7 yd at sea level and runs
to 300.1 yd on a fairway.

Ground (beta): on `GREEN` the impact model reproduces Penner's (2002)
published bounce heights, spin-back and drive runs (37.1 m and 17.9 m past
landing against his 37 m and 19 m), and sits at a scaled error of about 12%
from Biber et al.'s (2023) fit to natural-turf bounces. The PGA Tour
average drive runs 18 yd on `FAIRWAY`, between the USGA/R&A proposed bounce
model (17 yd at its 39° landing) and the ShotLink median (21 yd). Penner's
contact-angle law is faded in below 2 m/s of normal impact speed, which
ends the long train of millimeter hops the law otherwise sustains after
the first bounce.

Per-club results, sources, and known limitations are in
[docs/AERO.md](docs/AERO.md).

## License

Licensed under either of:

- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)

at your option.

## Disclaimer

dimple is an independent open-source project. It is not affiliated with,
endorsed by, or sponsored by TrackMan A/S, PGA TOUR, Inc., Acushnet Company,
Karsten Manufacturing Corporation, the United States Golf Association, or
The R&A. TrackMan is a trademark of TrackMan A/S. PGA TOUR and ShotLink are
trademarks of PGA TOUR, Inc. Titleist is a trademark of Acushnet Company. PING is a trademark of
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
