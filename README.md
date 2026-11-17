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
- **Ground**: 3D rigid-ball impacts on turf (a tilted contact plane, after
  Penner, with speed-dependent restitution, friction, and pivot friction),
  hops, then sliding and rolling with a Stimpmeter-set rolling resistance
  that rises with speed. The full spin vector is carried through every
  contact, so spin-back, the small kick of a curving ball off the turf, and
  drift from spin along the direction of travel all fall out of the
  physics. The contact angle is fitted to 693 measured bounces on natural
  tee turf at 2–39 m/s and landing angles of 16–90°, and the presets'
  compliance and friction are the medians of measured greens and fairways.
  Roll-out is within 3 yd of the USGA's 2019 range test at 20–36° (rms
  2 yd) and about 11 yd long at 12–18°; landings flatter than 16° are not
  measured.
- **Results**: carry, total, apex, descent angle, lateral offset, flight
  time, every sampled point, every impact, and the rest position, all in SI.

The aerodynamic model is calibrated against per-shot launch-monitor and
reference data and checked against published tour averages. See
[docs/AERO.md](docs/AERO.md) for the model, its sources, and its validation,
and [docs/FIT.md](docs/FIT.md) for how the ground model scores against each
dataset.

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
let firm = Surface { compliance: 0.9, ..Surface::FAIRWAY }.with_stimp(7.0);
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

`FAIRWAY`, `GREEN` and `FIRM_GREEN` carry measured constants: `FAIRWAY`
the medians of eight fairways (stimp 5 ft 8 in), `GREEN` the medians of
nine greens (stimp 10.5 ft), `FIRM_GREEN` one green dried to championship
firmness (stimp 10.5 ft). `compliance` is firmness: on greens it follows
the penetration depth of the USGA firmness meter,
`compliance ≈ 0.090 × depth (mm) + 0.25` over 8.8–14.3 mm.
`Surface::with_stimp(feet)` sets the rolling resistance from a stimpmeter
reading, and `Surface::stimp_feet()` reads it back. There are no rough or
bunker presets; see [Rough and sand](#rough-and-sand).

## Rough and sand

No published measurements exist that could set this model's constants for
rough or sand. No published test fires balls into rough, and none measures
a ball slowing through grass; the one published figure is that PGA TOUR
drives finishing in rough lose 8–9 yd of bounce and roll, with a weaker
dependence on landing angle (−0.45 yd per degree against −1.24; Quintavalla
& Rogers 2008, in USGA/R&A 2021). For sand, a ball dropped on loose sand
comes to rest part-buried (Uehara et al. 2003; Baker 2001), and the rolling
resistance of a ball on a granular bed is measured (Fukumoto et al. 2026;
De Blasio & Saeter 2009), but no measurement of a golf ball's speed or spin
leaving bunker sand exists.

So build them from a preset and tune by feel. The fields that matter:

- `plow` and `mu_roll` stop the ball sooner once it is on the ground:
  `plow` most at speed, `mu_roll` most as it slows (and on slopes, where
  it sets how steep a lie holds the ball).
- `absorption` makes a landing die: the share of the rebound the ground
  takes from a ball landing straight down, less for a shallow one. High
  values stop or "plug" the ball by its mark, as in sand.
- `compliance` softens the bounce: higher tilts the contact plane more, so
  the ball leaves steeper and slower and runs less.
- `mu_slide` sets how much the ball checks: higher friction takes more
  speed and spin out of each impact.

Starting values and sensible ranges. **These are chosen, not measured.**

| Field | Rough: start | Rough: range | Sand: start | Sand: range |
|---|---|---|---|---|
| `restitution` | `Penner2002` | `Penner2002` | `Constant(0.0)` | `Constant(0.0)`–`Constant(0.1)` |
| `compliance` | 1.15 | 1.05–1.5 | 0.0 | 0.0–0.3 |
| `mu_slide` | 0.18 | 0.15–0.4 | 0.40 | 0.2–0.6 |
| `mu_roll` | 0.3 | 0.15–0.5 | 0.3 | 0.25–0.5 |
| `plow` (1/s) | 2.0 | 0.5–4 | 0.36 | 0.36–2 |
| `absorption` | 0.0 | 0–0.3 | 0.9 | 0.8–0.95 |

```rust
use dimple::ground::{Restitution, Surface};

// Rough: a soft fairway's bounce that stops quickly once on the ground.
let rough = Surface {
    compliance: 1.15,
    mu_roll: 0.3,
    plow: 2.0,
    ..Surface::FAIRWAY
};
// Bunker sand: no bounce; the ball skips on, and the steeper it lands
// the less it keeps.
let sand = Surface {
    restitution: Restitution::Constant(0.0),
    compliance: 0.0,
    absorption: 0.9,
    mu_slide: 0.40,
    mu_roll: 0.3,
    ..Surface::FAIRWAY
};
```

The rough starts from `FAIRWAY` with `compliance` 1.15, inside the 0.86–1.55
measured across single fairways, so it bounces like a soft fairway; its
`mu_roll` and `plow` make it stop. With these values the preset driver runs
12.0 yd against `FAIRWAY`'s 21.1 yd, about the 8–9 yd that PGA TOUR drives
finishing in rough lose, and its run falls 0.76 yd per degree of landing
angle against `FAIRWAY`'s 1.03. The sand is the sand example of
[docs/AERO.md](docs/AERO.md#other-surfaces), where `absorption` and its
measured bounds are explained; what it does to wedge, iron and driver
landings is tabled there.

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
to 303.8 yd on a fairway.

Ground: the contact angle is fitted to 693 bounces on natural tee turf at
2–39 m/s and landing angles of 16–90° (Biber et al. 2023). `GREEN`
reproduces Cross's (2025) printed fits to 184 bounces on four greens at
17–27 m/s and 35–55°; `FIRM_GREEN` is fitted to 57 landings at 20–26 m/s
and 36–59° on one firm green. Roll-out is within 3 yd of the USGA's 2019
range test at 20–36° (rms 2 yd) and about 11 yd long at 12–18°, below the
16° to which impacts are measured. The PGA Tour average drive runs
23 yd on `FAIRWAY` against 21.7 yd on the PGA TOUR's mean roll-out line at
its 39° landing, and the R&A's driver landing bounces and rolls 24.8 yd
against 4–31 yd measured on five fairways.

Per-club results, sources, and known limitations are in
[docs/AERO.md](docs/AERO.md).

## References

Each entry gives what dimple takes from the source. Details of each model and
check are in [docs/AERO.md](docs/AERO.md).

### Ground: impact model

- Penner, A. R. (2002), "The run of a golf ball", *Canadian Journal of Physics*
  80(8), 931–940. [doi:10.1139/p02-035](https://doi.org/10.1139/p02-035)
  - Impact against a tilted effective plane, and the restitution against normal
    speed (eq. 5) that every preset uses.
- Penner, A. R. (2001), "The physics of golf: The optimum loft of a driver",
  *American Journal of Physics* 69(5), 563–568.
  [doi:10.1119/1.1344164](https://doi.org/10.1119/1.1344164)
  - The turf rebound measurements behind eq. 5.
- Daish, C. B. (1972), *The Physics of Ball Games*, English Universities Press,
  London.
  - The rigid-sphere grip-or-slide rebound, applied in the tilted frame.
- Biber, S. W., Jones, K. M., Champneys, A. R., Green, R. and Szalai, R. (2023),
  "Measurements and linearized models for golf ball bounce on a green",
  *Sports Engineering* 26, 50.
  [doi:10.1007/s12283-023-00442-4](https://doi.org/10.1007/s12283-023-00442-4)
  (open access, CC BY 4.0)
  - The contact-angle (tilt) law, fitted to its 693 bounces on natural tee turf
    at 2–39 m/s and 16–90°; a check on eq. 5; the mat example.
- Johnson, K. L. (1985), *Contact Mechanics*, Cambridge University Press.
  [doi:10.1017/CBO9781139171731](https://doi.org/10.1017/CBO9781139171731)
  - The torsional (pivot) friction lever of a Hertzian contact.

### Ground: measured bounces on turf

- Cross, R. (2025), "Bounce of a golf ball on the green", *Proceedings of the
  Institution of Mechanical Engineers, Part P: Journal of Sports Engineering
  and Technology* 239(3), 479–492.
  [doi:10.1177/17543371221150289](https://doi.org/10.1177/17543371221150289)
  - `GREEN` compliance and friction checked against its printed fits (eqs.
    14–16) over 184 bounces on four greens at 17–27 m/s and 35–55°; its
    low-speed chips (Table 1) mark where that fit ends; its rolling friction
    of 0.04–0.09 brackets the greens' `mu_roll`.
- Cross, R. (2021), "Oblique impact of a hard ball on a soft surface",
  *European Journal of Physics* 42, 065005.
  [doi:10.1088/1361-6404/ac221e](https://doi.org/10.1088/1361-6404/ac221e)
  - Bounce on mown, fairway-like grass at 31–85°: a check on `FAIRWAY`; its
    vertical drops check eq. 5.
- Johnson, A. D. G., Otto, S. R., Richardson, P. C. W. and Isaac, S. P.
  (2008), *An Investigation into Ball-Turf Impact Characteristics*, The R&A,
  Distance Insights library report R04.
  [PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/2887f275-0659-4812-a1be-30ab00a8472c/R04%20-%20An%20Investigation%20into%20Ball-Turf%20Impact%20Characteristics%20_Master%20Copy_%20_5.pdf)
  - `GREEN` and `FAIRWAY` compliance and sliding friction from five fairways
    and six greens (driver, 5-iron and sand wedge landings); `FAIRWAY`
    plowing from its roll-outs; `compliance` as firmness.
- USGA Research and Test Center (2005), *Preliminary Results of Turf Testing
  at Pinehurst #2, Spring 2005*, Distance Insights library report R45.
  [PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/df3cf55d-4c47-4dbd-8600-3c10d121d248/R45%20-%20Preliminary%20Results%20of%20Turf%20Testing%20at%20Pinehurst.pdf)
  - Per-shot rebounds on three fairways and three greens: `GREEN` and
    `FAIRWAY` compliance and sliding friction.
- USGA and The R&A (2007), *Second Report on the Study of Spin Generation*,
  Appendix D, "Results of Experiment on Golf Ball Impact with Green Surface"
  (test of 30 October 2006).
  [PDF](https://www.jga.or.jp/jga/html/rules/image/info_groove_second_report_spin_generation_en.pdf)
  - `FIRM_GREEN`: 57 impacts on a firm green at stimp 10–11 ft; the greens'
    stimp setting and plowing from its bounce-and-roll totals.

### Ground: run and roll

- USGA and The R&A (2020), *Agronomic Impacts on Bounce and Roll Distance*,
  Distance Insights library report R03 (tests of 2018).
  [PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/cdd48ad2-2ed1-4e67-ad53-1a6522110779/R03%20-%20Agronomic%20Impacts%20on%20Bounce%20and%20Roll%20Distance.pdf)
  - `FAIRWAY` rolling resistance from its fairway Stimpmeter reading
    (5 ft 8 in); drive run check (17.8 ± 13.3 yd over 624 shots).
- USGA and R&A Rules Limited (2021), *Proposed Bounce Model for Use in
  Evaluating Optimum Overall Distance*.
  [PDF](https://assets-us-01.kc-usercontent.com/c42c7bf4-dca7-00ea-4f2e-373223f80f76/250b3e34-9172-4ad1-991f-9d5560cbf832/Proposed%20Bounce%20Model%20for%20Use%20in%20Evaluating%20Optimum%20Overall%20Distance.pdf)
  - The 2019 low-speed range test as published there (roll-out
    88.3 − 1.79α yd at 10–36°): a check on the tilt law at shallow landings and
    on `FAIRWAY` plowing. Its proposed line (79.1 − 1.6α yd), the PGA TOUR
    runs and the Philadelphia Cricket Club lines it reports check the drive's
    run.
- Quintavalla, S. and Rogers, M. J. (2008), survey of bounce and roll at the
  2007 PGA TOUR, as cited in USGA and R&A Rules Limited (2021) above.
  - Why there is no rough preset: the one published figure, 8–9 yd less bounce
    and roll for drives finishing in rough.
- Burritt, A. and Henrikson, E. (2024), "Modeling golf ball roll with a
  variable coefficient of friction", *Sports Engineering* 27, 39.
  [doi:10.1007/s12283-024-00481-5](https://doi.org/10.1007/s12283-024-00481-5).
  Conference version: *Engineering of Sport 15*, Loughborough University
  repository, [doi:10.17028/rd.lboro.27051061.v1](https://doi.org/10.17028/rd.lboro.27051061.v1)
  (CC BY 4.0)
  - The form of `RollLaw::Exponential`; evidence that rolling resistance
    changes with speed.
- Hubbard and Alaways (1999), "Mechanical interaction of the golf ball with
  putting greens", in *Science and Golf III*.
  - The ~2 mm rest sinkage of a ball on a green, below the 5 mm bounce-to-roll
    height.

### Slopes and putting checks

- Penner, A. R. (2002), "The physics of putting", *Canadian Journal of
  Physics* 80(2), 83–96.
  [doi:10.1139/p01-137](https://doi.org/10.1139/p01-137)
  - Stimpmeter exit speed of 1.83 m/s.
- Holmes, B. W. (1986), "Dialogue concerning the Stimpmeter", *The Physics
  Teacher* 24(7), 401–404.
  [doi:10.1119/1.2342065](https://doi.org/10.1119/1.2342065)
  - The analysis of the device behind that exit speed.
- Brede, A. D. (1991), "Correction for slope in green speed measurement of golf
  course putting greens", *Agronomy Journal* 83(2), 425–426.
  [doi:10.2134/agronj1991.00021962008300020032x](https://doi.org/10.2134/agronj1991.00021962008300020032x)
  - Check: the 2ab/(a+b) slope correction of stimpmeter readings.
- Cross, R. (2024), "Trajectory of a golf ball on a sloping green", *European
  Journal of Physics* 45(3), 035002.
  [doi:10.1088/1361-6404/ad2c2d](https://doi.org/10.1088/1361-6404/ad2c2d)
  - Check: break angle across a slope independent of speed under constant
    rolling resistance.
- Whitlark, B. and Pringle, M. (2012), "Can't Stop Your Ball on the Green?",
  *USGA Green Section Record* 50(25).
  - Check: landing shots uphill and downhill on a 1.5% green.

### Sand

- Uehara, J. S., Ambroso, M. A., Ojha, R. P. and Durian, D. J. (2003),
  "Low-speed impact craters in loose granular media", *Physical Review
  Letters* 90, 194301.
  [doi:10.1103/PhysRevLett.90.194301](https://doi.org/10.1103/PhysRevLett.90.194301)
  - A ball dropped on loose sand comes to rest part-buried: no restitution.
- Baker, S. W. (2001), "True grit", *Greenkeeper International*, May 2001,
  34–36.
  - A golf ball fired into bunker sands stays in its mark.
- The R&A, *Golf Course Condition and Playability: Sustainability Guidelines
  and Support for Clubs*, §7 Bunkers.
  [PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/e6bdf791-b035-4cbc-89c9-a9ef12d81267/Course%20Condition%20and%20Playability%20Guidelines.pdf)
  - Bunker sand is chosen so a ball stays where it lands.
- Suo, B., Quillen, A. C. et al. (2023), "Subsurface pulse, crater and ejecta
  asymmetry from oblique impacts into granular media", *Icarus* 408, 115816.
  [doi:10.1016/j.icarus.2023.115816](https://doi.org/10.1016/j.icarus.2023.115816)
  - Speed kept after oblique impacts at 104 m/s: one bound of the sand
    example's `absorption`.
- Wright, E., Quillen, A. C. et al. (2020), "Ricochets on asteroids:
  Experimental study of low velocity grazing impacts into granular media",
  *Icarus* 351, 113963.
  [doi:10.1016/j.icarus.2020.113963](https://doi.org/10.1016/j.icarus.2020.113963)
  - Speed kept at 2–5 m/s, the other bound; depth driven by the normal
    velocity.
- Katsuragi, H. and Durian, D. J. (2007), "Unified force law for granular
  impact cratering", *Nature Physics* 3, 420–423.
  [doi:10.1038/nphys583](https://doi.org/10.1038/nphys583)
- Katsuragi, H. and Durian, D. J. (2013), "Drag force scaling for penetration
  into granular media", *Physical Review E* 87, 052208.
  [doi:10.1103/PhysRevE.87.052208](https://doi.org/10.1103/PhysRevE.87.052208)
  - Both: drag proportional to v², the basis of the `(1 − absorption)^(sin θ)`
    form.
- Fukumoto, T., Ebata, H., Sharma, I. and Katsuragi, H. (2026), "Depth and slip
  ratio dependencies of friction for a sphere rolling on a granular slope".
  [arXiv:2602.01652](https://arxiv.org/abs/2602.01652)
- De Blasio, F. V. and Saeter, M.-B. (2009), "Rolling friction on a granular
  medium", *Physical Review E* 79, 022301.
  [doi:10.1103/PhysRevE.79.022301](https://doi.org/10.1103/PhysRevE.79.022301)
  - Both: rolling resistance on sand, the example's `mu_roll` 0.3.

### Aerodynamics and flight

- Bearman, P. W. and Harvey, J. K. (1976), "Golf ball aerodynamics",
  *Aeronautical Quarterly* 27(2), 112–122.
  [doi:10.1017/S0001925900007617](https://doi.org/10.1017/S0001925900007617)
- Li, J., Tsubokura, M. and Tsunoda, M. (2017), "Numerical investigation of the
  flow past a rotating golf ball and its comparison with a rotating smooth
  sphere", *Flow, Turbulence and Combustion* 99(3–4), 837–864.
  [doi:10.1007/s10494-017-9859-1](https://doi.org/10.1007/s10494-017-9859-1)
- Lyu, B., Kensrud, J., Smith, L. and Tosaya, T. (2018), "Aerodynamics of golf
  balls in still air", *Proceedings* 2, 238 (12th Conference of the
  International Sports Engineering Association).
  [doi:10.3390/proceedings2060238](https://doi.org/10.3390/proceedings2060238)
- Crabill, J., Witherden, F. D. and Jameson, A. (2019), "High-order
  computational fluid dynamics simulations of a spinning golf ball", *Sports
  Engineering* 22, 9.
  [doi:10.1007/s12283-019-0300-y](https://doi.org/10.1007/s12283-019-0300-y)
  - These four: published drag and lift the fitted coefficients are compared
    with; Lyu et al. also give the spin decay rate.
- Tutelman, D. (2016), "Spin Decay for a Ball in Flight", with the TrackMan
  newsletter of October 2010 quoted there.
  [tutelman.com](https://www.tutelman.com/golf/ballflight/spinDecay.php)
  - Spin decay of about 4% per second.
- Justus, C. G. and Mikhail, A. (1976), "Height variation of wind speed and
  wind distributions statistics", *Geophysical Research Letters* 3(5),
  261–264. [doi:10.1029/GL003i005p00261](https://doi.org/10.1029/GL003i005p00261)
  - Shear exponent of `Wind::PowerLaw`.
- TrackMan, "Introducing updated tour averages" (PGA TOUR averages, 2024).
  [trackman.com](https://www.trackman.com/blog/golf/introducing-updated-tour-averages)
  - Held-out check of carry by club.

Fit quality against each dataset: see [docs/FIT.md](docs/FIT.md).

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
