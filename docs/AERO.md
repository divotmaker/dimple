# Golf Ball Aerodynamics Model

Aerodynamic model behind `dimple`, for any project that needs realistic ball
flight from launch monitor data.

## Overview

A shot is flown by one model value, an `AeroModel`: the ball and the air it
flies through. dimple's model, `Tour`, computes the drag coefficient (Cd)
and lift coefficient (Cl) of a modern tour ball as functions of spin ratio,
in a given air density and wind; dimple adds gravity and integrates the
trajectory with RK4 at a 1 ms fixed timestep, with exponential spin decay.

**Key inputs:** ball speed, launch angle, launch azimuth, backspin RPM,
sidespin RPM, and the air: a density (from altitude and temperature with
`standard_air_density`, or measured) and an optional wind.

**Key outputs:** carry distance, apex height, lateral deviation, flight time,
descent angle, full trajectory points; with a ground model, every impact,
the hops, slide and roll, and the rest position (see
[Ground Model](#ground-model)).

## Conventions

Internal units are SI: meters, seconds, kilograms, radians, rad/s. Results
are SI too, with the descent angle in degrees; `units` converts to yards and
feet.

Coordinate frame (right-handed, origin at the tee):

| Axis | Direction |
|---|---|
| +X | downrange, along the target line |
| +Y | up |
| +Z | right of the target line, from the golfer's view |

All lateral quantities are **positive = right** of the target line:

| Quantity | Positive means |
|---|---|
| launch azimuth | ball starts right of target |
| sidespin | ball curves right (fade/slice for a right-handed golfer) |
| lateral deviation (`lateral_m`) | ball lands right of target |
| `z` of any position / velocity | right of / moving right |

Spin is a single angular velocity vector ω (right-hand rule), from launch to
rest. With launch azimuth `az` and launch heading `d = (cos az, 0, sin az)`,
backspin is about the horizontal axis perpendicular to the start line,
`d × ŷ`, and positive sidespin (clockwise viewed from above) about -Y:

```
ω = backspin · (−sin az, 0, cos az) − sidespin · ŷ
```

For a straight launch this is `ω = (0, -sidespin, backspin)`. The Magnus
direction ω̂ × v̂ then stays in the vertical plane of the start line for
backspin at any azimuth (a push with no sidespin flies straight along its
start line) and points toward +Z (right) for positive sidespin. With `Tour`
the axis is fixed in flight and only the rate decays; `FlightResult`
reports the landing spin vector.

## Force Model

Three forces act on the ball at each timestep:

```
F_total = F_gravity + F_drag + F_magnus
```

- **Gravity:** −g in Y, with g the model's `AeroModel::gravity()`
  (standard gravity, 9.80665 m/s², unless a model overrides it). The ground
  uses the same g.
- **Drag:** opposes the velocity relative to the air. Fd = -Cd × q × A / m × v̂
- **Magnus lift:** perpendicular to that velocity, in the spin plane. Fl = Cl × q × A / m × (ω̂ × v̂)

Where q = ½ρv² (dynamic pressure, v the airspeed), A = cross-sectional
area, m = ball mass. Drag and lift are what `AeroModel::acceleration`
returns; gravity is added by dimple.

## Cd/Cl Model

The coefficients depend on the spin ratio alone, SR = ω × R / V (ω the spin
rate, R the radius, V the airspeed). They do not depend on the Reynolds
number: the reference data are fitted best without a drag crisis, so
neither the air's viscosity nor its temperature enters except through the
density.

### Drag coefficient

```
Cd = Cd_0 + Cd_spin × SR / (1 + SR / SR_sat)
```

Cd_0 is the drag at zero spin. Spin adds drag through a saturating term:
near-linear (slope Cd_spin) below SR ≈ SR_sat, capped at Cd_spin × SR_sat
at high spin ratio. Spin drag is dominated by induced drag from Magnus lift,
and since Cl saturates in SR, the drag it induces saturates too.

### Lift coefficient

```
Cl = Cl_max × (1 - exp(-SR / SR_scale)) + Cl_bump × u × exp(1 − u)

where:
  u = (SR / SR_bump)⁴
```

The first term is an exponential saturation curve in spin ratio. This means:
- **Low SR (driver, ~0.074):** Cl is well below the ceiling (0.141)
- **High SR (7-iron, ~0.29):** Cl approaches the ceiling (0.31)
- **Very high SR (wedge, ~0.45):** Cl is nearly saturated (0.34)

This spin-sensitive lift rewards the "high launch, low spin" fitting
strategy that dominates modern equipment.

The second term is a small low-spin lift term (**beta**). It is zero at zero
spin, peaks at Cl_bump where SR = SR_bump, and fades out by SR ≈ 0.045 (it
adds 1.7×10⁻⁴ there and under 10⁻⁹ from SR 0.06). It lifts fast, low-spin
shots (under about 1,500 rpm at 170 mph) toward the reference data, which
show a ramp and plateau in carry there that the saturation curve alone
misses. Lift rises with spin ratio everywhere, and goes to zero smoothly as
the spin does: a ball with no spin gets no lift.

## The `Tour` model

`Tour` is a modern tour ball (USGA conforming) in a given air:
`Tour::new(air_density, wind)`, or `Tour::SEA_LEVEL` for 1.225 kg/m³ and
still air. Its constants are fixed and not part of the API; the drag and
lift coefficients are **beta** and may change between releases.

| Parameter | Value | Role |
|---|---|---|
| mass | 0.04593 kg | USGA max (1.62 oz) |
| diameter | 0.04267 m | USGA min (1.68 in) |
| inertia factor | 0.4 | I = k·m·R², uniform solid sphere (ground phase only) |
| Cd_0 | 0.10 | Drag at zero spin |
| Cd_spin | 1.69 | Spin-drag slope (saturating) |
| SR_sat | 0.24 | Spin-drag saturation constant (cap ≈ 0.41) |
| Cl_max | 0.359 | Lift ceiling |
| SR_scale | 0.149 | Lift saturation rate |
| Cl_bump | 0.0185 | Height of the low-spin lift term (beta) |
| SR_bump | 0.027 | Spin ratio where the low-spin term peaks (beta) |
| spin decay | 0.04 s⁻¹ | see [Spin Decay](#spin-decay) |

Coefficients at a few spin ratios:

| SR | 0 | 0.03 | 0.074 (driver) | 0.10 | 0.29 (7-iron) | 0.45 (wedge) |
|---|---|---|---|---|---|---|
| Cd | 0.100 | 0.145 | 0.196 | 0.219 | 0.322 | 0.365 |
| Cl | 0 | 0.082 | 0.141 | 0.176 | 0.308 | 0.341 |

Published aerodynamic measurements for comparison:
- Li, Tsubokura & Tsunoda (2017), *Flow Turb. Combust.* 99(3) — LES CFD at
  Re=1.1×10⁵: Cd 0.217 (static), Cd 0.240 and Cl 0.135 spinning at SR 0.1.
  The model at SR 0.1: Cd 0.219, Cl 0.176.
- Lyu, Kensrud, Smith & Tosaya (2018), *ISEA Proc.* — trajectory fit on 13
  production balls, Cd about 0.20 on average at high Reynolds number
- Crabill (2019), *Sports Eng.* 22:1-9 — high-order CFD at Re=1.5×10⁵: Cd 0.247
- Bearman & Harvey (1976), *Aero. Q.* 27:112-122 — original wind tunnel study

### Calibration

The constants are fitted to per-shot carry, apex height and descent angle
in standard sea-level air (1.225 kg/m³), from two kinds of reference:

1. **A reference sweep at sea level:** 41 launches from 70 to 170 mph,
   straight and curving, with the carry, apex, descent angle and offline of
   each.
2. **Two per-shot launch-monitor sessions:** a hybrid (11 shots) and an
   8-iron (13 shots), 91–133 mph.

The low-spin lift term is fitted on 11 further reference shots at
150–185 mph and 9–15° launch below a spin ratio of 0.045. Every other set
below is held out.

### Fit quality

Carry and apex errors are model minus reference: bias / RMSE.

| Reference (shots) | Fitted? | Carry % | Apex ft | Descent ° |
|---|---|---|---|---|
| Reference sweep, sea level (41) | fit | −0.2 / 0.8, worst −3.6 | −0.8 / 2.2 | −0.1 / 0.5 |
| Launch monitor, hybrid (11) | fit | −0.3 / 0.3 | +4.1 / 4.1 | — |
| Launch monitor, 8-iron (13) | fit | +0.4 / 0.5 | +0.6 / 1.0 | — |
| Reference high-speed sweep 150–185 mph (27) | low-spin term only | +0.0 / 0.5 | −1.0 / 2.0 | −0.1 / 0.7 |
| Reference, all sea-level shots (112) | partly | −0.1 / 0.6, worst −3.6 | −0.9 / 1.9 | −0.2 / 0.6 |
| Reference sweep at 5,000 ft (41) | held out | −0.2 / 0.8 | −0.5 / 2.0 | — |
| Reference replay of the hybrid and 8-iron sessions (26) | held out | −0.2 / 0.4, −0.0 / 0.3 | −1.2 / 1.8, +0.1 / 0.9 | −0.5 / 0.7, −0.2 / 0.4 |
| Launch monitor, driver (12) | held out | +1.0 / 1.5 | +5.4 / 5.5 | — |
| Reference replay of the driver session (12) | held out | +0.0 / 0.3 | −1.5 / 1.6 | −0.4 / 0.5 |
| A second reference set, sea level (41) | held out | +0.5 / 0.9 | +3.0 / 6.4 | +0.1 / 1.3 |
| Per-shot validation set, 26–153 mph, 0.8k–14k rpm (1,178) | held out | +0.1 yd / 1.6 yd | +0.3 / 2.0 | +0.3 / 1.0 |

On the 1,178-shot set the offline error is −0.0 / 0.8 yd.

**Altitude.** From sea level to 5,000 ft on a standard day the reference
carries gain 5.7% on average and the model's 5.8%, from +3.2% for the
70 mph launches to +6.4% at 130 mph; per launch the gains differ by 0.19
points RMS (worst +0.45). The second reference set keeps sea-level temperature
at altitude; at its own settings it gains 3.69%, 6.92% and 10.08% at
2,500, 5,000 and 8,000 ft, the model 3.67%, 6.99% and 10.39%.

**Environmental response of the preset driver** (`ShotInput::driver`):
Denver (1,609 m, 17 °C) +7.5% carry (Titleist's ~6.1% assumes a standard
day); 90°F vs 50°F at sea level +9.1 yd (published ~8 yd, which includes
the ball's own response to temperature).

**Trackman 2024 PGA Tour averages** (25 °C sea level, held out): carry
bias +4.2 yd, RMSE 7.6 yd; descent angle bias −4.5°, RMSE 4.9°.

| Club | Model | Target | Δ | DA° (model) | DA° (TM) | DA Δ |
|---|---:|---:|---:|---:|---:|---:|
| Driver | 296.6y | 282.0y | +14.6y | 38.8 | 39.0 | -0.2 |
| 3-Wood | 264.4y | 249.0y | +15.4y | 42.2 | 44.0 | -1.8 |
| 5-Wood | 246.6y | 236.0y | +10.6y | 43.8 | 48.0 | -4.2 |
| Hybrid | 231.8y | 231.0y | +0.8y | 43.3 | 49.0 | -5.7 |
| 3-Iron | 226.0y | 218.0y | +8.0y | 41.7 | 48.0 | -6.3 |
| 4-Iron | 214.2y | 209.0y | +5.2y | 41.9 | 49.0 | -7.1 |
| 5-Iron | 202.5y | 199.0y | +3.5y | 42.9 | 50.0 | -7.1 |
| 6-Iron | 189.3y | 188.0y | +1.3y | 45.3 | 50.0 | -4.7 |
| 7-Iron | 173.9y | 176.0y | -2.1y | 46.2 | 51.0 | -4.8 |
| 8-Iron | 162.8y | 164.0y | -1.2y | 47.0 | 51.0 | -4.0 |
| 9-Iron | 151.3y | 152.0y | -0.7y | 47.6 | 52.0 | -4.4 |
| PW | 137.0y | 142.0y | -5.0y | 48.8 | 52.0 | -3.2 |

The per-shot references agree with each other to about 1% and disagree
with the tour table the same way: no setting of this model family fits both.
The tour averages mix competition and range shots across many events, and
the table is not smooth in ball speed (the hybrid at 149 mph carries 231 yd,
the 3-iron at 145 mph 218 yd), so the fit follows the per-shot data.

### Known limitations

1. **Trackman's tour table:** the model carries the woods 11–15 yd longer
   and lands the irons 3–7° shallower than its averages
   (see above). Above 135 mph the fit rests on the reference sweep; the one
   held-out driver session carries 1.0% short of the model with a 5 ft
   lower apex.
2. **Low spin (beta):** the low-spin lift term is fitted on a small set of
   150–185 mph shots at 9–15° launch below a spin ratio of 0.045, a range no
   shot in the 1,178-shot set reaches (its lowest launch spin ratio is
   0.048). Below a spin ratio of about 0.01 (under about 350 rpm at 170 mph) the
   model is extrapolated, and its zero-spin drag (0.10) is lower than wind
   tunnels measure for a non-spinning ball (about 0.2). A ball with no spin
   gets no lift, so an unmeasured spin must not be passed as zero.
3. **Slow, high-spin shots** (under 80 mph, or above a spin ratio of 0.4)
   carry about 1.5% short on average, up to 3.6%.
4. **Ball COR vs temperature unmodeled:** only the air density responds to
   temperature; cold-ball distance loss is understated.
5. **Single ball model.** Real balls vary ±5% in drag across dimple designs.
6. **Spin on the ground: modeled, lightly validated.** The full spin vector
   is carried through every impact, slide and roll (see
   [Ground Model](#ground-model)), so a curving ball kicks further off-line
   at the first bounces and spin along the direction of travel drifts the
   ball while it slides. Gaps: no published measurement of these lateral
   effects to validate against; the pivot (torsional friction) lever is
   uncalibrated; the fairway's compliance is chosen, not measured; landing
   shots on slopes are beta (see [Ground limitations](#limitations)).

### Your own model

Any type implementing `AeroModel` flies a shot: it returns the air's
acceleration on the ball (without gravity) from its state (time, position,
velocity, spin) and its height above the ground, advances the spin over a
flight step, and gives the radius and inertia factor the ground uses. A
model can also override `gravity()`; dimple then uses that value in flight
and in every ground contact. A ball with drag only:

```rust
use dimple::{AeroModel, BallState, ShotInput, Surface, Vec3, simulate_shot};

struct DragOnly { k: f64 }

impl AeroModel for DragOnly {
    fn acceleration(&self, ball: &BallState, _height: f64) -> Vec3 {
        ball.vel * (-self.k * ball.vel.length())
    }
    fn spin_after(&self, ball: &BallState, _dt: f64) -> Vec3 { ball.spin }
    fn radius_m(&self) -> f64 { 0.021_335 }
    fn inertia_factor(&self) -> f64 { 0.4 }
}

let shot = simulate_shot(&ShotInput::driver(), &DragOnly { k: 0.004 }, &Surface::FAIRWAY);
```

A model is a pure function of its arguments (the integrator calls
`acceleration` four times per step), so it has no interior mutability.

## Air

`standard_air_density(altitude_m, temperature_c)` gives the density from
the standard pressure at the altitude and the actual temperature there, by
the ideal gas law:

```
ρ = ρ₀ × (1 - L×h/T₀)^(g×M/(R×L)) × T₀ / T
```

Where ρ₀ = 1.225 kg/m³, L = 0.0065 K/m, T₀ = 288.15 K, g×M/(R×L) = 5.2559,
h the altitude (floored at sea level) and T the temperature in kelvin. The
first two factors are ρ₀ times the standard pressure ratio at h; the
temperature is used as given, not lowered with altitude. On a standard day
(15 °C at sea level, 5.1 °C at 1,524 m) this is the ISA density: 1.225 and
1.0555 kg/m³.

### Measured air (beta)

`Tour::new(air_density, wind)` takes the density as given, with no
altitude inference and no clamp at sea level, so cold, dense air is
modeled as dense. A launch monitor with a weather station can pass its own:
with pressure p (Pa) and temperature T (K), ρ = p / (287.05 · T). Humidity
lowers the density slightly (about 0.8% at 30 °C and 50% relative
humidity). A density of 0 gives ballistic flight under gravity alone.

```rust
use dimple::{Tour, Wind};

// 101.8 kPa and 2 °C from a weather station.
let tour = Tour::new(101_800.0 / (287.05 * (2.0 + 273.15)), Wind::None);
```

## Wind (beta)

A `Tour`'s wind is the wind the shot flies through; `Tour::SEA_LEVEL` is
`Wind::None`, still air. Wind is **beta**: the model below follows the
physics but has not been compared with measured shots in wind. A wind is a
horizontal air velocity in the public frame, `+x` blowing downrange (a
tailwind) and `+z` blowing to the right. The built-in winds ignore a
vertical component; a custom wind (below) keeps it.

```rust
use dimple::{ShotInput, Surface, Tour, Vec3, Wind, simulate_shot};

// 5 m/s into the golfer's face.
let tour = Tour::new(1.225, Wind::Uniform(Vec3::new(-5.0, 0.0, 0.0)));
let shot = simulate_shot(&ShotInput::driver(), &tour, &Surface::FAIRWAY);
```

Drag and Magnus lift act on the ball's velocity relative to the air,
`v − wind`, both in flight and in the hops between impacts (see
[Hops](#hops-and-the-bounce-to-roll-transition)). The impacts themselves are
unchanged, because the turf does not move. A uniform wind is not a frame
shift of the still-air trajectory: the airspeed at launch changes, and with
it the drag and the spin ratio, so the shot is integrated again rather than
translated.

`Wind::Uniform` is the same wind at every height. For a 150 mph drive at 13°
and 2700 rpm with `Tour` at sea level (`examples/wind.rs`):

| Wind | Carry | Total | Lateral at rest |
|---|---:|---:|---:|
| Still | 229.5 m | 245.9 m | 0.0 m |
| 5 m/s headwind | 213.0 m | 223.1 m | 0.0 m |
| 5 m/s tailwind | 241.3 m | 260.6 m | 0.0 m |
| 5 m/s left to right | 230.0 m | 246.3 m | 19.9 m right |

The headwind costs more than the same tailwind gains, because drag grows
with the square of the airspeed.

### Wind profile

A reported wind is usually measured 10 m above the ground, and the wind is
weaker near the ground and stronger above. `Wind::PowerLaw { at_10m }`
scales the 10 m reading by (h / 10)^α, where h is the ball's height above the
ground directly beneath it, with the shear exponent of Justus & Mikhail
(1976), α = 0.37 − 0.0881 ln(U₁₀) for the 10 m speed U₁₀ in m/s. Below
ground level the wind is zero, and at or below 0.1 m/s the wind is applied
uniformly. The exponent is floored at zero, so above about 66.7 m/s the wind
is applied uniformly too. The profile is sampled at
every RK4 stage of the flight and of the hops, so the hops, a few
centimeters above the turf, fly through nearly still air. It is a
neutral-stability assumption, not a measured wind field.

The same drive and winds as above, as 10 m readings:

| Wind (10 m) | Carry | Total | Lateral at rest |
|---|---:|---:|---:|
| 5 m/s headwind | 210.4 m | 220.6 m | 0.0 m |
| 5 m/s tailwind | 241.6 m | 259.9 m | 0.0 m |
| 5 m/s left to right | 230.2 m | 246.6 m | 21.7 m right |

The drive spends most of its flight above 10 m, where the profile is
stronger than the reading, so the headwind costs 2.6 m more carry than
uniform wind. The tailwind's total is shorter despite a longer carry,
because the hops lose the help of the wind.

Justus, C. G. & Mikhail, A. (1976), "Height variation of wind speed and
wind distributions statistics", *Geophysical Research Letters* 3(5),
261–264, doi:10.1029/GL003i005p00261.

### Your own wind

`Wind::custom` takes any function of the shot time, the ball's position
and its height above the ground (a `WindField`; closures implement it),
for gusts or a wind field from a course map. It must be a pure function of
those arguments. A headwind that starts two seconds into the shot:

```rust
use dimple::{ShotInput, Tour, Vec3, Wind, simulate_flight};

let gust = Wind::custom(|t: f64, _pos: Vec3, _height: f64| {
    if t > 2.0 { Vec3::new(-8.0, 0.0, 0.0) } else { Vec3::ZERO }
});
let flight = simulate_flight(&ShotInput::driver(), &Tour::new(1.225, gust));
```

A custom wind that returns a constant horizontal vector gives exactly the
shot of the same `Wind::Uniform`. A custom wind is used with all three
components, so it can blow up or down: a 3 m/s updraft raises the apex.

## Spin Decay

`Tour`'s spin decays exponentially during flight, about a fixed axis:

```
ω(t) = ω₀ × exp(-λ × t)
```

λ = 0.04 s⁻¹ (~21% loss over 6 seconds). Published data:
- Trackman (2010): ~4%/s → λ=0.04 → 78% remaining at 6s
- Lyu et al. (2018): 4%/s trajectory fits
- Tutelman: τ=30s → λ=0.033 → 82% at 6s

## Ground Model

**Status.** The impact model reproduces Penner's (2002) published bounce
heights and spin-back. On natural turf it is compared against the published
bounce measurements of Biber et al. (2023): on their natural tee `GREEN` sits at a scaled error of about 12% from their
fit to the data, where the paper's own physics models score 19–35% against the
data (see
[Validation against published bounce data](#validation-against-published-bounce-data)).
Rolling on slopes is checked against Brede's (1991) stimpmeter slope
correction and Cross's (2024) break-angle invariance (see
[Stimpmeter](#stimpmeter) and [Terrain and slopes](#terrain-and-slopes)).
A tour-average drive runs about 18 yd on `FAIRWAY`, between the USGA/R&A
proposed bounce model and the ShotLink median for its landing angle (see
[Run-out on a fairway](#run-out-on-a-fairway)).

**Beta:** landing shots on slopes, which agree with Whitlark & Pringle's
(2012) study of a 1.5% green only at the level of its published means; and
the run and total distances of full shots, the onset speed of the tilt
fade, `FAIRWAY`'s compliance, rolling resistance above 3 m/s, the pivot
lever and the lateral ground effects. Only the drive's run is compared
with measured fairway data, and no published measurement of wedge, pitch
or chip run on a fairway exists to compare against. The presets may be
recalibrated when such data exist.

Everything after landing is a rigid sphere with point contact, carrying the
full spin vector ω through every contact. The tracked point is the ball's
center, and the ball touches the ground when its center reaches the ground
height (the radius is not added), as in the flight's landing; on ground
level with the tee a ball touching the ground has `y = 0`. The ground is a
heightfield with slopes (see [Terrain and slopes](#terrain-and-slopes)).

### What spin can and cannot do on the ground

The contact point is at `r = −R n` (n the ground normal, R the radius) and
slips at `u = v_t + ω × r`; friction acts against `u`.

- **Spin about the normal can't move the contact point**, so it can't kick
  the ball; pivot (torsional) friction only decays it. On level ground this
  is the launch-monitor *sidespin* component, ω_y. Sidespin therefore does
  not steer a sliding or rolling ball on level ground. On a slope the
  normal is tilted, so part of ω_y lies in the contact plane and slips the
  contact point: on rigid gripping ground it kicks the ball by
  `(k/(1+k))·R·ω_y·(ŷ × n)`, along the line of play on a side slope and
  across it on an uphill or downhill one.
- **A curving ball kicks further off-line at its first bounces.** The spin
  axis is fixed in flight, perpendicular to the start line, but a curving
  ball lands heading further off-line than that. The backspin then slips
  the contact point partly across the heading, and friction removes more
  velocity along the start line than across it: a fade kicks right, a draw
  left. In the grip regime on rigid ground the heading after the bounce is
  exactly `atan((5/7)·v_z / ((5/7)·v_x − (2/7)·R·ω_b))` (backspin ω_b about
  +Z). The effect is **modest**: the driver preset with ±600 rpm sidespin
  lands heading 12.4° off-line and leaves its first bounce at 13.8°.
- **Spin along the direction of travel** ("rifle spin") slips the contact
  point sideways, so it drifts the ball laterally while it slides.
- **Penner's tilt adds a small throw the other way.** Each impact is against
  a tilted effective plane (below), so a fraction sin β of ω_y lies in the
  contact plane and slips the contact point sideways at impact: a ball
  spinning clockwise from above (fade sidespin) is thrown left, like side
  spin off a cushion. Pivot friction then clears ω_y, so the throw acts
  mainly at the landing impact. For the ±600 rpm driver on `FAIRWAY` at
  sea level, measured as lateral run across the landing heading, the
  heading-mismatch kick alone (the same landing with ω_y removed) adds
  about 0.9 m, the throw takes back about 0.1 m, and the net is about
  +0.8 m (0.9 yd) on a 17 m run. The fade comes to rest 32.1 m (35.1 yd)
  right of the target line, 4.3 m right of its landing point.

So dimple does not convert sidespin into a lateral nudge; it carries the
spin vector through the contact physics, and these effects follow from it.

### Impact

Each impact follows Penner (2002), "The run of a golf ball", *Can. J.
Phys.* 80, 931–940, in vector form. The ball deforms the turf and meets the
rising lip of its own crater, so the impact is treated as one off a rigid
plane inclined by an effective contact angle β toward the incoming ball.
With `v_n = −v·n` the approach speed, `v_t` the tangential velocity
(direction `t̂`), `v = |v|` and `φ = atan(|v_t| / v_n)` the incidence angle
from the normal:

```
θc   = 15.4° · (v / 18.6 m/s) · (φ / 44.4°)        (Penner eq. 8)
β    = min(compliance · θc · min(1, v_n / 2 m/s), φ)
n'   = n cos β − t̂ sin β                          (effective normal)
e    = restitution at |v·n'|                       (Penner eq. 5)
j_n  = (1 + e) · |v·n'|                            (normal impulse / m)
u    = (v − (v·n')n') + ω × (−R n')               (contact slip)
j_t  = min(|u| · k/(1+k), μ · j_n)                 (along −û)
v   += j_n n' − j_t û
ω   += (−R n') × (−j_t û) / (k R²)
Δω_n = −sign(ω_n) · min(|ω_n|, pivot · j_n / (k R²))   about n
```

`k` is the ball's inertia factor (`I = k m R²`, `AeroModel::inertia_factor`;
2/5, a uniform solid sphere, for `Tour`) and R its radius (`radius_m`).
Friction either brings the contact point to rest (grip: the ball leaves
rolling on the effective plane) or slides throughout at the Coulomb limit.
No impulse is applied when the ball is not
moving into the ground (grazing).

The factor `min(1, v_n / 2 m/s)` fades the tilt in for soft contacts; every
lofted landing (launch above about 12°) strikes faster than 2 m/s and gets
Penner's angle in full (see
[Hops and the bounce-to-roll transition](#hops-and-the-bounce-to-roll-transition)).

For a ball with no spin about the vertical traveling in a vertical plane,
this is exactly Penner's planar solution: the normal component
`v_n' = v_n cos β + v_t sin β` (eq. 6), Daish's rigid-surface rebound in
the effective frame (eqs. 2–4: grip `(5/7)·v_t' − (2/7)·R·ω_back`, or slip
with friction impulse `μ(1+e)·v_n'`), rotated back (eq. 7). The test suite
checks the vector form against the planar solution, with the fade,
impact by impact.

The coefficient of restitution falls with normal impact speed (Penner eq. 5,
his fit to his own turf measurements):

```
e = 0.510 − 0.0375·v_n' + 0.000903·v_n'²     v_n' ≤ 20 m/s
e = 0.120                                    v_n' > 20 m/s
```

`v_n'` is the normal speed against the effective plane, as in Penner. The two
branches meet within 0.0012 at 20 m/s. `Restitution::Constant(e)` uses a
fixed value instead.

The tilt steepens every forward rebound; with `e = 1, μ = 0` it rotates the
flat-ground mirror direction up by exactly 2β. β is capped at the incidence
angle, where the ball meets the effective plane head-on; this is a guard for
extreme `compliance` values, and the presets never reach it. Backspin can
send the ball backward. A rebound that points into the ground stays on the
ground and rolls.

**Pivot friction.** The friction moment about the surface normal is
`pivot` × normal force. It acts about the true normal `n`, not the tilted
`n'`, so it removes only spin about the vertical and never the spin along
the heading that drives a curving ball's kick. The presets use 1 mm: the
full-slip torsional couple of a Hertzian contact, (3π/16)·μ·a (Johnson,
*Contact Mechanics*, 1985), with μ = 0.4 and a contact radius a ≈ 4 mm. It
is uncalibrated.

Pivot friction never changes the velocity of the impact it acts in, nor a
slide or roll on level ground. It reaches the path only through spin about
the normal carried into a later impact, where the tilted plane turns part of
it into slip, and through the Magnus force that spin gives the hops in
between. With the 1 mm lever a full-shot landing removes all of it in one
impact: a driver landing on `FAIRWAY` can remove up to about 160 rad/s
(~1500 rpm), and a ±600 rpm-sidespin drive lands with 47 rad/s (450 rpm).
So on level ground the path is the same for every lever from 1 mm up. With
no pivot friction at all the spin survives to the later impacts, whose
throws move that drive's rest by about 11 cm; on rigid ground (no tilt) it
survives only into the hops, which it curves about 6 cm further the way the
sidespin curved the flight.

### Hops and the bounce-to-roll transition

As in Penner, the impact model applies to every impact, not just the
first. After the landing the ball makes a few progressively smaller hops,
then slides and rolls. Penner's run model, hops plus roll, agrees fairly
well with empirical drive runs (Cochran & Stobbs; Williams; his Fig. 5),
with the model's runs slightly longer.

**The tilt fade.** Penner's θc is a first-order extrapolation, linear in
speed and incidence angle, from a single measured green impact at
18.6 m/s, 44.4° from the vertical. Taken in full at every impact, it also
tilts the plane for the low, near-grazing hops after the first bounce:
arriving at 0.3–1 m/s of normal speed and 2–4 m/s along the ground, 75–85°
from the vertical, the ball gets 3–5° of tilt, which turns forward speed
into a rebound nearly as high as the one it arrived from. The hops then
settle at 5–30 mm and skim on for 20–40 impacts, losing 1–2% of their
forward speed each, until the ball slows to about 2 m/s. On a fairway they
cover most of the run, the fairway's rolling resistance barely matters, and
a fairway runs as far as a green. The tilt stands for the lip of the crater
the ball presses into the turf, and a contact that soft presses almost
none; Penner also expects his bounce model to break down as bounce heights
approach the ~2 mm a ball sits into a green (see below). So dimple fades the
tilt in linearly with the normal impact speed `v_n`:
`β = min(compliance · θc · min(1, v_n / 2 m/s), φ)`.

- Lofted landings (launch above about 12°) strike at 4 m/s (a low chip) to
  20 m/s (full shots) of normal speed and keep Penner's angle in full, so
  his published bounce heights and drive runs still reproduce (see
  [Surfaces](#surfaces)). A bump-and-run launched under about 12° lands
  inside the fade and skims 10–16 impacts before it rolls.
- On level ground the rebound of a low hop then leaves 0.5–0.8 times as
  fast as it arrived while the ball moves along the ground at 2–4 m/s, so
  the hops die out geometrically (see the slope analysis below): a full
  shot from the bag makes 4–10 impacts (a drive about 9, an iron, pitch or
  wedge 4–8), then slides and rolls, from about 2–4 m/s after a drive,
  iron, pitch or chip and about 1 m/s after a tour wedge.
- The 2 m/s onset is chosen, not measured. Runs barely depend on it between
  1.5 and 3 m/s (in 25 °C air the tour drive of
  [Run-out on a fairway](#run-out-on-a-fairway) runs 19.3 to 17.5 yd, its
  amateur lob wedge 6.2 to 6.0 yd), nor on the shape of the ramp: a
  quadratic or smoothstep fade moves them by about a yard or less. A fade
  weaker than linear, such as `√(v_n / 2 m/s)`, leaves most of the hop
  train in place.

The ball keeps bouncing while its rebound would rise at least 5 mm above the
ground, measured along the normal: `(v·n)² / (2·g·n_y) ≥ 5 mm` (on level
ground `v·n` ≥ ~0.313 m/s), Penner's transition height; he tested 1–10 mm
and found only a minor effect on the run. On a slope the hops must also die
out (below). Otherwise its normal velocity is dropped and it slides or
rolls. A safety bound of 200 impacts guards degenerate surfaces (e.g.
elastic and frictionless); real shots end by these rules well before it.

**Flight between impacts.** Each hop flies through the same air as the
shot: `GroundModel::settle` takes the shot's model, so the hops fly
through its air and wind, the wind sampled at the hop's height above the
terrain, and their drag and lift act on the velocity relative to the
moving air. The forces are the flight's: the model's `acceleration`
(`Tour`: drag and Magnus lift) plus gravity, the model's `gravity()`,
integrated with RK4 at 1 ms. The spin is held fixed over each hop: its axis
and rate set the lift, and the ball arrives at the next impact with exactly
the spin it left the last one with. `Tour`'s coefficients depend on the
spin ratio alone, so they hold at hop speeds (0.3–8 m/s) as in flight. A
ball leaving an impact is close to rolling, spin ratio about 1, which gives
Cd ≈ 0.43 and Cl ≈ 0.36. At sea level drag then decelerates the ball by
about 0.52 m/s² at 8 m/s, 0.13 m/s² at 4 m/s and 0.03 m/s² at 2 m/s, and the
forward topspin it leaves with pushes it down with about 0.44 m/s² (4.5% of
g) at 8 m/s. Small as these are, a drive spends nearly 3 s in its hops:
through sea-level air a `FAIRWAY` drive makes 9 impacts and runs 17.4 yd
past the landing point, against 11 impacts and 19.5 yd with ballistic hops
(`GREEN`: 23.2 yd; ballistic 25.5 yd). Irons and wedges, which hop little,
barely change. A model with no air (`Tour::new(0.0, Wind::None)`) gives
ballistic hops under gravity alone, the model Penner uses; in dimple they
lengthen a drive's run by about 2 yd.

The spin is not decayed during hops: dimple does not call the model's
`spin_after` between impacts. `Tour`'s flight decay rate (0.04 s⁻¹) is
calibrated at 50–75 m/s; an
aerodynamic spin-down torque with a moment coefficient linear in spin ratio
gives a rate proportional to speed, roughly ten times lower at hop speeds,
and no measurement of it at those speeds is used here.

Each hop step is checked for a ground crossing at 8 equal sub-steps and the
crossing is solved within the first sub-step that ends at or below the
ground, as the root of the clearance `y(τ) − height(x(τ), z(τ))` (regula
falsi with bisection, to 1e-12 m). The state τ into the step is the cubic
Hermite interpolant of the step's end states (position and velocity at both
ends, so exact at both ends of the step); under gravity alone it is the
exact parabola, which RK4 reproduces. The next impact uses the terrain's
normal there.

**Hops on a slope.** Downhill, gravity hands energy back between hops, and
the 5 mm height alone need not end them. At grazing incidence (φ ≈ 90°)
eq. 8 makes Penner's tilt proportional to the speed, `C·v` with
C = compliance · 15.4° · (90°/44.4°) / (18.6 m/s) (0.029 s/m for the
presets), and below 2 m/s of normal speed the fade scales it by w / 2 m/s,
w the normal speed: `β ≈ C·v·w / (2 m/s)`. Take a ball moving along the
ground at speed v and hopping with normal speed w. Against the tilted plane
the impact rebounds at `w' = e·w + (1 + e)·vβ` (normal restitution; the
grip friction loss is of higher order), so each impact scales a low hop's
normal speed by

```
q = w' / w ≈ e + (1 + e) · C·v² / (2 m/s)
```

With e ≈ 0.5 at these slow impacts, q < 1 below v_c ≈ 5 m/s (q ≈ 0.6 at
2 m/s and 0.8 at 4 m/s, as the full impact model gives). On level ground,
where every impact also slows the ball, the hops therefore shrink through
the 5 mm transition within a few impacts; a ball skimming faster than v_c
first hops higher, until its impacts have slowed it. Down a slope θ, each hop of time
`T = 2w/(g cos θ)` speeds the ball up by about `2w·tan θ`, and a slope that
keeps it moving at v_c or faster keeps its hops from shrinking at any
height: they hold or grow, and the ball never rolls. Drag at these speeds is
too weak to prevent it. Under the 5 mm rule alone, a drive landing on a 5°
downhill `FAIRWAY` is still hopping after 200 impacts (with ballistic hops,
already on 3.75°).

So on a slope, or in wind, the ball keeps hopping only while the train dies out. Each
hop is flown first and then judged. With `E = |v|²/2 + k·R²·|ω|²/2` per
unit mass, the energy the impact removed `ΔE = E_in − E_out` and the
energy the hop, flown through whatever air there is, handed back
`ΔU = E_next − E_out` (the potential energy released down the slope less
what drag took; the spin is held over the hop), it keeps hopping if

```
ΔU ≤ 0   or   ΔE − ΔU ≥ g · min_hop_height
```

Each impact must take at least the energy of a 5 mm hop more than the next
hop gives back. Otherwise the hop is discarded and the ball goes to the
ground at the impact, exactly as below the 5 mm height. On level ground in
still air the 5 mm height alone decides, so those results do not depend on
this rule, and a hop that climbs returns nothing.

Wind is the other source. Hops are flown through the wind, and a ball
moving slower than the air is pushed along by it, so on level ground too a
hop can hand back energy (`ΔU` then also holds the wind's work). Without
the rule a fast landing in a strong wind settles into a train of hops that
the wind refills as fast as the Penner impacts drain them (the tilted
contact plane turns the regained speed back into bounce), and the ball hops
on to the impact limit: hundreds of meters past the landing, with hops
meters high. Judged the same way as on a slope, such a train goes to the
ground once the impacts stop draining it; a wind strong enough at the turf
may still blow the ball some way along it.

This stays within Penner's transition. He notes that "the transition
between bouncing and rolling will, of course, take place over several
bounces and the 5 mm mark, which was used in the analysis, must be treated
as a simplification of the actual behaviour of a real golf ball" (p. 935),
and, since a ball at rest sits about 2 mm into a green (Hubbard & Alaways),
that "it is expected that the bounce model will break down as bounce
heights approach the 2 mm value" (pp. 934–935). A train the slope keeps
from shrinking is not a transition that decays through 5 mm; it is the
bounce model run indefinitely where it breaks down, at grazing impacts
where the linear θc law is least certain. Once the impacts no
longer clearly outrun the slope the ball rolls, and the rolling
resistance stands for the turf's dissipation at that speed.

### Slide and roll

On the ground the state is position, velocity and the spin vector,
integrated with RK4 at 1 ms within each smooth phase, with the terrain
sampled under the ball every step. With `n` the normal and g the model's
`gravity()`, the normal load is `m·g·n_y` and the slope pulls with the
tangential gravity `g_t = g − (g·n)n` (zero on level ground):

- **Slide** (the contact point slips): friction `−μ·m·g·n_y·û`, its torque
  `(−R n) × F`, and `g_t`. On level ground the slip direction is fixed and
  the slip speed falls linearly at `((1+k)/k)·μ·g` (3.5·μ·g for k = 2/5),
  so the slide is stepped exactly to the moment the slip reaches zero. On a
  slope `g_t` turns the slip direction (`du/dt = −((1+k)/k)·μ·g·n_y·û + g_t`),
  so û is re-evaluated at every RK4 stage. The moment the slip falls to
  `roll_threshold` (0.01 m/s) is solved within the step (regula falsi with
  bisection), and from there the slide finishes along that slip direction
  to zero slip, which is exact when the slip lies along the fall line.
  Backspin can reverse the ball.
- **Roll**: the spin tangential to the ground is locked to the velocity, and
  the ball decelerates at `μ_eff·g·n_y/(1+k)` against its velocity, with a
  speed-dependent rolling resistance, while the slope pulls with
  `g_t/(1+k)`. `mu_roll` is the resistance at zero speed and
  `Surface::roll_law` (a `RollLaw`) its speed dependence,
  `μ_eff = mu_roll · f(v)` with f(0) = 1:

  ```
  RollLaw::Constant                 f = 1
  RollLaw::Quadratic { v_ref }      f = 1 + (v / v_ref)²
  RollLaw::Exponential { v_ref }    f = exp(v / v_ref)
  ```

  `v_ref` must be nonzero (±∞ is the constant law); a zero or NaN `v_ref`
  rolls as `Constant`, and the stimpmeter helpers reject it (see
  [Stimpmeter](#stimpmeter)).

  The presets use `Quadratic { v_ref: 2.5 }`, dimple's own empirical law,
  not a fit to published data: 2.5 m/s is chosen so full-shot rollouts are
  plausible. At driver post-bounce speeds μ_eff is many times the base
  value; at stimpmeter speeds it is about 1.5×. Integrated, a ball that
  starts rolling at `v₀` stops after
  `(1+k)·v_ref² / (2·mu_roll·g) · ln(1 + v₀²/v_ref²)`: 3.48 m (stimp
  11.4 ft) on `GREEN` and 1.60 m (5.2 ft) on `FAIRWAY` for the
  stimpmeter's 1.83 m/s (see [Stimpmeter](#stimpmeter)).

  `Constant` is Penner's run model: his run paper uses a constant rolling
  friction, `a = −(5/7)·ρ_g·g` with ρ_g = 0.131 (eq. 9), which reads a
  stimp of about 6.0 ft (ρ = 7·v₀²/(10·g·d)).

  On slopes the published checks favor `Constant`: Brede's slope
  correction of stimpmeter readings is exact under it, and Cross's
  speed-independent break angle holds under it alone. Under the quadratic
  law both are approximate, and at the same stimp it lowers the grade at
  which a green holds a ball (`GREEN`: 5.5% against 6.9%). See
  [Stimpmeter](#stimpmeter) and [Terrain and slopes](#terrain-and-slopes).
  The presets keep the quadratic law, which keeps full-shot rollouts
  plausible.

  `Exponential` follows Burritt & Henrikson (2024), who fit an exponential
  function of speed to measured putts. They track whole ~20 ft putts with
  overhead cameras on four surfaces, solve for the rolling-friction
  coefficient frame by frame with Penner's roll model, and fit an
  exponential μ_roll(v) per surface. On average μ varies by about 40% over
  a putt. On three surfaces the ball resists more at speed; on one mat it
  resists less. dimple does not reproduce their exact parameterization or
  coefficients. Their fitted curves cover 0–3 m/s; values read from their
  Figure 1 (approximate, read from the figure, not tabulated values), with
  the `v_ref` that makes μ₀·exp(v/v_ref) pass through both readings:

  | Surface | μ at 0 m/s | μ at 3 m/s | equivalent `v_ref` |
  |---|---|---|---|
  | Mat A | ≈ 0.028 | ≈ 0.062 | ≈ +3.8 m/s |
  | Mat B | ≈ 0.057 | ≈ 0.030 | ≈ −4.7 m/s |
  | Practice green C | ≈ 0.027 | ≈ 0.175 | ≈ +1.6 m/s |
  | Practice fringe D | ≈ 0.155 | ≈ 0.285 | ≈ +4.9 m/s |

  So there is no single published rate: a surface fitted this way is
  `Surface { mu_roll: μ₀, roll_law: RollLaw::Exponential { v_ref }, .. }`.
  The presets' quadratic law gives μ(3 m/s)/μ(0) = 2.44, within that range
  (0.53–6.5) and next to Mat A's 2.2. Full shots reach the ground at up to
  about 4 m/s (see
  [Hops and the bounce-to-roll transition](#hops-and-the-bounce-to-roll-transition)),
  where the laws part: re-pinned to the same stimp, `Constant` lengthens the
  preset drive's run from 17.4 to 20.1 yd on `FAIRWAY` and from 23.2 to
  28.9 yd on `GREEN`; irons and wedges, which reach the ground slower, move
  by a few tenths of a yard. Re-pinned to the
  presets' stimp, `Exponential { v_ref: 7.0 }` gives `GREEN` μ₀ = 0.058
  and `FAIRWAY` μ₀ = 0.126; a roll from 5 m/s on `GREEN` runs 19.46 m
  (quadratic 13.05 m) and one from 10 m/s on `FAIRWAY` 23.18 m (quadratic
  10.53 m).
- **Spin about the normal** decays at `pivot·g·n_y / (k R²)` in both
  phases, saturating at zero.
- **Rest**: a rolling ball that is slowing and slower than one step's
  deceleration (or than `stop_speed`, default 0.01 m/s) finishes at
  constant deceleration and stops. At zero speed rolling resistance can
  hold the ball against the slope where `g·sin θ ≤ mu_roll·g·cos θ`, i.e.
  `tan θ ≤ mu_roll` (`GroundSample::holds_at_rest`, on the ground where
  the ball stops): up to 6.8° on `FAIRWAY` and 3.1° on `GREEN`. There it
  stays at rest; on steeper ground it starts rolling downhill again from
  zero speed (a ball rolling uphill stops once, then rolls back). The result reports why the ground phase
  ended: `Rest`, `MaxTime` (default 120 s after landing, which is also how a
  ball still rolling down a long slope ends) or `MaxImpacts`.

### Terrain and slopes

`Terrain::sample(x, z)` describes the ground as a heightfield: the height,
the unit normal (`n_y > 0`; for `y = h(x, z)` it is `(−∂h/∂x, 1, −∂h/∂z)`
normalized) and the surface at each horizontal position. `GroundPhysics`
samples it at the ball center's (x, z); the ball is on the ground when its
center is at that height.

- **Landing.** `simulate_shot` launches from `(0, height(0, 0), 0)`. Each
  1 ms flight step is checked for a crossing of the clearance
  `y − height(x, z)` at 8 equal sub-steps, and the crossing is solved
  within the first sub-step that ends at or below the ground (regula falsi
  with bisection, to 1e-12 m), so the ball lands where its flight first
  meets an elevated green, a downslope, curved ground or a ridge thinner
  than the step. A feature narrower than a sub-step (about 1 cm at 75 m/s)
  can still be missed. Carry is the horizontal distance from the tee to
  that point.
- **Impacts** use the sampled normal: the tilt, the slip and the pivot all
  act about the true normal of the slope.
- **On the ground** each step is taken in the tangent plane under the ball;
  the ball is then put back on the heightfield at its new (x, z), with its
  velocity projected onto the new tangent plane and, when rolling, its
  tangential spin re-locked to that velocity. On a plane this is exact. On
  curved ground the ball always follows the surface: over a crest where
  `v²·κ > g·n_y` (κ the curvature along the path) a real ball would leave
  the ground, which is not modeled, and entering a concave transition the
  velocity component along the new normal is removed, as in an inelastic
  contact.

Run past the landing point (yards) on planes rising (−) or falling (+)
downrange through the tee; `Rest` unless marked. `FAIRWAY`:

| Slope | −3° | −2° | −1° | 0° | +1° | +2° | +3° | +4° | +5° | +6° |
|---|---|---|---|---|---|---|---|---|---|---|
| driver | 15.5 | 16.1 | 16.7 | 17.4 | 18.6 | 19.9 | 21.0 | 22.2 | 24.6 | 27.6 |
| driver, ballistic hops | 17.0 | 17.4 | 18.3 | 19.5 | 20.7 | 21.7 | 22.9 | 24.1 | 25.6 | 28.1 |
| 7-iron (`seven_iron`) | 3.6 | 4.1 | 4.6 | 5.3 | 6.2 | 7.4 | 8.8 | 10.2 | 12.5 | 16.4 |
| wedge | 0.2 | 0.4 | 0.6 | 0.9 | 1.3 | 1.9 | 2.8 | 4.0 | 5.9 | 9.4 |

`GREEN`:

| Slope | −3° | −2° | −1° | 0° | +1° | +2° | +2.5° | +3° | +3.1° |
|---|---|---|---|---|---|---|---|---|---|
| driver | 18.2 | 19.5 | 21.0 | 23.2 | 26.7 | 32.2 | 37.4 | 50.2 | 60.1 |
| 7-iron (`seven_iron`) | 4.4 | 5.2 | 6.3 | 8.0 | 10.6 | 15.6 | 20.3 | 33.1 | 43.0 |
| wedge | 0.2 | 0.4 | 0.7 | 1.3 | 2.3 | 4.6 | 7.7 | 18.1 | 27.4 |

`simulate_shot` with `Tour::SEA_LEVEL` (hops through air); the "ballistic
hops" row is the same driver with `settle` given
`Tour::new(0.0, Wind::None)`. Downhill, the hops end by the slope rule
above and the ball rolls the rest of the way; the rule is evaluated on the
hop actually flown, so drag's loss over a hop
counts against the slope's return. On `FAIRWAY` the driver makes 8–10
impacts and the 7-iron 5 or 6 on every slope here; the runs grow with the
downslope mostly through the roll.

Lateral run (yards, + = right) on a plane falling to the right across the
line of play, all `Rest`; the same plane falling left gives exactly the
mirror image:

| Plane | driver | 7-iron (`seven_iron`) | wedge |
|---|---|---|---|
| `FAIRWAY` 2.5° | +3.7 | +1.9 | +1.0 |
| `FAIRWAY` 5° | +9.2 | +5.4 | +3.8 |
| `GREEN` 2.5° | +11.6 | +7.4 | +3.5 |

`GREEN` holds a ball at rest up to 3.1° (`tan θ ≤ 0.055`, 3.15°). Close to
that the net deceleration of a slow roll, `(5/7)·g·cos θ·(μ_eff − tan θ)`,
is tiny and the last meters take a long time: on a 3° downhill green a
drive stops 46 s after landing and a 7-iron 43 s; at 3.1° they take 82 s
and 78 s, and the wedge 71 s. The default 120 s `max_time` covers these.
Beyond the hold angle a ball on `GREEN` never stops, and on a long enough
slope it runs until `max_time`.

**Landing shots on a 1.5% green.** Whitlark & Pringle (2012) fired 857
shots from a pitching machine onto a bermudagrass green sloping 1.5%,
mowed to about 12 ft on the stimpmeter, uphill and downhill, at four
launcher settings chosen to bracket the fairway and rough landings of
their Table 1 (its angles read here as descent angles from horizontal).
Mean bounce and roll was 32.5 ft downhill and 13.5 ft uphill. dimple on `GREEN.with_stimp(12.0)` through sea-level air, for the
landings of their Table 1:

| Landings | uphill (ft) | downhill (ft) | difference | downhill / uphill |
|---|---|---|---|---|
| five fairway landings, mean | 22.3 | 35.7 | 13.4 | 1.60 |
| all ten, mean | 38.3 | 56.0 | 17.7 | 1.46 |
| measured (four settings) | 13.5 | 32.5 | 19.0 | 2.41 |

The fairway 7-iron (55 mph, 46°, 5100 rpm) runs 26.3 ft uphill, 33.4 ft
on level ground and 42.6 ft downhill. The uphill-to-downhill difference
is close to the measured one; the ratio is lower because dimple's uphill totals
are longer than the measured mean. Their four settings are not published,
so the comparison cannot separate the model from the shot mix. Turf firmer
by 0.08 in on the TruFirm ran 9–12 ft further at their uphill sites;
`compliance` 0.9 instead of 1.0 adds 7.8 ft to dimple's uphill fairway mean
(9.1 ft for the 7-iron), so `compliance` is dimple's firmness setting.
`FIRM_GREEN`'s lower `mu_slide` changes only the high-spin landings (the
sand wedge and 9-iron), because the friction limit binds only there.
Green speeds of 11 ft 3 in and 12 ft 3 in made no difference to their
bounce and roll; in dimple they change it by up to 5% (about 4% for the
7-iron). Their sand wedge
landing (49 mph, 51°, 9100 rpm) finishes 0.3–1.8 ft behind its landing
point ("in some cases this distance was negative"), and their shallow,
fast landings ran "over 100 feet on several occasions": dimple's shallow
rough landings run 83–87 ft downhill. Source: Whitlark, B. and Pringle, M.
(2012), "Can't Stop Your Ball on the Green?", *USGA Green Section Record*
50(25); summary values only.

**Break across a slope.** Under a constant rolling resistance a ball
rolled across an incline stops at a displacement angle from its launch line
that does not depend on its launch speed, at a distance proportional to v₀²
(Cross 2024): the motion scales exactly with v₀. dimple reproduces this
under `RollLaw::Constant` at `GREEN`'s stimp: launched 10° off straight
uphill on a 1° or 2° plane, the ball stops 1.125° or 2.035° off its launch
line at every speed from 0.5 to 3 m/s, at 0.832 or 0.694 m per (m/s)². Under
the presets' quadratic law, which resists a faster ball more, the angle
drifts with speed: from 1.36° to 1.12° on the 1° plane and from 2.41° to
2.04° on the 2° plane over 0.5–3 m/s. Source: Cross, R. (2024), "Trajectory
of a golf ball on a sloping green", *Eur. J. Phys.* 45, 035002,
[doi:10.1088/1361-6404/ad2c2d](https://doi.org/10.1088/1361-6404/ad2c2d).

### Surfaces

| Surface | Restitution | `mu_slide` | `compliance` | `mu_roll` | `roll_law` | `pivot` | Stimp |
|---|---|---|---|---|---|---|---|
| `Surface::FAIRWAY` | `Penner2002` | 0.40 | 1.0 | 0.12 | `Quadratic { v_ref: 2.5 }` | 1 mm | 5.2 ft |
| `Surface::GREEN` | `Penner2002` | 0.40 | 1.0 | 0.055 | `Quadratic { v_ref: 2.5 }` | 1 mm | 11.4 ft |
| `Surface::FIRM_GREEN` | `Penner2002` | 0.25 | 1.0 | 0.055 | `Quadratic { v_ref: 2.5 }` | 1 mm | 11.4 ft |

| Parameter | Source |
|---|---|
| eq. 5 restitution | Penner (2002), from Penner (2001), *Am. J. Phys.* 69, 563, measured on relatively soft turf |
| θc law (eq. 8), `compliance = 1.0` | Penner (2002): θc = 15.4° best fit to Haake's (1990) green impact at 18.6 m/s, 44.4° from vertical, 484 rad/s backspin; linear in v and φ as a first-order approximation |
| `mu_slide = 0.40` | Daish (1972), as used by Penner for a ball on a green |
| `FIRM_GREEN` `mu_slide = 0.25` | Penner (2002), Fig. 8: his firm (tournament) green, everything else unchanged |
| Stimpmeter exit speed 1.83 m/s | Penner (2002), "The physics of putting", *Can. J. Phys.* 80, 83, from Holmes (1986), *The Physics Teacher* 24, 401 |
| 5 mm transition height | Penner (2002), set above the ~2 mm equilibrium sinkage of a ball on a green (Hubbard & Alaways) |
| Tilt fade below 2 m/s of normal speed | Chosen, not measured (see [Hops and the bounce-to-roll transition](#hops-and-the-bounce-to-roll-transition)) |
| `mu_roll`, `roll_law` | Presets: dimple's empirical quadratic law, not fitted to published data: `GREEN`'s `mu_roll` gives stimp 11.4 ft; `v_ref = 2.5 m/s` is chosen for plausible full-shot rollouts. `RollLaw::Constant`: Penner (2002), eq. 9. `RollLaw::Exponential`: the form of Burritt, A. and Henrikson, E. (2024), "Modeling golf ball roll with a variable coefficient of friction", *Sports Engineering* 27:39, [doi:10.1007/s12283-024-00481-5](https://doi.org/10.1007/s12283-024-00481-5) (abstract only consulted). Open conference version: *Engineering of Sport 15* (ISEA 2024), Loughborough University repository, [doi:10.17028/rd.lboro.27051061.v1](https://doi.org/10.17028/rd.lboro.27051061.v1), CC BY 4.0, from which the surface values quoted here are read off Figure 1. |
| `pivot` | Hertzian full-slip torsional couple (Johnson 1985) with an assumed 4 mm contact radius; uncalibrated |
| `FAIRWAY` `compliance = 1.0` | Penner's green value, chosen for the fairway (see below) |

Penner uses one turf model throughout: eq. 8 unchanged, with his green
restitution and μ, both for drives compared against empirical drive runs
(Fig. 5) and for iron shots onto greens. He represents a firmer green by a
lower μ (0.25, Fig. 8) rather than by changing θc, noting that eq. 8 would in
general need adjusting for such greens; `FIRM_GREEN` is that green. He also
uses his green's rolling resistance for drives, noting that a fairway's
would be expected to be greater. `FAIRWAY` follows him: it takes his green's
impact model whole, `compliance` 1.0, since a fairway is no firmer than
his green, and differs from `GREEN` only in its rolling resistance, stimp
5.2 ft, which is chosen, not measured. It bounces exactly like `GREEN` and
rolls shorter. The run is sensitive to the compliance. Fairway run (yards
past carry, sea level) by `compliance`, other fields as `Surface::FAIRWAY`:

| `compliance` | driver | 7-iron (`seven_iron`) | wedge |
|---|---|---|---|
| 1.10 | 12.8 | 3.5 | 0.1 |
| 1.00 (`FAIRWAY`) | 17.4 | 5.3 | 0.9 |
| 0.90 | 23.2 | 7.4 | 2.0 |
| 0.85 | 26.1 | 8.4 | 2.7 |
| 0.75 | 32.3 | 10.5 | 4.0 |
| 0.50 | 47.9 | 15.5 | 7.4 |

With `GREEN`, the impact reproduces Penner's published drive results: a
landing at 26.8 m/s, 68.0° from vertical, 224 rad/s backspin rises 2.82 m on
the first bounce, and 26.2 m/s at 46.2° with 315 rad/s rises 1.18 m (Penner
reports 2.82 m and 1.17 m). With ballistic hops, as in Penner's model,
those drives make 16 and 8 impacts before rolling and stop 37.1 m and
17.9 m past the landing point (Penner's Fig. 4: about 37 m and 19 m);
through sea-level air they make 13 and 8 impacts and stop 32.0 m and
17.0 m past it. A 9-iron landing at 23.8 m/s, 34.0° bounces backward
with 809 rad/s backspin, and its backward run grows with backspin at 1000
and 1200 rad/s (Penner's Fig. 7). On `FIRM_GREEN` (μ = 0.25) the same
landings first bounce forward still carrying backspin, then run back
(Penner's Fig. 8).

### Run-out on a fairway

Run (yards past carry) for launches across the bag in Trackman's 25 °C
sea-level air, with the descent angle at landing:

| Launch (mph / launch° / rpm) | Lands at | `FAIRWAY` | `GREEN` |
|---|---|---|---|
| Lob wedge, tour (72 / 32 / 9000) | 46.6° | 1.9 | 2.9 |
| Lob wedge, amateur (65 / 32 / 6500) | 43.9° | 6.0 | 9.3 |
| Sand wedge, tour (80 / 30 / 9500) | 47.3° | 1.5 | 2.2 |
| Sand wedge, amateur (72 / 30 / 7500) | 44.3° | 4.7 | 7.3 |
| Pitching wedge, tour (102 / 24.2 / 9300) | 48.7° | 1.9 | 2.8 |
| Pitching wedge, amateur (85 / 26 / 7500) | 44.4° | 5.3 | 8.1 |
| 50-yard pitch (52 / 30 / 6000) | 37.4° | 7.7 | 11.7 |
| 50-yard pitch, less spin (52 / 32 / 4000) | 39.1° | 11.7 | 16.8 |
| Lob-wedge chip (28 / 36 / 3000) | 38.2° | 7.2 | 11.2 |
| Low chip (30 / 18 / 2500) | 19.0° | 15.7 | 21.4 |
| 7-iron, tour (120 / 16.3 / 7100) | 45.5° | 6.4 | 9.6 |
| 7-iron, amateur (105 / 17 / 5900) | 40.0° | 10.7 | 15.4 |
| Driver, tour (171 / 10.4 / 2545) | 38.8° | 18.4 | 24.2 |

Every launch in the table makes 4–10 impacts before it slides and rolls. Published measurements of run on a
fairway exist for drives:

- The PGA Tour average drive (Trackman's 171 mph, 10.4°, 2545 rpm, landing
  at 39°) lands at 38.8° in dimple and runs 18.4 yd. The USGA/R&A proposed
  bounce model, a linear fit to launcher tests on real fairways, gives
  `d = 79.1 − 1.6·α` yd for landing angle α: 16.7 yd at 39°. ShotLink drives
  of 2016–18 quoted there have a median run of 21 yd and a mean of 22 yd at
  39°, and the R&A/USGA agronomy study measured 17.8 ± 13.3 yd over 624
  fairway shots landing at 32–50° with about 1800 rpm.
- Trackman's tour averages land the 7-iron at 51° and the pitching wedge
  at 52°; dimple lands them at 45.5° and 48.7°. The proposed line reaches
  0 yd at about 49°, while its individual site fits give 0–21 yd at 50°;
  the tour 7-iron's 6.4 yd and the wedges' 1–2 yd lie in that spread.
- These fairway tests used driver-like impacts (26–32 m/s, about
  1800 rpm). A slower, higher-spin wedge or pitch should run less than the
  line at the same landing angle, as it does here (the amateur lob wedge:
  6.0 yd at 43.9°, line 8.9 yd). No published measurement of wedge, pitch
  or chip run on a fairway is known to dimple.

`FAIRWAY` runs 65–76% as far as `GREEN` for every launch, through its
higher rolling resistance alone, since the hops are the same. The rolling
resistance sets much of the run, so `mu_roll` (or `with_stimp`) is the
knob for a faster or slower fairway: the tour drive, amateur lob wedge and
less-spun pitch run 20.9, 7.4 and 13.8 yd at `mu_roll` 0.08 and 15.9, 4.6
and 9.4 yd at 0.25.

Sources: USGA and The R&A (2021), *Proposed Bounce Model for Use in
Evaluating Optimum Overall Distance*
([PDF](https://assets-us-01.kc-usercontent.com/c42c7bf4-dca7-00ea-4f2e-373223f80f76/250b3e34-9172-4ad1-991f-9d5560cbf832/Proposed%20Bounce%20Model%20for%20Use%20in%20Evaluating%20Optimum%20Overall%20Distance.pdf)),
which also reports the ShotLink runs; R&A and USGA (2020), Distance
Insights Project, *R03: Agronomic Impacts on Bounce and Roll Distance*
([PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/cdd48ad2-2ed1-4e67-ad53-1a6522110779/R03%20-%20Agronomic%20Impacts%20on%20Bounce%20and%20Roll%20Distance.pdf));
Trackman PGA Tour averages
([trackman.com](https://www.trackman.com/blog/golf/introducing-updated-tour-averages)).

### Validation against published bounce data

Biber et al. (2023) filmed 1,023 bounces of a modern multi-piece ball fired
from a launcher onto two surfaces, neither of them a green: an artificial
teeing mat 30 mm thick on a rigid wooden base (Campaign A, 330 bounces)
and a natural teeing area (Campaign B, 693 bounces). They publish the
measurements as figures and as fits: rigid-bounce and Penner-model
parameters (their Table 3), and a piecewise-affine map from landing to
lift-off for each surface, with one branch for bounces that enter rolling
and one for bounces that slip throughout (their eqs (8)–(9), fitting the
data to 6.6% and 7.8%, Table 2). dimple is compared against those fits, not
the raw measurements: `impact` and each affine fit on an equally weighted
grid of landings inside the campaign's measured range (natural tee: 120
landings, 5–30 m/s, 20–80° descent, −2000 to 12,000 rpm; mat: 96 landings,
10–40 m/s, 25–70°, −2000 to 8000 rpm), with the paper's error measures in
its scaled units (velocities in ball radii per second, spin in rad/s): δ
the mean of |p_F − P_F| / |p₀| over the landings, and Δ the same per
component. Through dimple, the paper's own Penner-model parameters for the
tee (r 0.222, μ 0.999, k_P 2.973e-4) score 23.5% against the fit, against
the 21.3% the paper reports against the data, so the fit stands in for the
data to within a few points.

Natural tee (Campaign B):

| Against the tee fit | δ | Δẋ | Δẏ | Δω |
|---|---|---|---|---|
| `GREEN` and `FAIRWAY` | 12.2% | 7.4% | 3.0% | 7.7% |
| `GREEN` on a rigid plane (`compliance` 0) | 20.6% | 14.8% | 6.7% | 10.1% |
| `GREEN` with `mu_slide` 1.0 | 22.2% | 11.5% | 3.2% | 17.2% |

For comparison, against the data the paper's rigid-bounce fit scores
35.0%, its Penner-model fit 21.3% and its fixed-angle fit 19.2%.
`GREEN`'s mean normal rebound ratio |ẏ_F / ẏ₀| is 0.388 against the fit's
0.408, and of the 80 backspin landings 26 leave with topspin against the
fit's 23. Most of the remaining error is tangential: dimple keeps about a
third of the forward speed the fit keeps (mean ẋ_F / ẋ₀ 0.107 against
0.313). Balls that enter rolling leave the ground slipping with excess
topspin, contact-point velocity ẋ + ω < 0 (the paper's Fig. 6); dimple's
tilted contact plane gives that sign in the ground frame for 56 of the 57
landings on the fit's rolling branch, with about 70% of the fitted mean
(−41 against −58 radii/s), but not the spread of the measurements, which the authors put down to the natural
variability of turf; they suggest that the tangential and spin error left by
rigid-contact models needs a tangential compliance of the ground, which dimple
does not model. The fits'
μ ≈ 1 is not adopted: it grips every landing, while the slipping bounces
leave with their backspin, and in dimple it nearly doubles the error.
Their typical driver landing (28.5 m/s, 37.3°, 34.9 rev/s) rebounds from
`GREEN` with a normal ratio of 0.354 against the tee fit's 0.342 and a
forward ratio of 0.214 against 0.277, leaving with topspin as in the fit.

In the paper's units Penner's θc law is β = k_P·s·φ, with s the speed in
ball radii per second and φ in radians, and `compliance` 1.0 is
k_P = 3.98e-4. The paper's fitted k_P are 1.543e-4 for the mat and
2.973e-4 for the tee, i.e. `compliance` 0.39 and 0.75. On the tee
`compliance` 0.75 or 0.8 scores 11.7% or 11.6% against the presets' 12.2%,
a difference well inside the fit's own error; `GREEN` and `FAIRWAY` keep
Penner's 1.0, with which `GREEN` reproduces his published results.

Artificial mat (Campaign A): `Penner2002` restitution is far too low for a
mat on a rigid base. `GREEN`'s mean normal rebound ratio is 0.32 against
the fit's 0.58 (the paper fits r = 0.42–0.54 for the mat, Table 3), and
`GREEN` scores δ 30.5% with Δẏ 18.1%; no change of `compliance` or
`mu_slide` alone brings Δẏ down. A constant restitution does: the mat
example under [Other surfaces](#other-surfaces) scores δ 17.2% with Δẏ
4.6%. Use `Restitution::Constant` for mats.

Source: Biber, S. W., Jones, K. M., Champneys, A. R., Green, R. and
Szalai, R. (2023), "Measurements and linearized models for golf ball bounce
on a green", *Sports Engineering* 26:50,
[doi:10.1007/s12283-023-00442-4](https://doi.org/10.1007/s12283-023-00442-4),
CC BY 4.0. The fitted coefficients quoted here are from its Table 3 and
eqs (8)–(9); dimple's figures are its own runs against those fits.

### Stimpmeter

Green speed is quoted as a stimpmeter reading: how far, in feet, a ball
released down the USGA stimpmeter rolls on level green. The ball leaves the
ramp already rolling at `STIMPMETER_SPEED` = 1.83 m/s (6.0 ft/s), the value
Penner (2002, "The physics of putting") takes from Holmes's (1986) analysis
of the device. As a cross-check from the geometry (release notch 30 in up
the ramp, released at about 20°), a solid ball rolling 30 in down a 20°
plane without loss would leave at about 1.91 m/s, and about 1.88 m/s in
the 145° V-groove, whose rails touch the ball off its lowest point; dimple
cites no direct measurement of the exit speed.

With the rolling deceleration `μ_eff·g/(1+k)` above, a ball leaving at
`v₀` rolls `d = (1+k)·I / (2·mu_roll·g)`, with `I = 2∫₀^v₀ v/f(v) dv`, so a
surface reads stimp `d` when `mu_roll = (1+k)·I / (2·g·d)`. With
s₀ = v₀/v_ref, `I` has a closed form for each law:

```
Constant:     I = v₀²
Quadratic:    I = v_ref² · ln(1 + s₀²)
Exponential:  I = 2·v_ref² · (1 − e^(−s₀)·(1 + s₀))
```

Both parameterized forms tend to `v₀²` as `v_ref` → ±∞. `Exponential` works
for a negative `v_ref` unchanged.

`Surface::with_stimp(feet)` sets `mu_roll` by this closed form, keeping
every other field (including `roll_law`), and `Surface::stimp_feet()` is
the inverse. Both use k = 2/5, the uniform solid sphere of
`Tour`. `with_stimp` returns the surface unchanged for a
reading that is not finite and positive, or a law whose `v_ref` is zero or
NaN, for which `stimp_feet` is NaN.

| Stimp (ft) | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|
| `GREEN.with_stimp(…).mu_roll` | 0.105 | 0.0785 | 0.0628 | 0.0523 | 0.0449 |

A simulated stimpmeter roll (rolling at 1.83 m/s on level ground) on each
of these surfaces stops within 1e-9 (relative) of the requested distance.

**On a slope.** Brede (1991) showed that the readings a up and b down a
uniform slope combine as 2ab/(a+b), which reads the level value on slopes
from 0 to 5.6% (a laboratory carpet runway) and on field greens up to 6%;
the plain average of the two does not. Constant rolling resistance
decelerates the ball at g(μ cos θ ± sin θ)/(1+k), so under
`RollLaw::Constant` the correction is exact apart from the 1/cos θ of
measuring along the surface (0.005% at 1%, 0.18% at 6%). Under the presets'
quadratic law the up- and downhill rolls cover different speeds and it is
approximate:

| Grade | 1% | 2% | 3% | 4% | 5% | 6% |
|---|---|---|---|---|---|---|
| `GREEN` up / down (ft) | 9.95 / 13.41 | 8.82 / 16.26 | 7.93 / 20.73 | 7.20 / 28.94 | 6.59 / 51.44 | neither rests |
| `GREEN` 2ab/(a+b) vs level 11.42 ft | +0.04% | +0.17% | +0.43% | +0.95% | +2.3% | – |
| `FAIRWAY` 2ab/(a+b) vs level 5.23 ft | +0.01% | +0.05% | +0.11% | +0.20% | +0.33% | +0.49% |

On `GREEN` the plain average is already 10% long at 2%. `GREEN` holds a
ball at rest only to a 5.5% grade (`tan θ ≤ mu_roll`), so a downhill
stimpmeter roll on a steeper `GREEN` never stops, and an uphill one stops
once and rolls back; with `RollLaw::Constant` at the same stimp
(`mu_roll` 0.0687) it holds to 6.9%.

Source: Brede, A. D. (1991), *Agron. J.* 83:425–426,
[doi:10.2134/agronj1991.00021962008300020032x](https://doi.org/10.2134/agronj1991.00021962008300020032x).

### Other surfaces

Beyond the presets above, dimple adds a surface preset only when its
parameters come from a published source. Rough, bunker sand and firm or
soft fairways have no published set of this model's parameters
(restitution, impact friction, contact angle and rolling resistance
together), so there are no presets for them. Build them from a preset with struct update
syntax and, for the roll, a stimp reading:

```rust
// A soft, slow fairway: rolling like stimp 4 ft, rebounding a little softer.
let soft_fairway = Surface { compliance: 1.1, ..Surface::FAIRWAY }.with_stimp(4.0);
// A firm, fast fairway.
let firm_fairway = Surface { compliance: 0.75, ..Surface::FAIRWAY }.with_stimp(7.0);
// A soft fairway whose rolling resistance grows exponentially with speed.
let exp_fairway = Surface {
    compliance: 1.1,
    roll_law: RollLaw::Exponential { v_ref: 3.8 },
    ..Surface::FAIRWAY
}
.with_stimp(4.0);
```

An artificial mat on a rigid base rebounds much more than turf, which the
`Penner2002` restitution does not give (see
[Validation against published bounce data](#validation-against-published-bounce-data)).
This example scores δ 17.2% against Biber et al.'s mat fit (`GREEN`: 30.5%).
It is built from one mat in one study and is not a preset:

```rust
// A hitting mat on a rigid base.
let mat = Surface {
    restitution: Restitution::Constant(0.45),
    compliance: 0.4,
    mu_slide: 0.25,
    ..Surface::GREEN
};
```

`compliance` moves the bounce and the run most (see the compliance table
above); `mu_slide` sets how much backspin an impact checks; `with_stimp`
sets the roll at stimpmeter speed and `roll_law` how it changes with
speed.

### Plugging in your own ground

The ground API has three layers:

1. **`Surface`** — flat ground of one material, also a complete ground
   model: `simulate_shot(&input, &Tour::SEA_LEVEL, &Surface::FAIRWAY)`, or with
   a field changed, `&Surface { mu_roll: 0.09, ..Surface::FAIRWAY }`.
2. **`GroundPhysics<T: Terrain>`** — the physics above over any terrain.
   `Terrain::sample(x, z)` returns a `GroundSample { height, normal,
   surface }`; `FnTerrain(|x, z| ...)` wraps a closure, so per-point
   surfaces need no new types:

   ```rust
   let range = GroundPhysics::new(FnTerrain(|x: f64, z: f64| GroundSample {
       height: 0.0,
       normal: Vec3::Y,
       surface: if on_green(x, z) { Surface::GREEN } else { Surface::FAIRWAY },
   }));
   simulate_shot(&input, &Tour::SEA_LEVEL, &range);
   ```

   Its public fields are the tunables: `dt` (1 ms), `min_hop_height`
   (5 mm), `stop_speed` and `roll_threshold` (0.01 m/s), `max_impacts`
   (200), `max_time` (120 s).
3. **`GroundModel`** — replace the physics: implement `height(x, z)` and
   `settle(landing, model)` (`simulate_shot` passes the shot's
   `AeroModel`: its `acceleration` for the hops, `radius_m` and
   `inertia_factor` for the contacts, and `gravity()` throughout; an
   implementation uses as much of it as it models), and
   optionally `landing(a, b)` (by default the
   solved crossing of the straight segment with the heightfield; an engine
   can substitute its own sweep).
   It is object-safe, so `&dyn GroundModel` works. The building blocks
   `ground::impact` and `GroundPhysics::step` are public (both take the
   model), and the result
   types (`GroundResult`, `Contact`, `GroundPoint`) have `new`
   constructors. `BallState` (time, position, velocity, spin) is a plain
   struct.

`simulate_shot` launches from `(0, height(0, 0), 0)`, ends the flight at
`landing()`, and hands the landing state to `settle()`. `simulate_flight`
always lands on the plane y = 0. Both place the landing point exactly on
the ground at the solved crossing.

**Obstacles.** `ground::impact` takes any outward normal, not only ground
with `normal.y > 0` (that restriction belongs to `Terrain` and
`GroundPhysics`, which need a heightfield): a wall, a tree trunk or the
underside of an overhang. The spin acts on these bounces as on the ground,
through the friction impulse: backspin checks a ball off a cart path and
sidespin kicks it sideways off a wall. For a hard surface use
`compliance: 0.0` (a rigid flat contact) and a `Restitution::Constant`;
dimple has no published values for cart paths, walls or bark. A
`GroundModel` whose `landing()` also reports where a flight segment meets
an obstacle can bounce the ball in `settle()` with `impact` and the
obstacle's normal, then fly on with `simulate_shot_from(&rebound, model,
self)`; the crate's `tests/obstacles.rs` does this for a tree trunk. The
pattern keeps only the rebound point of the flight between the obstacle and
the next landing: `GroundResult` has no slot for flight samples, so the
points of the continued flight are not in the shot's result.

### Limitations

- **Natural-turf data only.** The restitution fit and the θc law come from
  Penner's green and are checked against a natural tee (Biber et al. 2023).
  An artificial mat on a rigid base needs `Restitution::Constant` (see
  [Other surfaces](#other-surfaces)). `FAIRWAY` differs from `GREEN` only
  by rolling resistance; neither its compliance (Penner's green value) nor
  its stimp has a published fairway source, and the run depends strongly
  on both.
- **One impact measurement.** θc rests on a single Haake impact and is
  extrapolated linearly in speed and incidence angle. Spin and turf
  firmness are not inputs. In Biber et al.'s units Penner's green
  (`compliance` 1.0) has a higher k_P than either of their fits (the mat
  0.39, the tee 0.75).
- **The tilt fade is a modeling choice.** Fading θc in below 2 m/s of
  normal speed ends the near-grazing hop train that the linear law
  sustains; no measurement of soft, grazing turf contacts sets its onset or
  shape. Runs change little for onsets of 1.5–3 m/s. A bump-and-run
  launched under about 12° lands inside the fade and skims 10–16 impacts
  before it rolls.
- **Tangential compliance is not modeled.** Measured bounces leave the ground slipping, with more topspin than rolling,
  and the authors suggest that the tangential and spin error of rigid-contact
  models needs a tangential compliance of the ground (Biber et al. 2023).
  dimple's tilted rigid contact plane gives the sign for balls that grip, with
  about 70% of the fitted mean on natural turf, but not the spread of the
  measurements, and keeps about a third of the forward speed the tee fit keeps.
- **Lateral ground effects are unvalidated.** The kick of a curving ball,
  the sidespin throw at impact and rifle-spin drift follow from the
  rigid-sphere contact model; no published measurement is compared against.
  The pivot lever is uncalibrated.
- **No rolling-resistance law is validated above 3 m/s.** Only the
  low-speed `mu_roll` is tied to a measurement (the stimpmeter). The
  presets' quadratic form and `v_ref` are empirical, and the only published
  speed-dependent fits (Burritt & Henrikson 2024) cover putts, 0–3 m/s.
  `RollLaw::Exponential` takes such a per-surface fit; extrapolated past
  3 m/s it grows quickly (`v_ref` = 1.6 m/s gives about 150 × `mu_roll` at
  8 m/s, so a fast roll stops almost at once), and with a negative `v_ref`
  a ball has no terminal speed on a slope steeper than its hold angle.
- **The quadratic law exceeds sliding friction at speed.** On `FAIRWAY`
  μ_eff passes `mu_slide` (0.40) at 3.8 m/s and reaches 1.35 at 8 m/s,
  more than the turf could exert on a sliding ball. Drives, pitches and
  chips reach the ground at up to about 4 m/s, so this touches the start of
  their roll, and fast rolls started by the caller and long downhill rolls;
  `RollLaw` offers gentler laws.
- **Spin is not decayed during hops.** The flight's decay rate is
  calibrated at flight speeds and overstates the decay at hop speeds (see
  [Hops and the bounce-to-roll transition](#hops-and-the-bounce-to-roll-transition)).
- **Heightfield only.** One height per (x, z): no walls, lips or cup. A
  step in the height is crossed as if the ball were lifted onto it.
  Crossings are found at 8 sub-steps per 1 ms flight or hop step, so a
  ridge or wall narrower than a sub-step (about 1 cm at full speed) can be
  missed.
- **Curved ground is followed, not flown over.** Each ground step is taken
  in the tangent plane and the ball is put back on the surface: no launch
  off crests, and the normal component of the velocity is removed at
  concave transitions. Exact on planes; a first-order approximation on
  curved ground.
- **Slopes: rolling checked, landings beta.** Rolling on a slope (normal
  load, downhill pull, hold where `tan θ ≤ mu_roll`) is checked against
  Brede's stimpmeter slope correction and Cross's break-angle invariance:
  exact under `RollLaw::Constant`, approximate under the presets' quadratic
  law. The impact model on a slope is checked only against the 1.5%-green
  means of Whitlark & Pringle (2012), whose shot mix is unpublished; dimple
  knows of no published bounce-and-spin measurements on steeper slopes.
  Penner's impact model was fitted on level turf, and he attributes its
  slightly long drive runs to its assumption of level ground.
- **The slope hop rule is a modeling choice.** It ends hop trains that
  Penner's impact model would sustain down a slope (see
  [Hops and the bounce-to-roll transition](#hops-and-the-bounce-to-roll-transition)).
  Its margin, one 5 mm hop's energy per impact, carries over Penner's
  transition height; it is not measured.
- **Slopes beyond the hold angle never come to rest.** Steeper than
  `atan(mu_roll)` (3.15° on `GREEN`) the ball rolls on until the ground
  changes or `max_time` ends it with `MaxTime`. Just inside the hold angle
  the last meters of a slow roll take over a minute (82 s for a drive on a
  3.1° `GREEN`; see [Terrain and slopes](#terrain-and-slopes)); a
  `max_time` set below the 120 s default can cut such rolls short.

## Implementation

The aerodynamic model lives in `src/aero.rs`, the air in `src/air.rs`, the
ground model in `src/ground.rs`. Key types:

- **`AeroModel`** — the air's effect on one ball: acceleration without
  gravity, spin over a flight step, radius, inertia factor, and gravity
  (standard by default).
- **`Tour`** — dimple's model: the tour ball in an air density and a wind
  (`Tour::new`, `Tour::SEA_LEVEL`).
- **`standard_air_density`** — density from altitude and temperature.
- **`Wind`** — none, uniform, height-dependent (`PowerLaw`) or your own
  (`Wind::custom`, any `WindField`).
- **`Surface`**, **`Terrain`**, **`GroundPhysics`**, **`GroundModel`** — the
  ground (above).

Flight simulation: `simulate_flight(input, model) -> FlightResult`

RK4 integration at 1ms timestep, gravity plus the model's acceleration at
every stage, the spin advanced by the model once per step, trajectory
sampled every 10ms. Stops at the solved crossing of ground level.

Full shot: `simulate_shot(input, model, ground) -> ShotResult` runs the
flight from the tee on `ground` until it crosses the terrain, then
`ground.settle(landing, model)`: impacts, hops through the same air
and wind, slide and roll to rest (`ShotResult.ground`, with absolute
times).

From any state: `simulate_flight_from(start, model)` and
`simulate_shot_from(start, model, ground)` fly the same way from a
`BallState` (time, position, velocity, spin), for a ball that flies on after
hitting a tree, a wall or a flagstick. Time continues from `start.time`,
positions stay in the frame above, carry, lateral offset and total are
measured from the start's horizontal position and apex from its height, and
a start below the ground lands at once, placed on it. `simulate_flight` and
`simulate_shot` are these calls from the launch state.

## Test Suite

Unit tests live beside the code in each module. Test categories:

- **Structural invariants** — energy conservation, origin start, axis normalization
- **Sign conventions** — azimuth, slice/hook, push-draw, mirror symmetry, rest position
- **Cd/Cl** — no lift without spin, lift and drag rise with spin, saturating spin drag, coefficients independent of speed at a given spin ratio, driver and 7-iron values, Li et al. (2017)
- **Low-spin lift term** — zero and flat at zero spin, no step anywhere, peak height and position, lift rising with spin ratio through it, fading out above a spin ratio of 0.045
- **Descent angle** — iron DA within expected ranges; the reported descent angle is the landing velocity's
- **Altitude/temperature** — ISA density at 1,524 m, Denver gain, hot-cold delta, directional checks, the preset driver's 5,000 ft gain
- **Spin decay** — exponential decay rate against Trackman measurements, about a fixed axis
- **Wind profile** — the power law and its 0.1 m/s threshold, and the hops sampling the profile
- **Wind** — a headwind shortens carry and costs more than the same tailwind gains, a crosswind drifts the ball downwind, the built-in winds ignore a vertical component, the wind carries through the hops, strong wind does not keep the ball bouncing, and a custom ground model receives the shot's model, wind included
- **Custom models and winds** — a drag-only model flies and lands through `&dyn AeroModel`; a model with the default gravity, or standard gravity given explicitly, flies bit for bit as `Tour`; lower gravity flies and rolls farther; a constant `Wind::custom` equals `Wind::Uniform` bit for bit; a gust that starts mid-flight shortens the carry less than the same wind throughout; a custom updraft raises the apex and a downdraft lowers it
- **Flights from any state** — `simulate_flight`/`simulate_shot` equal the `_from` calls from the launch state bit for bit (flat ground, a slope, wind); a flight split at any sample and resumed lands within 1e-6 m of the whole one (still air and a time-dependent wind), and a split shot comes to rest at the same point; a start below the ground lands at once on it, flat and on terrain; a start on the ground moving into it lands within the first step
- **Obstacles** — `impact` with any outward normal: backspin checks a ball off a rigid cart path, sidespin kicks it sideways off a wall and the kick flips with the spin's sign, a spinless bounce returns the normal velocity times the restitution for walls, a trunk, an overhang and a cart path; a custom `GroundModel` with a tree trunk bounces a drive off the bark and continues with `simulate_shot_from`
- **Launch spin axis** — a pushed shot with no sidespin stays on its start line; exact landing-spin projections
- **Ground: planar regression** — with ballistic hops, the 3D impact reproduces the planar Penner solution with the tilt fade impact by impact (preset driver, 7-iron and wedge on both presets) to float tolerance
- **Ground: physics** — Penner (2002) restitution, drive bounce heights and runs (Fig. 4, ballistic hops), 9-iron spin-back (Figs. 7–8), 5 mm bounce-to-roll transition, compliance tilt (2β rotation, cap), the tilt fade below 2 m/s of normal speed and the low skimming hops it shrinks, `FAIRWAY` bouncing like `GREEN` and rolling shorter from the same landings, grip-regime kick formula, no kick from spin about the normal, tilt-throw magnitude (grip and slip), rifle-spin drift, grazing impacts, pivot friction about the true normal and its independence from the level-ground path, general inertia factor (k = 2/3), exact slide-to-roll time, stimpmeter closed form, `with_stimp`/`stimp_feet` round trip under every roll law and simulated stimpmeter rolls (6–14 ft) to 1e-9, preset stimp readings, `FIRM_GREEN` as Penner's Fig. 8 green
- **Ground: roll laws** — every law is `mu_roll` at zero speed, exponential roll-outs against the closed form (positive and negative `v_ref`), Penner's ρ_g = 0.131 reads 6.0 ft under the constant law, an infinite `v_ref` is the constant law, a rising exponential stops a fast roll sooner than the constant law at the same stimp, zero or NaN `v_ref` rejected by the stimp helpers and rolled as the constant law
- **Ground: hops through air** — at zero density hops are bit-for-bit ballistic (on the parabola, `Surface` and `GroundPhysics` identical); through air a hop arrives earlier, nearer and slower with its spin unchanged; the run lost to drag shrinks with altitude; topspin pushes a hop down and backspin holds it up; sidespin curves a hop right (fade spin) or left, mirrored
- **Ground: invariants** — energy never increases at an impact or along the ground phase, with ballistic hops and through air; mirror (±sidespin, ±azimuth) and rotation about Y give mirrored and rotated results; fades rest right, draws rest left
- **Ground: terrain** — rest on slopes (hold below atan `mu_roll`, roll downhill above, both presets; one stop and no chatter rolling uphill), long downhill reports `MaxTime`, near-critical (3.1°) `GREEN` comes to rest within the default `max_time`, slide-to-roll on a plane against the closed form, slip turning across a slope, vertical spin kicks on slopes only, uphill landings check up and downhill landings release, elevated green (shorter carry, lands on it), tee on the terrain, landing and hop crossings on curved terrain to 1e-9 m (ballistic hops on the parabola; hops through air on the surface), first crossing of thin walls and ridges in a landing segment and a hop step, ball kept on the surface; slope hop rule, downhill `FAIRWAY` up to 6° (through air and ballistic, the driver run growing with the slope either way, shorter through air) and `GREEN` up to 2.5° end at rest, side slopes rest and mirror exactly, rest judged on the surface where the ball stops, slide-to-roll on curved ground without flicker; mirror and rotation about Y with a sloped plane through air, and energy on slopes at both densities
- **Ground: validation against published data** — Brede's (1991) slope correction exact under the constant law and within 0.5–1.5% to 4% under the presets' law, `GREEN`'s 5.5% hold grade; Penner's θc in Biber et al.'s (2023) units and their fitted k_P as `compliance`; agreement with their natural-tee fit (δ, Δẏ, normal rebound, spin reversal, lift-off slip sign; a rigid plane and μ = 1 do worse); their mat needs a higher constant restitution; their typical driver landing; Whitlark & Pringle's (2012) 1.5% green (uphill/downhill asymmetry, the sand wedge stopping by its pitch mark, green speed, `compliance` as firmness); Cross's (2024) speed-independent break angle under the constant law
- **Run-out** — the hop train dies out (at most 10 impacts, each hop after the first at most 0.7× the one before, at most two below 1 cm) from wedges to drives on both presets; `FAIRWAY` roll and total for a lob wedge, a pitch and the PGA Tour average drive; the preset drive's carry, run and total; the drive runs shorter on `FAIRWAY` than on `GREEN`
- **Ground: API** — per-point surfaces (fairway-to-green roll), `MaxTime`/`MaxImpacts` termination, tunables, interpolated landing, `&dyn GroundModel`

Run: `cargo test`
