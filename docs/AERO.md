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

A model that wants the viscosity implements `AeroModel` with a drag that
depends on the Reynolds number. `examples/viscosity.rs` does this with a
smooth drag crisis (illustrative constants, not fitted) and the viscosity
from Sutherland's law at the air temperature, and compares carries at 0 °C
and 35 °C in air of the same density.

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
   uncalibrated; landing shots on slopes are checked only against one
   1.5% green (see [Ground limitations](#limitations)).

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

**Status.** Every impact constant of the presets comes from measured turf.
The contact angle is fitted to 693 measured bounces on natural tee turf at
2–39 m/s and landing angles of 16–90° (Biber et al. 2023), and the
presets' `compliance` and `mu_slide` are the medians of measured greens and
fairways (R&A 2008; USGA 2005); `GREEN` is checked against Cross's (2025)
printed fits to 184 bounces on four greens and `FAIRWAY` against Cross's
(2021) fairway-like grass (see [Impact](#impact) and [Surfaces](#surfaces);
every dataset's score is tabulated in [FIT.md](FIT.md)).
Rolling resistance is set by Stimpmeter readings and its rise with speed
(`plow`) by measured roll-outs. Roll-out is within 3 yd of the USGA's 2019
range test at 20–36° (rms 2 yd) and about 11 yd long at 12–18°, a
tour-average drive runs about
23 yd on `FAIRWAY` against 21.7 yd on the PGA TOUR's mean line at 39°
(see [Run-out on a fairway](#run-out-on-a-fairway)), and rolling on slopes
is checked against Brede's (1991) stimpmeter slope correction and Cross's
(2024) break-angle invariance (see [Stimpmeter](#stimpmeter) and
[Terrain and slopes](#terrain-and-slopes)). What no measurement sets: the
impact below a 16° landing, rough, sand, the pivot lever and the lateral
ground effects (see [Limitations](#limitations)).

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
  lands heading 12.4° off-line and leaves its first bounce at 13.7°.
- **Spin along the direction of travel** ("rifle spin") slips the contact
  point sideways, so it drifts the ball laterally while it slides.
- **The tilt adds a small throw the other way.** Each impact is against
  a tilted effective plane (below), so a fraction sin β of ω_y lies in the
  contact plane and slips the contact point sideways at impact: a ball
  spinning clockwise from above (fade sidespin) is thrown left, like side
  spin off a cushion. Pivot friction then clears ω_y, so the throw acts
  mainly at the landing impact. For the ±600 rpm driver on `FAIRWAY` at
  sea level, measured as lateral run across the landing heading, the
  heading-mismatch kick alone (the same landing with ω_y removed) adds
  about 0.95 m, the throw takes back about 0.1 m, and the net is about
  +0.85 m (0.9 yd) on a 20 m run. The fade comes to rest 32.9 m (36.0 yd)
  right of the target line, 5.1 m right of its landing point.

So dimple does not convert sidespin into a lateral nudge; it carries the
spin vector through the contact physics, and these effects follow from it.

### Impact

Each impact follows the construction of Penner (2002), "The run of a golf
ball", *Can. J. Phys.* 80, 931–940, in vector form. The ball deforms the
turf, so the impact is treated as one off a rigid plane inclined by an
effective contact angle β toward the incoming ball. With `v_n = −v·n` the
approach speed, `v_t` the tangential velocity (direction `t̂`, speed
`|v_t|`), `v = |v|`, `φ = atan(|v_t| / v_n)` the incidence angle from the
normal and `θ = 90° − φ` the landing angle from the surface:

```
tan β₀ = (|v_t| / 50 m/s) · (1 − exp(−(θ / 20°)³))   (contact angle)
β    = min(compliance · β₀, φ)
n'   = n cos β − t̂ sin β                          (effective normal)
e    = restitution at |v·n'|                       (Penner eq. 5)
j_n  = (1 + e) · |v·n'|                            (normal impulse / m)
u    = (v − (v·n')n') + ω × (−R n')               (contact slip)
j_t  = min(|u| · k/(1+k), μ · j_n)                 (along −û)
v   += j_n n' − j_t û
ω   += (−R n') × (−j_t û) / (k R²)
Δω_n = −sign(ω_n) · min(|ω_n|, pivot · j_n / (k R²))   about n
kept = (1 − absorption)^(sin θ)                    (sin θ = v_n / v)
v'   = kept · v
ω'   = kept · ω
```

`k` is the ball's inertia factor (`I = k m R²`, `AeroModel::inertia_factor`;
2/5, a uniform solid sphere, for `Tour`) and R its radius (`radius_m`).
Friction either brings the contact point to rest (grip: the ball leaves
rolling on the effective plane) or slides throughout at the Coulomb limit.
No impulse is applied when the ball is not
moving into the ground (grazing).

The last three lines take the share of the rebound a bed such as loose
sand absorbs: `absorption` for a ball landing straight down, less the
shallower it lands (θ is the landing angle from the surface), and nothing
for a grazing touch. Every preset has `absorption` 0, so `kept` is exactly 1
and the rebound is the turf's; a value above 1 counts as 1, which stops the
ball on its mark (see [Other surfaces](#other-surfaces)).

**The contact angle** is fitted to 693 measured bounces on natural tee turf
at 2–39 m/s and landing angles of 16–90° (Biber et al. 2023; the 95%
intervals are 49.2–50.9 m/s for the speed scale and 19.8–21.8° for the
angle scale). Measured from each bounce's impulse, the tilt follows the
tangential speed alone: on steep landings `tan β₀ = |v_t| / 50 m/s`, and
the normal speed adds nothing. It starts to fall below a landing of about
26° and nothing lands below 16°. The turf resists the ball along its path
with a force that is a fraction `|v_t| / 50 m/s` of the normal force; a
ball that lands flat presses little turf ahead of itself, which the factor
in θ gives: 0.97 at 30°, 0.63 at 20°, 0.03 at 6°, and zero as `v_n` goes
to zero. No impact is measured below 16°. Roll-out after landings at
12–36° is compared with the USGA's 2019 range test: within 3 yd at 20–36°
(rms 2 yd), about 11 yd long at 12–18° (see
[Run-out on a fairway](#run-out-on-a-fairway)). The power 3 is the least
certain constant: the bounces allow about 2–4 and the 2019 roll-out line
about 2–3.5, and 3 is the lowest round value in that band at which no
grazing landing leaves faster along the normal than it arrived.

Rebound is 0.07–0.15 livelier than measured for landings flatter than
about 35° (Biber et al. 2023; Cross 2021): the tilted-plane rebound keeps
too much normal speed there. Applying the restitution to the true normal
speed instead was tested and fits the bounces worse.

The same law fits greens. Cross (2025) prints three fits to 184 measured
bounces on four greens at 17–27 m/s and 35–55°: the rebound ratio
(eq. 14), the ratio of horizontal to vertical impulse (eq. 15) and the
offset of the normal force ahead of the ball's center (eq. 16).
`GREEN` reproduces all three: model minus fit is −0.012 / 0.030 for the
rebound ratio, −0.012 / 0.041 for the impulse ratio and −0.51 / 0.64 mm for
the offset (mean / rms). The offset is the level-ground form of the tilt,
`D / R = (sin(β + ϕ) − sin ϕ) / cos(β + ϕ)` with `ϕ = atan mu_slide`: 7.9 mm
mean against 8.2 mm, rising 0.36 mm per m/s against 0.30 and falling
0.14 mm per degree against 0.14. Fitting `compliance` and `mu_slide` to
those three equations alone gives 1.12 and 0.23, against `GREEN`'s 1.10 and
0.22. Cross's coefficient of friction, 0.6–0.8 on the level surface, is the
same impulse as `mu_slide` 0.22 on the tilted plane: for a bounce that
slides throughout the two are tied by `COF = tan(β + atan mu_slide)`. On
Cross's (2021) mown, fairway-like grass, 51 throws at 12.5 m/s and
31–85°, `FAIRWAY`'s rebound ratio is within 0.01 of the measured one from
40° up and 0.07 high at 30–40°.

For a ball with no spin about the vertical traveling in a vertical plane,
this is exactly Penner's planar construction: the normal component
`v_n' = v_n cos β + v_t sin β` (eq. 6), Daish's rigid-surface rebound in
the effective frame (eqs. 2–4: grip `(5/7)·v_t' − (2/7)·R·ω_back`, or slip
with friction impulse `μ(1+e)·v_n'`), rotated back (eq. 7). The test suite
checks the vector form against an independent planar solution, impact by
impact.

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
angle, where the ball meets the effective plane head-on; at the presets'
compliance this binds only on landings faster than about 45 m/s along the
normal. Backspin can
send the ball backward. A rebound that points into the ground stays on the
ground and rolls.

**Pivot friction.** The friction moment about the surface normal is
`pivot` × normal force. It acts about the true normal `n`, not the tilted
`n'`, so it removes only spin about the vertical and never the spin along
the heading that drives a curving ball's kick. The presets use 1 mm: the
full-slip torsional couple of a Hertzian contact, (3π/16)·μ·a (Johnson,
*Contact Mechanics*, 1985), with a contact radius a ≈ 4 mm. At the presets'
sliding friction (0.17–0.22) that formula gives 0.4–0.5 mm; the presets'
1 mm is its value for μ ≈ 0.4, the upper end of an estimate that is
uncalibrated.

Pivot friction never changes the velocity of the impact it acts in, nor a
slide or roll on level ground. It reaches the path only through spin about
the normal carried into a later impact, where the tilted plane turns part of
it into slip, and through the Magnus force that spin gives the hops in
between. With the 1 mm lever a full-shot landing removes all of it in one
impact: a driver landing on `FAIRWAY` can remove up to about 130 rad/s
(~1300 rpm), and a ±600 rpm-sidespin drive lands with 47 rad/s (450 rpm).
So on level ground the path is the same for every lever from 1 mm up. With
no pivot friction at all the spin survives to the later impacts, whose
throws move that drive's rest by about 13 cm; on rigid ground (`FAIRWAY`
with no tilt) it survives only into the hops, which it curves about 45 cm
further the way the sidespin curved the flight.

### Hops and the bounce-to-roll transition

As in Penner, the impact model applies to every impact, not just the
first. After the landing the ball makes a few progressively smaller hops,
then slides and rolls.

**Low hops die out.** The low hops after a landing meet nearly flat ground:
a hop arriving at 0.3–1.5 m/s of normal speed while the ball moves along
the ground at 2–4 m/s lands at 4–37°, where the factor in θ of the contact
angle is small, and for the lowest hops vanishes. On level ground each
rebound then leaves about half as fast along the normal as it arrived, so
the hops die out geometrically (see the slope analysis below). A full shot
from the bag makes 4–7 impacts (a drive 7, an iron 5, a wedge 4–5), then
slides and rolls, from about 4.5 m/s after a drive, 2 m/s after an iron
and under 1 m/s after a tour wedge.

**Shallow, fast landings skim.** At a landing of 10° the factor in θ is
0.12: a thin shot landing at 35 m/s and 10° with 2000 rpm of backspin gets
4.9° of tilt on `FAIRWAY`, leaves at 31.9 m/s and rises 1.1 m. Over landings
at 2–20° and 20–46 m/s with 0–4000 rpm of backspin on `FAIRWAY` the highest
first hop is 8.3 m (3.2 m for landings at 2–10°) and no later hop rises
above the first. A grazing landing with backspin hops higher than one
without, through the lift of the spin it leaves with, and how high is set
by `mu_slide`: the fairways' measured 0.12–0.33 applied to landings nobody
has measured.

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
g) at 8 m/s. Small as these are, they add up over a drive's hops:
through sea-level air a `FAIRWAY` drive makes 7 impacts and runs 21.1 yd
past the landing point, against 7 impacts and 23.1 yd with ballistic hops
(`GREEN`: 24.8 yd; ballistic 26.9 yd). Irons and wedges, which hop little,
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
the 5 mm height alone need not end them. Take a ball moving along the
ground at speed v and hopping with normal speed w, so that it lands at
θ ≈ w/v. For such a flat landing the contact angle is
`β ≈ compliance · (v / 50 m/s) · (θ / 20°)³`, and against the tilted plane
the impact rebounds at `w' = e·w + (1 + e)·vβ` (normal restitution; the
grip friction loss is of higher order), so each impact scales a low hop's
normal speed by

```
q = w' / w ≈ e + (1 + e) · C · w² / v,   C = compliance / (50 m/s · (20°)³)
```

with the angle in radians (C ≈ 0.49 s/m for `FAIRWAY`). With e ≈ 0.5 at
these slow impacts the tilt term is small for every low hop (0.03 for a hop
at 0.3 m/s on a ball moving at 2 m/s), so q stays near e. On level ground,
where every impact also slows the ball, the hops therefore shrink through
the 5 mm transition within a few impacts. Down a slope θ, each hop of time
`T = 2w/(g cos θ)` speeds the ball up by about `2w·tan θ`, which the next
impact partly returns as bounce, and on a steep enough slope the hops need
not shrink. Drag at these speeds is too weak to prevent it.

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
the wind refills as fast as the impacts drain them (the tilted
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
below the 16° the contact angle is measured to. Once the impacts no
longer clearly outrun the slope the ball rolls, and the rolling
resistance stands for the turf's dissipation at that speed.

### Slide and roll

On the ground the state is position, velocity and the spin vector,
integrated with RK4 at 1 ms within each smooth phase, with the terrain
sampled under the ball every step. With `n` the normal and g the model's
`gravity()`, the normal load is `m·g·n_y` and the slope pulls with the
tangential gravity `g_t = g − (g·n)n` (zero on level ground):

- **Slide** (the contact point slips): friction `−μ·m·g·n_y·û`, its torque
  `(−R n) × F`, `g_t`, and the plowing force `−plow·m·v` (below), which acts
  through the center and has no torque. On level ground without plowing the
  slip direction is fixed and the slip speed falls linearly at
  `((1+k)/k)·μ·g` (3.5·μ·g for k = 2/5), so the slide is stepped exactly to
  the moment the slip reaches zero. On a slope `g_t` turns the slip
  direction (`du/dt = −((1+k)/k)·μ·g·n_y·û + g_t − plow·v`), and plowing,
  which follows the velocity rather than the slip, makes the slip fall at a
  rate that changes with the speed, so û is re-evaluated at every RK4 stage.
  The moment the slip falls to `roll_threshold` (0.01 m/s) is solved within
  the step (regula falsi with bisection), and from there the slide finishes
  along that slip direction to zero slip, which is exact when the slip lies
  along the fall line; with plowing the moment of zero slip is solved for
  the same way. Backspin can reverse the ball.
- **Roll**: the spin tangential to the ground is locked to the velocity, and
  the ball decelerates at `(μ_eff·g·n_y + plow·v)/(1+k)` against its
  velocity, with a speed-dependent rolling resistance and plowing, while the
  slope pulls with `g_t/(1+k)`. `mu_roll` is the resistance at zero speed and
  `Surface::roll_law` (a `RollLaw`) its speed dependence,
  `μ_eff = mu_roll · f(v)` with f(0) = 1:

  ```
  RollLaw::Constant                 f = 1
  RollLaw::Exponential { v_ref }    f = exp(v / v_ref)
  ```

  `v_ref` must be nonzero (±∞ is the constant law); a zero or NaN `v_ref`
  rolls as `Constant`, and the stimpmeter helpers reject it (see
  [Stimpmeter](#stimpmeter)).

  Every preset uses `Constant` with plowing, `Surface::plow`: while the
  ball is on the ground, sliding or rolling, the turf resists it with a
  force per unit mass of `plow` times its velocity, on top of the sliding
  friction and the rolling resistance. A sliding ball slowed by it alone
  travels `v / plow`, a rolling one `(1+k)·v / plow`. It does not act
  during an impact or a hop, and as a drag it does not scale with gravity.
  With `c = plow` a ball that starts rolling at `v₀` under the constant law
  stops after `(1+k)·[v₀/c − (μg/c²)·ln(1 + c·v₀/(μg))]`, μ = `mu_roll`.

  Why a fairway needs it: the Stimpmeter fixes the resistance at 1.8 m/s,
  and rolls that start at 2.3–6.4 m/s (R&A 2008) need about twice that.
  Without plowing a 12° landing runs 143 yd where the USGA's 2019 range
  test measured 67. A constant resistance of 0.28 with no plowing fits the
  R&A roll-outs as well, but reads 2.8 ft on the Stimpmeter and misses the
  2019 line by 12.7 yd rms (8.3 yd with plowing). Only the size of `plow`
  is fitted; its form, a resistance proportional to speed, is assumed.
  `FAIRWAY`'s 0.36 /s fits the R&A bounce-and-roll totals on five fairways
  (each site's first rebound matched; 0.2–0.6 fit about as well) and the
  2019 line at 20–36°. On greens the totals of 57 landings on one firm green
  (USGA and R&A 2007) reject a constant resistance alone and accept a mild rise;
  the greens' 0.2 /s fits them best (0.1–0.3 fit about as well).

  `Constant` is Penner's run model: his run paper uses a constant rolling
  friction, `a = −(5/7)·ρ_g·g` with ρ_g = 0.131 (eq. 9), which reads a
  stimp of about 6.0 ft (ρ = 7·v₀²/(10·g·d)).

  On slopes the published checks favor `Constant`: Brede's slope
  correction of stimpmeter readings is exact under it, and Cross's
  speed-independent break angle holds under it alone. With plowing both
  are approximate (see [Stimpmeter](#stimpmeter) and
  [Terrain and slopes](#terrain-and-slopes)).

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
  Their fits cover putts only; past 3 m/s the exponential grows quickly.
- **Spin about the normal** decays at `pivot·g·n_y / (k R²)` in both
  phases, saturating at zero.
- **Rest**: a rolling ball that is slowing and slower than one step's
  deceleration (or than `stop_speed`, default 0.01 m/s) finishes at
  constant deceleration and stops. At zero speed rolling resistance can
  hold the ball against the slope where `g·sin θ ≤ mu_roll·g·cos θ`, i.e.
  `tan θ ≤ mu_roll` (`GroundSample::holds_at_rest`, on the ground where
  the ball stops): up to 5.45° on `FAIRWAY` and 2.91° on `GREEN`. There it
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

| Slope | −3° | −2° | −1° | 0° | +1° | +2° | +3° | +4° | +5° |
|---|---|---|---|---|---|---|---|---|---|
| driver | 19.5 | 19.6 | 20.2 | 21.1 | 22.4 | 24.0 | 26.3 | 28.8 | 33.2 |
| driver, ballistic hops | 21.4 | 21.4 | 22.1 | 23.1 | 24.5 | 26.3 | 28.2 | 31.4 | 34.9 |
| 7-iron (`seven_iron`) | 4.8 | 5.2 | 5.8 | 6.6 | 7.7 | 9.1 | 10.9 | 13.5 | 17.0 |
| wedge | 2.4 | 2.9 | 3.3 | 3.6 | 4.1 | 4.7 | 5.6 | 7.1 | 10.0 |

`GREEN`:

| Slope | −2.5° | −2° | −1° | 0° | +1° | +2° | +2.5° | +2.8° | +2.9° |
|---|---|---|---|---|---|---|---|---|---|
| driver | 20.6 | 21.1 | 22.6 | 24.8 | 28.3 | 34.1 | 39.0 | 43.5 | 46.0 |
| 7-iron (`seven_iron`) | 3.9 | 4.2 | 5.0 | 6.3 | 8.3 | 12.1 | 15.7 | 19.6 | 21.8 |
| wedge | 0.7 | 1.1 | 1.4 | 1.5 | 2.0 | 3.0 | 4.4 | 6.4 | 7.8 |

`simulate_shot` with `Tour::SEA_LEVEL` (hops through air); the "ballistic
hops" row is the same driver with `settle` given
`Tour::new(0.0, Wind::None)`. Downhill, the hops end by the slope rule
above and the ball rolls the rest of the way; the rule is evaluated on the
hop actually flown, so drag's loss over a hop
counts against the slope's return. On `FAIRWAY` the driver makes 6 or 7
impacts and the 7-iron 5 or 6 on every slope here; the runs grow with the
downslope mostly through the roll. Steeper than the hold angle (5.45° on
`FAIRWAY`, 2.91° on `GREEN`, uphill or down) a ball does not come to rest.

Lateral run (yards, + = right) on a plane falling to the right across the
line of play, all `Rest`; the same plane falling left gives exactly the
mirror image:

| Plane | driver | 7-iron (`seven_iron`) | wedge |
|---|---|---|---|
| `FAIRWAY` 2.5° | +4.5 | +2.2 | +1.6 |
| `FAIRWAY` 5° | +12.0 | +7.2 | +6.0 |
| `GREEN` 2.5° | +10.9 | +5.4 | +2.5 |

`GREEN` holds a ball at rest up to 2.91° (`tan θ ≤ 0.0509`). Close to
that the net deceleration of a slow roll, `(5/7)·g·cos θ·(μ_eff − tan θ)`,
is tiny and the last meters take a long time: on a 2.8° downhill green a
drive stops 29 s after landing and a 7-iron 24 s; at 2.9° they take 43 s
and 38 s, and the wedge 30 s. The default 120 s `max_time` covers these.
Beyond the hold angle a ball on `GREEN` never stops, and on a long enough
slope it runs until `max_time`.

**Landing shots on a 1.5% green.** Whitlark & Pringle (2012) fired 857
shots from a pitching machine onto a bermudagrass green sloping 1.5%,
mowed to about 12 ft on the stimpmeter, uphill and downhill, at four
launcher settings chosen to bracket the fairway and rough landings of
their Table 1 (its angles read here as descent angles from horizontal).
Mean bounce and roll was 32.5 ft downhill and 13.5 ft uphill, a 19 ft
difference. dimple on `GREEN.with_stimp(12.0)` through sea-level air runs
their fairway 7-iron landing (55 mph, 46°, 5100 rpm) 23.7 ft uphill,
31.0 ft on level ground and 41.1 ft downhill: a 17.4 ft difference, close to
the measured one, with uphill and downhill totals longer than the measured
means. Their four settings are not published, so the comparison cannot
separate the model from the shot mix. Turf firmer by 0.08 in (2 mm) on the
USGA firmness meter ran 9–12 ft further at their uphill sites; by the green
mapping of `compliance` (see [Other surfaces](#other-surfaces)) 2 mm is
0.18 of compliance, which adds 17.2 ft to dimple's uphill 7-iron, the same
order. `FIRM_GREEN`'s `mu_slide` 0.17 instead of 0.22 moves that landing
2.3 ft and matters more for high-spin landings, where the friction limit
binds. Green speeds of 11 ft 3 in and 12 ft 3 in made no difference to
their bounce and roll; in dimple the uphill 7-iron moves 0.8 ft. Their
sand wedge landing (49 mph, 51°, 9100 rpm) finished near its landing point
("in some cases this distance was negative"); dimple's stops 1.3 ft behind
it uphill, 3.3 ft past it on level ground and 6.4 ft past it downhill.
Source: Whitlark, B. and Pringle, M. (2012), "Can't Stop Your Ball on the
Green?", *USGA Green Section Record* 50(25); summary values only.

**Break across a slope.** Under a constant rolling resistance a ball
rolled across an incline stops at a displacement angle from its launch line
that does not depend on its launch speed, at a distance proportional to v₀²
(Cross 2024): the motion scales exactly with v₀. dimple reproduces this
under `RollLaw::Constant` without plowing (`mu_roll` 0.055): launched 10°
off straight uphill on a 1° or 2° plane, the ball stops 1.368° or 2.428°
off its launch line at every speed from 0.5 to 3 m/s, at 0.990 or 0.801 m
per (m/s)². `GREEN`'s plowing, which resists a faster ball more, makes the
angle drift with speed: from 1.37° to 1.06° on the 1° plane and from 2.43°
to 1.93° on the 2° plane over 0.5–3 m/s. Source: Cross, R. (2024), "Trajectory
of a golf ball on a sloping green", *Eur. J. Phys.* 45, 035002,
[doi:10.1088/1361-6404/ad2c2d](https://doi.org/10.1088/1361-6404/ad2c2d).

### Surfaces

| Surface | Restitution | `mu_slide` | `compliance` | `mu_roll` | `roll_law` | `plow` | `pivot` | `absorption` | Stimp |
|---|---|---|---|---|---|---|---|---|---|
| `Surface::FAIRWAY` | `Penner2002` | 0.18 | 1.05 | 0.0955 | `Constant` | 0.36 /s | 1 mm | 0 | 5 ft 8 in |
| `Surface::GREEN` | `Penner2002` | 0.22 | 1.10 | 0.0509 | `Constant` | 0.2 /s | 1 mm | 0 | 10.5 ft |
| `Surface::FIRM_GREEN` | `Penner2002` | 0.17 | 1.0 | 0.0509 | `Constant` | 0.2 /s | 1 mm | 0 | 10.5 ft |

| Parameter | Value | Source and measured range |
|---|---|---|
| Contact angle | V0 50 m/s, θ0 20°, power 3 | Fitted to 693 bounces on natural tee turf at 2–39 m/s and 16–90° (Biber et al. 2023); θ0 and the power checked on the USGA's 2019 roll-out line at 20–36° (USGA/R&A 2021). Ranges: V0 49.2–50.9 m/s, θ0 19.8–21.8°, power about 2–3.5. Confirmed on greens by Cross (2025), eq. 16 |
| eq. 5 restitution | | Penner (2002), from Penner (2001), *Am. J. Phys.* 69, 563, measured on relatively soft turf. Confirmed by Biber et al. (2023), Cross (2025, eq. 14) and Cross's (2021) vertical drops; the scale each measured set prefers is 0.73–1.23, mean 0.99 |
| `FAIRWAY` `compliance`, `mu_slide` | 1.05, 0.18 | Medians of eight fairways (R&A 2008: a 64 mph, 36° drive landing on five; USGA 2005: 252 landings at 23–48° on three at Pinehurst No. 2); single fairways range 0.86–1.55 and 0.12–0.33. Cross's (2021) fairway-like grass, a ninth site, fits 1.25 and 0.22 |
| `FAIRWAY` `mu_roll` | stimp 5 ft 8 in | A Stimpmeter reading on a creeping bentgrass fairway mown at 0.400 in (USGA/R&A 2018 tests; a dense tee mown at 0.300 in read 3 ft 4 in) |
| `FAIRWAY` `plow` | 0.36 /s | Fitted to the R&A (2008) bounce-and-roll totals on five fairways and the USGA's 2019 roll-out line; 0.2–0.6 fit about as well |
| `GREEN` `compliance`, `mu_slide` | 1.10, 0.22 | Medians of nine greens (R&A 2008: 5-iron and sand wedge landings on six; USGA 2005: 4-iron to pitching wedge on three); range 1.05–1.52 and 0.17–0.26. Fitted to Cross's (2025) printed fits to 184 bounces on four greens alone: 1.12 and 0.23 |
| `FIRM_GREEN` `compliance`, `mu_slide` | 1.0, 0.17 | 57 impacts at 20–26 m/s and 36–59° on one green dried to championship firmness, Stimpmeter 10–11 ft (USGA and R&A 2007). Cross (2025) gives a friction of 0.6 on a firm green; `FIRM_GREEN` reads 0.56 there |
| Greens `mu_roll` | stimp 10.5 ft | The reading of that green (10–11 ft): the one green with both a Stimpmeter reading and measured roll-outs. A setting, independent of the impact constants (`with_stimp`). Cross (2025) quotes rolling friction of 0.04–0.09 for fast to slow greens |
| Greens `plow` | 0.2 /s | Fitted to that green's bounce-and-roll totals; 0.1–0.3 fit about as well |
| Stimpmeter exit speed | 1.83 m/s | Penner (2002), "The physics of putting", *Can. J. Phys.* 80, 83, from Holmes (1986), *The Physics Teacher* 24, 401 |
| 5 mm transition height | | Penner (2002), set above the ~2 mm equilibrium sinkage of a ball on a green (Hubbard & Alaways) |
| `RollLaw::Constant`, `RollLaw::Exponential` | | `Constant`: Penner (2002), eq. 9. `Exponential`: the form of Burritt, A. and Henrikson, E. (2024), "Modeling golf ball roll with a variable coefficient of friction", *Sports Engineering* 27:39, [doi:10.1007/s12283-024-00481-5](https://doi.org/10.1007/s12283-024-00481-5) (abstract only consulted). Open conference version: *Engineering of Sport 15* (ISEA 2024), Loughborough University repository, [doi:10.17028/rd.lboro.27051061.v1](https://doi.org/10.17028/rd.lboro.27051061.v1), CC BY 4.0, from which the surface values quoted here are read off Figure 1 |
| `pivot` | 1 mm | Hertzian full-slip torsional couple (Johnson 1985) with an assumed 4 mm contact radius; an estimate, uncalibrated |
| `absorption` | 0 | Turf returns the whole rebound |

**`GREEN`'s friction range.** `mu_slide` is fitted to landings at
17–27 m/s. One green measured at 9–10 m/s (Cross 2025, Table 1) slides at
0.29; a chip with 3,300 rpm of backspin leaves that green 1.3 m/s slower
than this preset gives. The same law fits those chips with that green's own
constants, so the gap is friction that falls with speed, which dimple does
not model.

The run is sensitive to `compliance`. Fairway run (yards past carry, sea
level) by `compliance`, other fields as `Surface::FAIRWAY`:

| `compliance` | driver | 7-iron (`seven_iron`) | wedge |
|---|---|---|---|
| 1.25 | 12.5 | 2.6 | 1.0 |
| 1.15 | 16.8 | 4.4 | 2.3 |
| 1.05 (`FAIRWAY`) | 21.1 | 6.6 | 3.6 |
| 0.95 | 25.5 | 9.1 | 5.2 |
| 0.85 | 30.1 | 11.8 | 7.0 |
| 0.75 | 34.7 | 14.7 | 9.1 |
| 0.50 | 44.0 | 22.2 | 14.9 |

Sources: Biber et al. (2023), Cross (2025) and Cross (2021), cited under
[Validation against published bounce data](#validation-against-published-bounce-data);
Johnson, Otto, Richardson and Isaac (2008), *An Investigation into
Ball-Turf Impact Characteristics*, The R&A
([PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/2887f275-0659-4812-a1be-30ab00a8472c/R04%20-%20An%20Investigation%20into%20Ball-Turf%20Impact%20Characteristics%20_Master%20Copy_%20_5.pdf));
USGA Research and Test Center (2005), *Preliminary Results of Turf Testing
at Pinehurst #2, Spring 2005*
([PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/df3cf55d-4c47-4dbd-8600-3c10d121d248/R45%20-%20Preliminary%20Results%20of%20Turf%20Testing%20at%20Pinehurst.pdf));
USGA and The R&A (2007), *Second Report on the Study of Spin Generation*,
appendix on golf ball impact with a green surface (test of 30 October 2006)
([PDF](https://www.jga.or.jp/jga/html/rules/image/info_groove_second_report_spin_generation_en.pdf));
R&A and USGA (2020), *R03: Agronomic Impacts on Bounce and Roll Distance*
(tests of 2018), and USGA and The R&A (2021), *Proposed Bounce Model for Use
in Evaluating Optimum Overall Distance*, cited under
[Run-out on a fairway](#run-out-on-a-fairway).

### Run-out on a fairway

Run (yards past carry) for launches across the bag in Trackman's 25 °C
sea-level air, with the descent angle at landing:

| Launch (mph / launch° / rpm) | Lands at | `FAIRWAY` | `GREEN` |
|---|---|---|---|
| Lob wedge, tour (72 / 32 / 9000) | 46.6° | 6.7 | 4.2 |
| Lob wedge, amateur (65 / 32 / 6500) | 43.9° | 9.3 | 9.2 |
| Sand wedge, tour (80 / 30 / 9500) | 47.3° | 6.2 | 3.6 |
| Sand wedge, amateur (72 / 30 / 7500) | 44.3° | 8.6 | 7.3 |
| Pitching wedge, tour (102 / 24.2 / 9300) | 48.7° | 5.2 | 3.0 |
| Pitching wedge, amateur (85 / 26 / 7500) | 44.4° | 8.7 | 7.8 |
| 50-yard pitch (52 / 30 / 6000) | 37.4° | 12.5 | 12.7 |
| 50-yard pitch, less spin (52 / 32 / 4000) | 39.1° | 14.8 | 17.1 |
| Lob-wedge chip (28 / 36 / 3000) | 38.2° | 9.5 | 11.5 |
| Low chip (30 / 18 / 2500) | 19.0° | 24.7 | 30.7 |
| 7-iron, tour (120 / 16.3 / 7100) | 45.5° | 8.9 | 8.7 |
| 7-iron, amateur (105 / 17 / 5900) | 40.0° | 14.6 | 15.8 |
| Driver, tour (171 / 10.4 / 2545) | 38.8° | 23.0 | 26.8 |

Every launch in the table makes 5–7 impacts before it slides and rolls.
Published measurements of run on a fairway:

- **The USGA's 2019 range test** fired balls at 75–90 ft/s and 1,800 rpm
  onto a fairway at landing angles of 10–36°; their roll-out follows
  `88.3 − 1.79·α` yd. Averaging landings at 75 and 90 ft/s, `FAIRWAY` is
  within 3 yd of that line at 20–36° (rms 2 yd) and about 11 yd long at
  12–18° (13.2 and 11.0 yd at 12° and 16°), below the 16° to which the
  contact angle is measured; at 75 ft/s alone it is about 5 yd long there.
- **The R&A's driver landing** (2008; 64 mph, 36°, 2,250 rpm) bounced and
  rolled 4–31 yd on five fairways by firmness (median 24 yd); `FAIRWAY`
  gives 24.8 yd.
- **The PGA Tour average drive** (Trackman's 171 mph, 10.4°, 2545 rpm,
  landing at 39°) lands at 38.8° in dimple and runs 23.0 yd. PGA TOUR
  drives of 2016–18 roll out 21.7 yd at 39° on the mean line over about 80
  events (10–60 yd from event to event), and the R&A/USGA agronomy study
  measured 17.8 ± 13.3 yd over 624 fairway shots landing at 32–50° with
  about 1800 rpm.
- Wedge, pitch and chip landings (39–49° here) lie inside the measured
  impacts, but no published roll-out covers them on a fairway.

High-spin wedges run a little farther on `FAIRWAY` than on `GREEN`: the
greens grip harder and press deeper, so those landings leave a green
slower, and the greens' lower rolling resistance makes up the difference
only for landings that roll farther. `mu_roll` (or `with_stimp`) and
`plow` are the knobs for a faster or slower fairway: the tour drive,
amateur lob wedge and less-spun pitch run 23.7, 9.6 and 15.3 yd at
`mu_roll` 0.08 and 19.3, 8.3 and 12.4 yd at 0.25.

Sources: USGA and The R&A (2021), *Proposed Bounce Model for Use in
Evaluating Optimum Overall Distance*
([PDF](https://assets-us-01.kc-usercontent.com/c42c7bf4-dca7-00ea-4f2e-373223f80f76/250b3e34-9172-4ad1-991f-9d5560cbf832/Proposed%20Bounce%20Model%20for%20Use%20in%20Evaluating%20Optimum%20Overall%20Distance.pdf)),
which also reports the ShotLink runs; R&A and USGA (2020), Distance
Insights Project, *R03: Agronomic Impacts on Bounce and Roll Distance*
([PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/cdd48ad2-2ed1-4e67-ad53-1a6522110779/R03%20-%20Agronomic%20Impacts%20on%20Bounce%20and%20Roll%20Distance.pdf));
Johnson, Otto, Richardson and Isaac (2008), *An Investigation into
Ball-Turf Impact Characteristics*, The R&A
([PDF](https://assets.randa.org/c42c7bf4-dca7-00ea-4f2e-373223f80f76/2887f275-0659-4812-a1be-30ab00a8472c/R04%20-%20An%20Investigation%20into%20Ball-Turf%20Impact%20Characteristics%20_Master%20Copy_%20_5.pdf));
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
component. dimple's contact angle is fitted to Campaign B's 693 bounces
themselves (with compliance 1.0 and `mu_slide` 0.35, the tee constants);
the scores below measure the model against the paper's own summary of the
same bounces.

Natural tee (Campaign B):

| Against the tee fit | δ | Δẋ | Δẏ | Δω |
|---|---|---|---|---|
| tee constants (`compliance` 1.0, `mu_slide` 0.35) | 10.5% | 5.5% | 2.7% | 7.2% |
| `GREEN` | 16.8% | 6.2% | 2.7% | 14.2% |
| `FAIRWAY` | 20.6% | 7.1% | 2.9% | 17.8% |
| tee constants on a rigid plane (`compliance` 0) | 21.9% | 15.5% | 6.7% | 11.4% |
| tee constants with `mu_slide` 1.0 | 22.0% | 9.9% | 3.1% | 18.2% |

For comparison, against the data the paper's rigid-bounce fit scores
35.0%, its Penner-model fit 21.3% and its fixed-angle fit 19.2%. The
presets score worse than the tee constants on this turf through their
lower friction, which is measured on greens and fairways. With the tee
constants the mean normal rebound ratio |ẏ_F / ẏ₀| is 0.371 against the
fit's 0.408, and of the 80 backspin landings 19 leave with topspin against
the fit's 23. Most of the remaining error is tangential: dimple keeps about
60% of the forward speed the fit keeps (mean ẋ_F / ẋ₀ 0.191 against
0.313). Balls that enter rolling leave the ground slipping with excess
topspin, contact-point velocity ẋ + ω < 0 (the paper's Fig. 6); dimple's
tilted contact plane gives that sign in the ground frame for 55 of the 57
landings on the fit's rolling branch, with half the fitted mean (−29
against −58 radii/s), but not the spread of the measurements, which the
authors put down to the natural variability of turf; they suggest that the
tangential and spin error left by rigid-contact models needs a tangential
compliance of the ground, which dimple does not model. The fits' μ ≈ 1 is
not adopted: it grips every landing, while the slipping bounces leave with
their backspin, and in dimple it doubles the error. Their typical driver
landing (28.5 m/s, 37.3°, 34.9 rev/s) rebounds from `GREEN` with a normal
ratio of 0.356 against the tee fit's 0.342 and a forward ratio of 0.236
against 0.277, leaving with topspin as in the fit.

Artificial mat (Campaign A): `Penner2002` restitution is far too low for a
mat on a rigid base. `GREEN`'s mean normal rebound ratio is 0.34 against
the fit's 0.58 (the paper fits r = 0.42–0.54 for the mat, Table 3), and
`GREEN` scores δ 29.4% with Δẏ 17.1%; no change of `compliance` or
`mu_slide` alone brings Δẏ down. A constant restitution does: the mat
example under [Other surfaces](#other-surfaces) scores δ 16.8% with Δẏ
3.7%. Use `Restitution::Constant` for mats.

Source: Biber, S. W., Jones, K. M., Champneys, A. R., Green, R. and
Szalai, R. (2023), "Measurements and linearized models for golf ball bounce
on a green", *Sports Engineering* 26:50,
[doi:10.1007/s12283-023-00442-4](https://doi.org/10.1007/s12283-023-00442-4),
CC BY 4.0. The fitted coefficients quoted here are from its Table 3 and
eqs (8)–(9); dimple's figures are its own runs against those fits.
Cross, R. (2025), "Bounce of a golf ball on the green", *Proc. IMechE
Part P* 239(3), 479–492,
[doi:10.1177/17543371221150289](https://doi.org/10.1177/17543371221150289);
Cross, R. (2021), "Oblique impact of a hard ball on a soft surface",
*Eur. J. Phys.* 42, 065005,
[doi:10.1088/1361-6404/ac221e](https://doi.org/10.1088/1361-6404/ac221e).

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
Exponential:  I = 2·v_ref² · (1 − e^(−s₀)·(1 + s₀))
```

The exponential form tends to `v₀²` as `v_ref` → ±∞ and works for a
negative `v_ref` unchanged.

With plowing, `c = plow` > 0, the deceleration is
`(mu_roll·f(v)·g + c·v)/(1+k)` and
`d = (1+k)·∫₀^v₀ v/(mu_roll·f(v)·g + c·v) dv`: for `Constant`,
`d = (1+k)·[v₀/c − (μg/c²)·ln(1 + c·v₀/(μg))]` with μ = `mu_roll`; for the
other laws the integral is taken numerically. `d` falls as `mu_roll`
grows, and `mu_roll` is found by bisection. Plowing alone stops the ball
within `(1+k)·v₀/c`, so a reading at or beyond that (23.3 ft for
`FAIRWAY`'s 0.36 /s, 42.0 ft for the greens' 0.2 /s) has no `mu_roll`.

`Surface::with_stimp(feet)` sets `mu_roll` by this closed form (by
bisection with plowing), keeping every other field (including `roll_law`
and `plow`), and `Surface::stimp_feet()` is the inverse. Both use k = 2/5,
the uniform solid sphere of `Tour`. `with_stimp` returns the surface
unchanged for a reading that is not finite and positive, a law whose
`v_ref` is zero or NaN, for which `stimp_feet` is NaN, or a reading no
`mu_roll` ≥ 0 gives.

| Stimp (ft) | 6 | 8 | 10 | 12 | 14 |
|---|---|---|---|---|---|
| `GREEN.with_stimp(…).mu_roll` | 0.1064 | 0.0740 | 0.0546 | 0.0418 | 0.0327 |

A simulated stimpmeter roll (rolling at 1.83 m/s on level ground) on each
of these surfaces stops within 1e-7 (relative) of the requested distance,
and within 1e-9 without plowing.

**On a slope.** Brede (1991) showed that the readings a up and b down a
uniform slope combine as 2ab/(a+b), which reads the level value on slopes
from 0 to 5.6% (a laboratory carpet runway) and on field greens up to 6%;
the plain average of the two does not. Constant rolling resistance
decelerates the ball at g(μ cos θ ± sin θ)/(1+k), so under
`RollLaw::Constant` the correction is exact apart from the 1/cos θ of
measuring along the surface (0.005% at 1%, 0.18% at 6%). With the presets'
plowing the up- and downhill rolls cover different speeds and it is
approximate:

| Grade | 1% | 2% | 3% | 4% | 5% | 6% |
|---|---|---|---|---|---|---|
| `GREEN` up / down (ft) | 9.25 / 12.16 | 8.26 / 14.47 | 7.47 / 17.92 | 6.82 / 23.78 | 6.27 / 38.23 | down does not rest |
| `GREEN` 2ab/(a+b) vs level 10.50 ft | +0.04% | +0.16% | +0.40% | +0.92% | +2.6% | – |
| `FAIRWAY` 2ab/(a+b) vs level 5.67 ft | +0.01% | +0.05% | +0.12% | +0.22% | +0.36% | +0.55% |

On `GREEN` the plain average is already 8% long at 2%. `GREEN` holds a
ball at rest only to a 5.09% grade (`tan θ ≤ mu_roll`), so a downhill
stimpmeter roll on a steeper `GREEN` never stops, and an uphill one stops
once and rolls back.

Source: Brede, A. D. (1991), *Agron. J.* 83:425–426,
[doi:10.2134/agronj1991.00021962008300020032x](https://doi.org/10.2134/agronj1991.00021962008300020032x).

### Other surfaces

Beyond the presets above, dimple adds a surface preset only when its
parameters come from a published source. Build other surfaces from a
preset with struct update syntax and, for the roll, a stimp reading.

**Firmness is `compliance`.** On the R&A's six greens (2008) the fitted
`compliance` follows the penetration depth of the USGA firmness meter:
`compliance ≈ 0.090 × depth (mm) + 0.25` over 8.8–14.3 mm (r = 0.96;
residual 0.04), so each millimeter firmer is 0.09 less compliance and
about 2° less tilt on a 5-iron landing. On their five fairways the
relation is indicative only (`0.074 × depth + 0.23` over 8.9–15.4 mm,
r = 0.70; residual 0.17). Softer greens also rebound less, which is not
modeled.

```rust
// A soft, slow fairway: rolling like stimp 4 ft, penetration about 15 mm.
let soft_fairway = Surface { compliance: 1.35, ..Surface::FAIRWAY }.with_stimp(4.0);
// A firm, fast fairway: penetration about 9 mm.
let firm_fairway = Surface { compliance: 0.9, ..Surface::FAIRWAY }.with_stimp(7.0);
// A firm green: penetration about 8 mm (championship firmness).
let firm_green = Surface { compliance: 0.97, ..Surface::GREEN }.with_stimp(12.0);
```

**Rough and sand** have no published measurements this model's
constants could be fitted to: no published test fires balls into rough,
and none measures a ball slowing through grass. The one published figure
is that PGA TOUR drives finishing in rough lose 8–9 yd of bounce and roll,
with a weaker dependence on landing angle (−0.45 yd per degree against
−1.24; Quintavalla & Rogers 2008, in USGA/R&A 2021). The README's "Rough
and sand" section gives starting values to tune by feel; the sand example
below shows what `absorption` does.

An artificial mat on a rigid base rebounds much more than turf, which the
`Penner2002` restitution does not give (see
[Validation against published bounce data](#validation-against-published-bounce-data)).
This example scores δ 16.8% against Biber et al.'s mat fit (`GREEN`: 29.4%).
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
sets the roll at stimpmeter speed and `plow` how much more a fast ball is
resisted.

**Sand.** A ball does not bounce on sand; it skips. Dropped on loose
sand it comes to rest part-buried whatever its own restitution (Uehara et
al. 2003), and bunker sand is chosen so that a ball stays where it lands
"without bouncing" (The R&A, *Golf Course Condition and Playability*). A
ball that lands at an angle keeps part of its motion along the surface and
travels on, with most of its speed and spin gone, and the steeper it lands
the less it keeps: 0.80 of its speed at 10° from the surface, 0.30 at 30°
and 0.008 at 50° for a small projectile at 104 m/s (Suo et al. 2023); about
0.14 for a marble at 2–5 m/s, which stops in its crater above about 50°
(Wright et al. 2020). A golf ball fired into bunker sands at 71° stays in
its mark, half to fully buried (Baker 2001).

```rust
// Bunker sand: no bounce. The ball skips on along the sand, and the steeper it lands the less it keeps.
let sand = Surface {
    restitution: Restitution::Constant(0.0),
    compliance: 0.0,
    absorption: 0.9,
    mu_slide: 0.40,
    mu_roll: 0.3,
    ..Surface::FAIRWAY
};
```

- `restitution` 0 returns no normal speed and `compliance` 0 keeps the
  contact plane flat, so the ball does not hop: it leaves along the surface.
- `absorption` is the share of a rebound the sand takes from a ball landing
  straight down. The sand resists the ball's motion through it with a drag
  proportional to v² (Katsuragi & Durian 2007, 2013), which leaves
  `e^(−s/d1)` of the velocity after a path `s`, and the ball drives in as
  deep as the normal part of its velocity takes it (Wright et al. 2020). A
  landing at angle θ from the surface therefore keeps
  `(1 − absorption)^(sin θ)` of its velocity and spin. With 0.9, a ball with
  no spin leaves with 0.61, 0.37, 0.21, 0.13 and 0.08 of its speed at 10°,
  20°, 30°, 40° and 50°, between the two experiments from 20° to 40°.
- `mu_roll` 0.3 is the resistance of a sphere of the ball's density rolling
  on a granular bed (Fukumoto et al. 2026: about 0.25; De Blasio & Saeter
  2009: 0.45–0.5 for glass). `mu_slide` 0.40 is chosen; `roll_law`
  (constant), `plow` (0.36 /s) and `pivot` are `FAIRWAY`'s. No source gives
  values for a ball on sand.

On this surface no landing hops:

| Landing | Rest past the mark | Speed out | Spin in → out |
|---|---|---|---|
| Sand wedge, 49 mph, 51°, 9100 rpm | 0.20 m | 1.17 m/s | 953 → 26 rad/s |
| 7-iron, 55 mph, 46°, 5100 rpm | 0.68 m | 1.91 m/s | 534 → 56 rad/s |
| Preset driver, 28.5 m/s, 40° | 2.03 m | 3.30 m/s | 211 → 147 rad/s |
| Chip, 8 m/s, 30°, 3000 rpm | 0.41 m | 1.68 m/s | 314 → 40 rad/s |
| Lob, 15 m/s, 65°, 6000 rpm | 0.00 m | 0.11 m/s | 628 → 1 rad/s |
| Thin, 35 m/s, 10°, 2000 rpm | 33 m | 21.5 m/s | 209 → 51 rad/s |

Full shots leave a short track, a little longer the faster and shallower
they land. A steep landing trickles a short way (the wedge) or stops on its
mark (the lob); it does not plug, and dimple reports where the ball stops,
not how deep it sits. A thin shot keeps most of its speed and skips out of a bunker: the
33 m above is on endless sand. With `absorption` 0.95 the sand is softer
(wedge 0.07 m, 7-iron 0.27 m, driver 0.90 m); with 0.8 it is firmer
(0.57, 1.71 and 4.45 m). It is an example, not a preset: no published
measurement gives a golf ball's rebound, spin or run on bunker sand, the
two oblique experiments used other projectiles, and 0.9 is a choice
between them.

Sources: Uehara, Ambroso, Ojha and Durian
(2003), "Low-speed impact craters in loose granular media", *Phys. Rev.
Lett.* 90, 194301,
[doi:10.1103/PhysRevLett.90.194301](https://doi.org/10.1103/PhysRevLett.90.194301);
The R&A, *Golf Course Condition and Playability: Sustainability Guidelines
and Support for Clubs*, §7 Bunkers; Suo et al. (2023), "Subsurface
pulse, crater and ejecta asymmetry from oblique impacts into granular
media", *Icarus*,
[doi:10.1016/j.icarus.2023.115816](https://doi.org/10.1016/j.icarus.2023.115816);
Wright, Quillen et al. (2020), "Ricochets on asteroids:
Experimental study of low velocity grazing impacts into granular media",
*Icarus* 351, 113963,
[doi:10.1016/j.icarus.2020.113963](https://doi.org/10.1016/j.icarus.2020.113963);
Baker, S. W. (2001), "True grit", *Greenkeeper International*, May 2001,
34–36; Katsuragi and Durian (2007), "Unified force law for
granular impact cratering", *Nature Physics* 3, 420,
[doi:10.1038/nphys583](https://doi.org/10.1038/nphys583), and (2013), "Drag
force scaling for penetration into granular media", *Phys. Rev. E* 87,
052208, [doi:10.1103/PhysRevE.87.052208](https://doi.org/10.1103/PhysRevE.87.052208);
Fukumoto, Ebata, Sharma and Katsuragi (2026), "Depth and
slip ratio dependencies of friction for a sphere rolling on a granular
slope", [arXiv:2602.01652](https://arxiv.org/abs/2602.01652); De Blasio
and Saeter (2009), "Rolling friction on a granular medium",
*Phys. Rev. E* 79, 022301,
[doi:10.1103/PhysRevE.79.022301](https://doi.org/10.1103/PhysRevE.79.022301).

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

- **Natural-turf data only.** The contact angle is fitted on natural tee
  turf and the presets' constants on natural greens and fairways. An
  artificial mat on a rigid base needs `Restitution::Constant` (see
  [Other surfaces](#other-surfaces)).
- **No impact is measured below a 16° landing.** The contact angle there is
  the fitted law carried below its data. Roll-out after landings at
  12–18° runs about 11 yd longer than the USGA's 2019 range test (about
  5 yd if those landings were at its slowest speed, 75 ft/s), and the power of
  the landing-angle factor is bounded (about 2–3.5) rather than pinned.
- **The rebound is too lively when shallow.** It is 0.07–0.15 livelier than
  measured for landings flatter than about 35° (Biber et al. 2023; Cross
  2021). Applying the restitution to the true normal speed was tested and
  fits worse.
- **Friction does not fall with speed.** One green measured at 9–10 m/s
  slides at 0.29 (Cross 2025, Table 1), and Cross's eq. 15 rises 0.0095 per
  m/s where the model rises 0.0150; `GREEN`'s 0.22 is fitted at 17–27 m/s,
  so chips with heavy backspin leave the green faster than measured
  (1.3 m/s at 10 m/s and 3,300 rpm). On the USGA's
  firm green high-spin landings that stop by their mark roll on about
  5 ft in dimple, and rolls over 50 ft come up about 9 ft short.
- **One site differs from the next more than the model errs.** Single
  fairways range 0.86–1.55 in `compliance` and 0.12–0.33 in `mu_slide`,
  greens 1.05–1.52 and 0.17–0.26; the presets are medians, and the firmness
  mapping of `compliance` is measured on greens and only indicative on
  fairways.
- **Plowing's form is assumed.** The resistance proportional to speed is a
  form; only its size is fitted, to roll-outs that start at 2–6 m/s. The
  fairway's Stimpmeter reading comes from a different fairway than the
  roll-outs. The published speed-dependent fits of the rolling resistance
  itself (Burritt & Henrikson 2024) cover putts, 0–3 m/s;
  `RollLaw::Exponential` takes such a per-surface fit, and extrapolated
  past 3 m/s it grows quickly (`v_ref` = 1.6 m/s gives about 150 ×
  `mu_roll` at 8 m/s, so a fast roll stops almost at once), and with a
  negative `v_ref` a ball has no terminal speed on a slope steeper than its
  hold angle.
- **Tangential compliance is not modeled.** Measured bounces leave the ground slipping, with more topspin than rolling,
  and the authors suggest that the tangential and spin error of rigid-contact
  models needs a tangential compliance of the ground (Biber et al. 2023).
  dimple's tilted rigid contact plane gives the sign for balls that grip, with
  about half the fitted mean on natural turf, but not the spread of the
  measurements, and keeps about 60% of the forward speed the tee fit keeps.
- **Lateral ground effects are unvalidated.** The kick of a curving ball,
  the sidespin throw at impact and rifle-spin drift follow from the
  rigid-sphere contact model; no published measurement is compared against.
  The pivot lever is uncalibrated.
- **A fast plowing roll is held rolling.** Rolling against `plow` needs a
  friction of k/(1+k) · `plow` · v, which passes `mu_slide` · g above
  (1+k)/k · `mu_slide` · g / `plow`: about 17 m/s on `FAIRWAY` and 38 m/s
  on `GREEN`. Faster than that a real ball would slide with overspin;
  dimple keeps the ball rolling from its first zero slip.
- **Sand has no stop angle, no skip through the air and an assumed effect
  on spin.** Loose sand stops a projectile that lands steeper than about
  50° (Wright et al. 2020; Suo et al. 2023); dimple's absorption only
  shrinks toward a straight-down landing, so a steep ball with little
  backspin still runs a fraction of a meter. A shallow fast landing skips
  along the surface, slowed by `mu_slide` and `plow`, where a real ball may
  leave the sand for a moment. The absorption does not depend on speed,
  and it takes the same share of the spin as of the velocity, which no
  published measurement covers.
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
- **Slopes: rolling checked, landings checked on one green.** Rolling on a
  slope (normal load, downhill pull, hold where `tan θ ≤ mu_roll`) is
  checked against Brede's stimpmeter slope correction and Cross's
  break-angle invariance: exact under `RollLaw::Constant` without plowing,
  approximate with the presets' plowing. The impact model on a slope is
  checked only against the 1.5%-green means of Whitlark & Pringle (2012),
  whose shot mix is unpublished; dimple knows of no published
  bounce-and-spin measurements on steeper slopes. The impact constants are
  fitted on level turf.
- **The slope hop rule is a modeling choice.** It ends hop trains that
  the impact model would sustain down a slope (see
  [Hops and the bounce-to-roll transition](#hops-and-the-bounce-to-roll-transition)).
  Its margin, one 5 mm hop's energy per impact, carries over Penner's
  transition height; it is not measured.
- **Slopes beyond the hold angle never come to rest.** Steeper than
  `atan(mu_roll)` (2.91° on `GREEN`, 5.45° on `FAIRWAY`) the ball rolls on
  until the ground changes or `max_time` ends it with `MaxTime`. Just
  inside the hold angle the last meters of a slow roll take most of a
  minute (43 s for a drive on a 2.9° `GREEN`; see
  [Terrain and slopes](#terrain-and-slopes)); a `max_time` set below the
  120 s default can cut such rolls short.

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

- **Fit doc** — `examples/fit_doc.rs` recomputes the rows of [FIT.md](FIT.md) built from printed numbers (Cross 2025 eqs. 14–16 and Table 1, the USGA 2007 fits, the Stimpmeter readings, the 2019, Philadelphia and PGA TOUR roll-out lines) and fails when they differ from `docs/fit-summary.csv` by more than rounding
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
- **Ground: planar regression** — with ballistic hops, the 3D impact reproduces an independent planar solution of the contact-angle law impact by impact (preset driver, 7-iron and wedge on both presets) to float tolerance, and the rest point with each surface's roll law and plowing
- **Ground: physics** — Penner (2002) restitution, 5 mm bounce-to-roll transition, the contact angle (14.6°, 20.1°, 20.7° and 14.9° at the four landings Penner (2002) publishes at compliance 1.0, `atan(v_t / 50 m/s)` on steep landings whatever the normal speed, 4.63° at 35 m/s and 10°, scaled by `compliance`, zero as the normal speed vanishes), compliance tilt (2β rotation, incidence cap), the low skimming hops it shrinks, a 46 m/s landing at 15° on all three presets that never turns backward and whose first hop stays under 10 m, shallow landings at 2–10° on all three presets skipping forward with no rebound faster than its arrival and each hop no higher than the last, `FAIRWAY` bouncing lower than `GREEN` and rolling shorter from the same landings, backspin reducing forward speed up to the friction bound, a high-spin wedge spinning back on `GREEN`, grip-regime kick formula, no kick from spin about the normal, tilt-throw magnitude (grip and slip), rifle-spin drift, grazing impacts, pivot friction about the true normal and its independence from the level-ground path, general inertia factor (k = 2/3), exact slide-to-roll time, stimpmeter closed form (the constant law without plowing, and with plowing on both presets), `stimp_feet` with plowing against its closed form and the longest reading plowing allows, a fast roll and a fast slide on `FAIRWAY` that never speed up or reverse and rest within 1% of the closed form, slides with heavy plowing or sliding friction that end at the exact zero-slip time, `with_stimp`/`stimp_feet` round trip under every roll law and simulated stimpmeter rolls (6–14 ft) to 1e-9 without plowing and 1e-7 with it, preset stimp readings (5 ft 8 in and 10.5 ft)
- **Ground: roll laws** — every law is `mu_roll` at zero speed, exponential roll-outs against the closed form (positive and negative `v_ref`), Penner's ρ_g = 0.131 reads 6.0 ft under the constant law, an infinite `v_ref` is the constant law, a rising exponential stops a fast roll sooner than the constant law at the same stimp, zero or NaN `v_ref` rejected by the stimp helpers and rolled as the constant law
- **Ground: hops through air** — at zero density hops are bit-for-bit ballistic (on the parabola, `Surface` and `GroundPhysics` identical); through air a hop arrives earlier, nearer and slower with its spin unchanged; the run lost to drag shrinks with altitude; topspin pushes a hop down and backspin holds it up; sidespin curves a hop right (fade spin) or left, mirrored
- **Ground: invariants** — energy never increases at an impact or along the ground phase, with ballistic hops and through air; mirror (±sidespin, ±azimuth) and rotation about Y give mirrored and rotated results; fades rest right, draws rest left
- **Ground: terrain** — rest on slopes (hold below atan `mu_roll`, roll downhill above, both presets; one stop and no chatter rolling uphill), long downhill reports `MaxTime` at its plowing terminal speed, near-critical (2.9°) `GREEN` comes to rest within the default `max_time`, slide-to-roll on a plane against the closed form, slip turning across a slope, vertical spin kicks on slopes only, uphill landings check up and downhill landings release, elevated green (shorter carry, lands on it), tee on the terrain, landing and hop crossings on curved terrain to 1e-9 m (ballistic hops on the parabola; hops through air on the surface), first crossing of thin walls and ridges in a landing segment and a hop step, ball kept on the surface; slope hop rule, downhill `FAIRWAY` up to 5° (through air and ballistic, the driver run growing with the slope either way, shorter through air) and `GREEN` up to 2.5° end at rest, side slopes rest and mirror exactly, rest judged on the surface where the ball stops, slide-to-roll on curved ground without flicker; mirror and rotation about Y with a sloped plane through air, and energy on slopes at both densities
- **Ground: validation against published data** — Brede's (1991) slope correction exact under the constant law without plowing, within 0.4% to 3% and 0.8% to 4% with a green's plowing and within 0.4% to 6% with a fairway's, a 5.5% hold grade at `mu_roll` 0.055 and the plowing terminal speed past it; agreement of the tee constants with their natural-tee fit (δ, Δẏ, normal rebound, spin reversal, lift-off slip sign; a rigid plane and μ = 1 do worse; the presets within 25%); their mat needs a higher constant restitution; their typical driver landing; Whitlark & Pringle's (2012) 1.5% green (uphill/downhill asymmetry, the sand wedge stopping near its pitch mark and behind it uphill, green speed, `compliance` as firmness through the green penetration mapping); Cross's (2024) speed-independent break angle under the constant law; the USGA's 2019 roll-out line on `FAIRWAY` within 3 yd at 20–36° (rms 2 yd); the R&A's 2008 driver landing bouncing and rolling within the five fairways' 4–31 yd; the USGA's 2006 firm-green printed fits on `FIRM_GREEN` (backspin out within 1,250 rpm, first bounce within 7 ft)
- **Ground: absorption** — the presets absorb nothing; `impact` with `absorption` equals the turf rebound times `(1 − absorption)^(sin θ)` bit for bit; full absorption stops the ball on its mark; on the sand example no landing hops, every landing ends at rest after one contact with no energy gained, steeper landings keep less speed and stop sooner, the sand wedge keeps under 5% of its spin, full shots leave a 0.3–3 m track, and a thin shot leaves faster than 15 m/s
- **Run-out** — the hop train dies out (at most 10 impacts, each hop after the first at most 0.7× the one before, at most two below 1 cm) from wedges to drives on both presets; `FAIRWAY` roll and total for a lob wedge, a pitch and the PGA Tour average drive; the preset drive's carry, run and total; the drive runs shorter on `FAIRWAY` than on `GREEN`; a topped drive (130 mph, 2°) skims forward on `FAIRWAY` and stays below 0.5 m
- **Ground: API** — per-point surfaces (fairway-to-green roll), `MaxTime`/`MaxImpacts` termination, tunables, interpolated landing, `&dyn GroundModel`

Run: `cargo test`
