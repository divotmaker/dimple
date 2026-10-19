# Golf Ball Aerodynamics Model

Aerodynamic model behind `dimple`, for any project that needs realistic ball
flight from launch monitor data.

## Overview

The model computes drag coefficient (Cd) and lift coefficient (Cl) as functions
of Reynolds number and spin ratio, then integrates the trajectory using RK4 at
a 1ms fixed timestep. It handles altitude, temperature, and spin decay.

**Key inputs:** ball speed, launch angle, backspin RPM, sidespin RPM, altitude,
temperature.

**Key outputs:** carry distance, apex height, lateral deviation, flight time,
full trajectory points; with a ground model, every impact, the hops, slide
and roll, and the rest position (see [Ground Model](#ground-model)).

## Conventions

Internal units are SI: meters, seconds, kilograms, radians, rad/s.

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
| lateral deviation (`lateral_m`, `lateral_yards`) | ball lands right of target |
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
start line) and points toward +Z (right) for positive sidespin. The axis is
fixed in flight; only the rate decays. `FlightResult` reports the landing
spin vector and its exact projections on the launch backspin axis
(`landing_backspin_rads`) and on -Y (`landing_sidespin_rads`).

## Force Model

Three forces act on the ball at each timestep:

```
F_total = F_gravity + F_drag + F_magnus
```

- **Gravity:** constant -9.80665 m/s² in Y
- **Drag:** opposes velocity. Fd = -Cd × q × A / m × v̂
- **Magnus lift:** perpendicular to velocity, in the spin plane. Fl = Cl × q × A / m × (ω̂ × v̂)

Where q = ½ρv² (dynamic pressure), A = cross-sectional area, m = ball mass.

## Cd/Cl Model

The aerodynamic coefficients use a sigmoid transition at the dimpled-ball drag
crisis (~Re 70,000), with saturating spin-dependent drag and exponential lift
saturation.

### Drag coefficient

```
Cd = Cd_base(Re) + Cd_spin × SR / (1 + SR / SR_sat)

where:
  Cd_base = Cd_sub + (Cd_super - Cd_sub) × σ(Re)
  σ(Re)   = 1 / (1 + exp(-(Re - Re_crit) / Re_width))    [sigmoid]
  SR      = ω × R / V                                      [spin ratio]
  Re      = ρ × V × d / μ                                  [Reynolds number]
```

The base drag transitions from subcritical to supercritical at the drag crisis.
Spin adds drag through a saturating term: near-linear (slope Cd_spin) below
SR ≈ SR_sat, capped at Cd_spin × SR_sat at high spin ratio. Spin drag is
dominated by induced drag from Magnus lift, and since Cl saturates in SR, the
drag it induces saturates too. (Set SR_sat = ∞ for a purely linear term, as in
the PATENT_1997 model.)

### Lift coefficient

```
Cl = Cl_max(Re) × (1 - exp(-SR / SR_scale))

where:
  Cl_max = Cl_sub + (Cl_super - Cl_sub) × σ(Re)
```

Lift follows an exponential saturation curve in spin ratio. This means:
- **Low SR (driver, ~0.07):** Cl is well below the ceiling (~0.15)
- **High SR (iron, ~0.29):** Cl approaches the ceiling (~0.30)
- **Very high SR (wedge, ~0.45):** Cl is fully saturated (~0.32)

This spin-sensitive lift is a key feature of modern ball design — it rewards
the "high launch, low spin" fitting strategy that dominates modern equipment.

## Ball Models

### `BallModel::TOUR` — Modern Tour Ball (2020s)

| Parameter | Value | Role |
|---|---|---|
| mass_kg | 0.04593 | Ball mass (USGA max: 1.62 oz) |
| diameter_m | 0.04267 | Ball diameter (USGA min: 1.68 in) |
| inertia_factor | 0.4 | I = k·m·R², uniform solid sphere (ground phase only) |
| cd_sub | 0.295 | Subcritical base drag |
| cd_super | 0.186 | Supercritical base drag (zero spin) |
| cd_spin | 0.85 | Spin-drag slope (saturating) |
| sr_sat | 0.25 | Spin-drag saturation constant (cap ≈ 0.21) |
| cl_sub | 0.203 | Subcritical lift ceiling |
| cl_super | 0.324 | Supercritical lift ceiling |
| sr_scale | 0.115 | Lift saturation rate |
| re_crit | 69,600 | Drag crisis Reynolds number |
| re_width | 106 | Sigmoid transition width |

Aerodynamic sources:
- Li, Tsubokura & Tsunoda (2017), *Flow Turb. Combust.* 99(3) — LES CFD, Cd/Cl at Re=1.1×10⁵
- Lyu, Kensrud, Smith & Tosaya (2018), *ISEA Proc.* — trajectory fit on 13 production balls
- Crabill (2019), *Sports Eng.* 22:1-9 — high-order CFD at Re=1.5×10⁵
- Bearman & Harvey (1976), *Aero. Q.* 27:112-122 — original wind tunnel study
- Trackman PGA Tour Averages (2024) — carry and descent angle validation

### Calibration

The TOUR constants are jointly calibrated against two references:

1. **2024 Trackman PGA Tour Averages** — 12 clubs, carry + descent angle,
   normalized to 25°C sea level. Real-world radar data; authoritative for
   descent angles.
2. **A 1,178-shot per-shot validation set** spanning 26–153 mph ball speed
   and 0.8k–14k rpm spin, fitting carry, peak height, descent angle, and
   offline for every shot — providing breadth (amateur speeds, sidespin and
   offline behavior) that tour averages can't.

The saturating spin-drag term is what lets both references fit at once:
steep iron descent angles (47–52°) come from heavy late-flight induced drag
on high-spin balls, while subcritical lift stays at a wind-tunnel-plausible
cl_sub=0.203 — which in turn produces correct sidespin curvature and
altitude response.

### Fit quality

Trackman 2024 (25°C sea level): carry RMSE **2.8 yds**, DA RMSE **1.3°**,
9 of 12 clubs within ±3y (3-wood, hybrid and PW fall outside).
Per-shot validation set (1,178 shots): carry RMSE **3.2 yds** (bias
−1.0 yds), peak **2.6 ft**, descent **2.7°**, offline RMSE **1.2 yds**
(bias −0.3 yds).

Per-club results:

| Club | Model | Target | Δ | DA° (model) | DA° (TM) | DA Δ |
|---|---:|---:|---:|---:|---:|---:|
| Driver | 279.6y | 282.0y | -2.4y | 39.7 | 39.0 | +0.7 |
| 3-Wood | 253.4y | 249.0y | +4.4y | 46.6 | 44.0 | +2.6 |
| 5-Wood | 238.7y | 236.0y | +2.7y | 48.6 | 48.0 | +0.6 |
| Hybrid | 226.0y | 231.0y | -5.0y | 48.3 | 49.0 | -0.7 |
| 3-Iron | 220.6y | 218.0y | +2.6y | 46.7 | 48.0 | -1.3 |
| 4-Iron | 210.3y | 209.0y | +1.3y | 46.9 | 49.0 | -2.1 |
| 5-Iron | 199.8y | 199.0y | +0.8y | 47.8 | 50.0 | -2.2 |
| 6-Iron | 188.1y | 188.0y | +0.1y | 49.9 | 50.0 | -0.1 |
| 7-Iron | 173.9y | 176.0y | -2.1y | 50.5 | 51.0 | -0.5 |
| 8-Iron | 163.5y | 164.0y | -0.5y | 50.9 | 51.0 | -0.1 |
| 9-Iron | 152.1y | 152.0y | +0.1y | 51.2 | 52.0 | -0.8 |
| PW | 137.3y | 142.0y | -4.7y | 52.2 | 52.0 | +0.2 |

Environmental response (density-only): Denver +6.9% driver carry (published
~6.1%), 90°F-vs-50°F delta ~9y (published ~8y).

### Known limitations

1. **PW and hybrid carry ~5y short, 3-wood ~4y long.** Residual structural
   limit of the single-sigmoid model.
2. **Fast low-spin shots (140+ mph, <3k rpm) run a few yards short** with
   slightly high peaks relative to the per-shot validation set.
3. **Ball COR vs temperature unmodeled** — only air density responds to
   temperature; cold-ball distance loss is understated.
4. **Single ball model.** Real balls vary ±5% in drag across dimple designs.
5. **Spin on the ground: modeled, lightly validated.** The full spin vector
   is carried through every impact, slide and roll (see
   [Ground Model](#ground-model)), so a curving ball kicks further off-line
   at the first bounces and spin along the direction of travel drifts the
   ball while it slides. Gaps: no published measurement of these lateral
   effects to validate against; the pivot (torsional friction) lever is
   uncalibrated; the fairway's compliance is chosen, not measured; drag is
   ignored between hops; the ground is level (no slopes).

### `BallModel::PATENT_1997` — 1990s Tour Ball

| Parameter | Value | Role |
|---|---|---|
| mass_kg | 0.04593 | Same physical ball (USGA limits unchanged) |
| diameter_m | 0.04267 | Same |
| inertia_factor | 0.4 | Same |
| cd_sub | 0.19 | Lower subcritical base |
| cd_super | 0.29 | Much higher supercritical base drag |
| cd_spin | 0.15 | Spin-dependent drag |
| sr_sat | ∞ | Linear spin drag (no saturation in patent data) |
| cl_sub | 0.12 | Subcritical lift |
| cl_super | 0.285 | Similar lift ceiling |
| sr_scale | 0.01 | **Near-instant saturation** — Cl independent of SR |
| re_crit | 100,000 | Higher drag crisis point |
| re_width | 8,000 | Gradual transition |

Source: US6186002B1 (Quintavalla/Acushnet, filed 1997) — force-balance wind
tunnel measurements on production golf balls (Table 1, 7 data points). Model
fitted with RMSE: Cd=0.008, Cl=0.006.

Key differences from modern balls:
- **Base Cd much higher** (0.29 vs 0.186 + spin term; effective driver Cd
  0.30 vs 0.24) — dominant factor for carry gap
- **Cl saturates instantly** (sr_scale=0.01) — lift is the same for driver and
  wedge spin rates, unlike modern balls where Cl is spin-sensitive
- **Lower subcritical Cl** (0.12 vs 0.203) and a gradual drag crisis
  (re_width=8000 vs 106)

## Environment Model

Air density is computed from altitude and temperature using the barometric
formula for the troposphere:

```
ρ = ρ₀ × (1 - L×h/T₀)^(g×M/(R×L)) × (T_standard / T_actual)
```

Where ρ₀=1.225 kg/m³, L=0.0065 K/m, T₀=288.15 K. The temperature correction
scales by the ratio of standard to actual temperature at altitude (ideal gas
law, ρ ∝ 1/T at constant pressure).

Dynamic viscosity of air is constant at 1.81×10⁻⁵ Pa·s (adequate for 0–40°C).

## Spin Decay

Spin decays exponentially during flight:

```
ω(t) = ω₀ × exp(-λ × t)
```

Default λ = 0.04 s⁻¹ (~21% loss over 6 seconds). Published data:
- Trackman (2010): ~4%/s → λ=0.04 → 78% remaining at 6s
- Lyu et al. (2018): 4%/s trajectory fits
- Tutelman: τ=30s → λ=0.033 → 82% at 6s

## Ground Model

Everything after landing is a rigid sphere with point contact, carrying the
full spin vector ω through every contact. The tracked point is the ball's
centre, and the ball touches the ground when its centre reaches the ground
height (the radius is not added), as in the flight's landing; on ground
level with the tee a ball touching the ground has `y = 0`. The ground is
level in this version; slopes are not yet modeled.

### What spin can and cannot do on the ground

The contact point is at `r = −R n` (n the ground normal, R the radius) and
slips at `u = v_t + ω × r`; friction acts against `u`.

- **Spin about the normal can't move the contact point**, so it can't kick
  the ball; pivot (torsional) friction only decays it. On level ground this
  is the launch-monitor *sidespin* component, ω_y. Sidespin therefore does
  not steer a sliding or rolling ball.
- **A curving ball kicks further off-line at its first bounces.** The spin
  axis is fixed in flight, perpendicular to the start line, but a curving
  ball lands heading further off-line than that. The backspin then slips
  the contact point partly across the heading, and friction removes more
  velocity along the start line than across it: a fade kicks right, a draw
  left. In the grip regime on rigid ground the heading after the bounce is
  exactly `atan((5/7)·v_z / ((5/7)·v_x − (2/7)·R·ω_b))` (backspin ω_b about
  +Z). The effect is **modest**: the driver preset with ±600 rpm sidespin
  lands heading 12.4° off-line and leaves its first bounce at 13.9°.
- **Spin along the direction of travel** ("rifle spin") slips the contact
  point sideways, so it drifts the ball laterally while it slides.
- **Penner's tilt adds a small throw the other way.** Each impact is against
  a tilted effective plane (below), so a fraction sin β of ω_y lies in the
  contact plane and slips the contact point sideways at impact: a ball
  spinning clockwise from above (fade sidespin) is thrown left, like side
  spin off a cushion. Pivot friction then clears ω_y, so the throw acts
  mainly at the landing impact. For the ±600 rpm driver on `FAIRWAY`,
  measured as lateral run across the landing heading, the heading-mismatch
  kick alone (the same landing with ω_y removed) adds about 1.5 m, the
  throw takes back about 0.1 m, and the net is about +1.4 m (1.5 yd) on a
  32 m run. The fade comes to rest 35.5 m (38.8 yd) right of the target
  line, 8.3 m right of its landing point.

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
β    = min(compliance · θc, φ)
n'   = n cos β − t̂ sin β                          (effective normal)
e    = restitution at |v·n'|                       (Penner eq. 5)
j_n  = (1 + e) · |v·n'|                            (normal impulse / m)
u    = (v − (v·n')n') + ω × (−R n')               (contact slip)
j_t  = min(|u| · k/(1+k), μ · j_n)                 (along −û)
v   += j_n n' − j_t û
ω   += (−R n') × (−j_t û) / (k R²)
Δω_n = −sign(ω_n) · min(|ω_n|, pivot · j_n / (k R²))   about n
```

`k` is the ball's inertia factor (`I = k m R²`, 2/5 for a uniform sphere,
`BallModel::inertia_factor`). Friction either brings the contact point to
rest (grip: the ball leaves rolling on the effective plane) or slides
throughout at the Coulomb limit. No impulse is applied when the ball is not
moving into the ground (grazing).

For a ball with no spin about the vertical travelling in a vertical plane,
this is exactly Penner's planar solution: the normal component
`v_n' = v_n cos β + v_t sin β` (eq. 6), Daish's rigid-surface rebound in
the effective frame (eqs. 2–4: grip `(5/7)·v_t' − (2/7)·R·ω_back`, or slip
with friction impulse `μ(1+e)·v_n'`), rotated back (eq. 7). The test suite
checks the vector form against the planar solution impact by impact.

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
it into slip. With the 1 mm lever a full-shot landing removes all of it in
one impact: a driver landing on `FAIRWAY` can remove up to about 155 rad/s
(~1500 rpm), and a ±600 rpm-sidespin drive lands with 47 rad/s (450 rpm).
So on level ground the path is the same for every lever from 1 mm up, and
on rigid ground (no tilt) for every lever including zero. With no pivot
friction at all the spin survives to the later impacts and their throws
move that drive's rest by about 0.1 m.

### Hops and the bounce-to-roll transition

As in Penner, eq. 8 applies to every impact, not just the first. After the
first bounce the ball makes a run of progressively smaller hops; each
shallow impact tilts the effective plane and turns forward speed into a
small rebound. For a drive these hops cover most of the run, and the ball
rolls only the last few metres (Penner's Fig. 4). Penner's run model, hops
plus roll, agrees fairly well with empirical drive runs (Cochran & Stobbs;
Williams; his Fig. 5), with the model's runs slightly longer.

The ball keeps bouncing while its rebound would rise at least 5 mm above the
ground, `(v·n)² / (2g) ≥ 5 mm` (`v·n` ≥ ~0.313 m/s), Penner's transition
height; he tested 1–10 mm and found only a minor effect on the run. Its
normal velocity is then dropped and it slides or rolls. A safety bound of
200 impacts guards degenerate surfaces (e.g. elastic and frictionless);
real shots end by the 5 mm rule well before it. Between impacts the ball
follows a gravity-only (drag-free) path, integrated with RK4 at 1 ms; the
next ground crossing is solved exactly within the final step.

### Slide and roll

On the ground the state is position, velocity and the spin vector,
integrated with RK4 at 1 ms within each smooth phase, with the surface
sampled under the contact point every step:

- **Slide** (the contact point slips): friction `−μ·m·g·û` and its torque
  `(−R n) × F`. The slip speed falls linearly at `((1+k)/k)·μ·g`
  (3.5·μ·g for k = 2/5), so the slide is stepped exactly to the moment the
  slip reaches zero. Backspin can reverse the ball.
- **Roll**: the spin tangential to the ground is locked to the velocity, and
  the ball decelerates at `μ_eff·g/(1+k)` with the speed-dependent rolling
  resistance of Kensrud et al. (2024):

  ```
  μ_eff = mu_roll · (1 + (v / roll_v_ref)²)
  ```

  At driver post-bounce speeds this is many times the base value; at
  stimpmeter speeds it is about 1.5×. Integrated, a ball that starts rolling
  at `v₀` stops after `(1+k)·v_ref² / (2·mu_roll·g) · ln(1 + v₀²/v_ref²)`:
  3.48 m (stimp 11.4 ft) on `GREEN` and 1.60 m (5.2 ft) on `FAIRWAY` for
  the stimpmeter's 1.83 m/s.
- **Spin about the normal** decays at `pivot·g / (k R²)` in both phases,
  saturating at zero.
- **Rest**: a rolling ball slower than one step's deceleration (or than
  `stop_speed`, default 0.01 m/s) finishes at constant deceleration and
  stops. The result reports why the ground phase ended: `Rest`, `MaxTime`
  (default 30 s after landing) or `MaxImpacts`.

### Surfaces

| Surface | Restitution | `mu_slide` | `compliance` | `mu_roll` | `roll_v_ref` | `pivot` |
|---|---|---|---|---|---|---|
| `Surface::FAIRWAY` | `Penner2002` | 0.40 | 0.9 | 0.12 | 2.5 m/s | 1 mm |
| `Surface::GREEN` | `Penner2002` | 0.40 | 1.0 | 0.055 | 2.5 m/s | 1 mm |

| Parameter | Source |
|---|---|
| eq. 5 restitution | Penner (2002), from Penner (2001), *Am. J. Phys.* 69, 563, measured on relatively soft turf |
| θc law (eq. 8), `compliance = 1.0` | Penner (2002): θc = 15.4° best fit to Haake's (1990) green impact at 18.6 m/s, 44.4° from vertical, 484 rad/s backspin; linear in v and φ as a first-order approximation |
| `mu_slide = 0.40` | Daish (1972), as used by Penner for a ball on a green |
| 5 mm transition height | Penner (2002), set above the ~2 mm equilibrium sinkage of a ball on a green (Hubbard & Alaways) |
| `mu_roll`, `roll_v_ref` | Kensrud et al. (2024), "Modeling golf ball roll with a variable coefficient of friction", *Sports Engineering*: speed-dependent form; `GREEN`'s `mu_roll` gives stimp 11.4 ft |
| `pivot` | Hertzian full-slip torsional couple (Johnson 1985) with an assumed 4 mm contact radius; uncalibrated |
| `FAIRWAY` `compliance = 0.9` | Chosen, not measured (see below) |

Penner uses one turf model throughout: eq. 8 unchanged, with his green
restitution and μ, both for drives compared against empirical drive runs
(Fig. 5) and for iron shots onto greens. He represents a firmer green by a
lower μ (0.25, Fig. 8) rather than by changing θc, noting that eq. 8 would in
general need adjusting for such greens. `FAIRWAY` keeps eq. 8 at a chosen
`compliance = 0.9`, slightly below the green's 1.0, so fairway rebounds
come off slightly firmer and run further. The run is highly sensitive to
this value. Fairway run (yards past carry) by `compliance`, other fields as
`Surface::FAIRWAY`:

| `compliance` | driver | 7-iron (`seven_iron`) | wedge |
|---|---|---|---|
| 1.00 | 25.6 | 3.8 | 0.1 |
| 0.90 (`FAIRWAY`) | 34.6 | 5.9 | 0.8 |
| 0.85 | 39.7 | 7.9 | 1.3 |
| 0.75 | 51.3 | 12.3 | 2.4 |
| 0.50 | 95.1 | 28.2 | 5.6 |

With `GREEN`, the impact reproduces Penner's published drive results: a
landing at 26.8 m/s, 68.0° from vertical, 224 rad/s backspin rises 2.82 m on
the first bounce, and 26.2 m/s at 46.2° with 315 rad/s rises 1.18 m (Penner
reports 2.82 m and 1.17 m). Those drives make 43 and 36 impacts before
rolling and stop 39.0 m and 20.5 m past the landing point (Penner's Fig. 4:
about 37 m and 19 m). A 9-iron landing at 23.8 m/s, 34.0° bounces backward
with 809 rad/s backspin, and its backward run grows with backspin at 1000
and 1200 rad/s (Penner's Fig. 7). With μ = 0.25 the same landings first
bounce forward still carrying backspin, then run back (Penner's Fig. 8).

### Plugging in your own ground

The ground API has three layers:

1. **`Surface`** — flat ground of one material, also a complete ground
   model: `simulate_shot(&input, &ball, &env, &Surface::FAIRWAY)`, or with
   a field changed, `&Surface { mu_roll: 0.09, ..Surface::FAIRWAY }`.
2. **`GroundPhysics<T: Terrain>`** — the physics above over any terrain.
   `Terrain::sample(x, z)` returns a `GroundSample { height, normal,
   surface }`; `FnTerrain(|x, z| ...)` wraps a closure, so per-point
   surfaces need no new types:

   ```rust
   let range = GroundPhysics::new(FnTerrain(|x: f64, z: f64| GroundSample {
       height: 0.0,
       normal: DVec3::Y,
       surface: if on_green(x, z) { Surface::GREEN } else { Surface::FAIRWAY },
   }));
   simulate_shot(&input, &ball, &env, &range);
   ```

   Its public fields are the tunables: `dt` (1 ms), `min_hop_height`
   (5 mm), `stop_speed` and `roll_threshold` (0.01 m/s), `max_impacts`
   (200), `max_time` (30 s).
3. **`GroundModel`** — replace the physics: implement `height(x, z)` and
   `settle(landing, ball)`, and optionally `landing(a, b)` (by default an
   interpolated height crossing; an engine can substitute its own sweep).
   It is object-safe, so `&dyn GroundModel` works. The building blocks
   `ground::impact` and `GroundPhysics::step` are public, and the result
   types (`GroundResult`, `Contact`, `GroundPoint`) have `new`
   constructors. `BallState` (time, position, velocity, spin) is a plain
   struct.

`simulate_shot` launches from `(0, height(0, 0), 0)`, ends the flight at
`landing()`, and hands the landing state to `settle()`. `simulate_flight`
always lands on the plane y = 0. Both interpolate the landing point to the
exact ground crossing.

### Limitations

- **Green data only.** The restitution fit and the θc law both come from
  green-like turf. `FAIRWAY` differs from `GREEN` by compliance and rolling
  resistance; its compliance value has no published source, and the run
  depends strongly on it.
- **One impact measurement.** θc rests on a single Haake impact and is
  extrapolated linearly in speed and incidence angle. Spin and turf
  firmness are not inputs.
- **Lateral ground effects are unvalidated.** The kick of a curving ball,
  the sidespin throw at impact and rifle-spin drift follow from the
  rigid-sphere contact model; no published measurement is compared against.
  The pivot lever is uncalibrated.
- **No drag between hops.** A drive's first hops leave at up to ~8 m/s and
  stay airborne for over a second, so drag is not negligible; ignoring it
  lengthens the run somewhat (an estimated few yards on a drive).
- **Level ground.** Every impact is on a level surface. `Terrain` must be
  level (normal +Y, the same height everywhere); its surfaces are honoured
  point by point. Slopes are not modeled. Penner attributes his model's
  slightly long drive runs to the same assumption.

## Implementation

The aerodynamic model lives in `src/aero.rs`, the ground model in
`src/ground.rs`. Key types:

- **`BallModel`** — physical dimensions, inertia factor, aero coefficients.
  Presets: `TOUR`, `PATENT_1997`.
- **`AeroParams`** — per-shot spin state (axis, rate, decay), from a spin
  vector (`AeroParams::from_spin_vector`, with
  `ShotInput::launch_spin_rads`) or straight-launch backspin/sidespin RPM
  (`AeroParams::from_spin`).
- **`Environment`** — altitude + temperature → air density.
- **`Surface`**, **`Terrain`**, **`GroundPhysics`**, **`GroundModel`** — the
  ground (above).

Flight simulation: `simulate_flight(input, ball, env) -> FlightResult`

RK4 integration at 1ms timestep, spin decays each step, trajectory sampled
every 10ms. Stops at the interpolated crossing of ground level.

Full shot: `simulate_shot(input, ball, env, ground) -> ShotResult` runs the
flight onto `ground`, then `ground.settle()`: impacts, hops, slide and roll
to rest (`ShotResult.ground`, with absolute times).

## Test Suite

Unit tests live beside the code in each module. Test categories:

- **Structural invariants** — energy conservation, origin start, axis normalization
- **Sign conventions** — azimuth, slice/hook, push-draw, mirror symmetry, rest position
- **Patent wind tunnel** — 7 Cd/Cl data points against PATENT_1997 model
- **Trackman carry** — 8 shots validated against published PGA Tour Averages (±10-12y)
- **Descent angle** — iron DA within expected ranges
- **Altitude/temperature** — Denver gain, hot-cold delta, directional checks
- **Spin decay** — exponential decay rate against Trackman measurements
- **Era comparison** — TOUR vs PATENT_1997 carry gaps, Cd reduction, Cl sensitivity
- **Launch spin axis** — a pushed shot with no sidespin stays on its start line; exact landing-spin projections
- **Ground: planar regression** — the 3D impact reproduces the planar Penner solution impact by impact (preset driver, 7-iron and wedge on both presets) to float tolerance
- **Ground: physics** — Penner (2002) restitution, drive bounce heights and runs (Fig. 4), 9-iron spin-back (Figs. 7–8), 5 mm bounce-to-roll transition, compliance tilt (2β rotation, cap), grip-regime kick formula, no kick from spin about the normal, tilt-throw magnitude (grip and slip), rifle-spin drift, grazing impacts, pivot friction about the true normal and its independence from the level-ground path, general inertia factor (k = 2/3), exact slide-to-roll time, stimpmeter closed form
- **Ground: invariants** — energy never increases at an impact or along the ground phase; mirror (±sidespin, ±azimuth) and rotation about Y give mirrored and rotated results; fades rest right, draws rest left
- **Ground: API** — per-point surfaces (fairway-to-green roll), `MaxTime`/`MaxImpacts` termination, tunables, interpolated landing, `&dyn GroundModel`

Run: `cargo test`
