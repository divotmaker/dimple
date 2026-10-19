# CLAUDE.md — dimple conventions

`dimple` is a pure-math golf ball-flight library: aerodynamics, RK4 flight,
the ground phase (impacts, hops, slide, roll), and shot statistics. Single
library crate at the repo root.

## Layout

| Path | Contents |
|---|---|
| `src/aero.rs` | `BallModel` (Cd/Cl model + presets), `AeroParams` (spin state), force model |
| `src/trajectory.rs` | `ShotInput`, `simulate_flight`, `simulate_shot`, `interpolate_trajectory` |
| `src/ground.rs` | `Surface`, `Terrain`, `GroundPhysics`, `GroundModel`: 3D impact (Penner), hops, slide/roll |
| `src/ground/` | Ground test module and planar regression reference data |
| `src/environment.rs` | Altitude + temperature → air density |
| `src/rk4.rs` | Fixed-step RK4 integrators (`step` for flight, generic `step_state`) |
| `src/stats.rs` | `ShotStats`, `SessionStats` |
| `src/units.rs` | Constants and unit conversions |
| `docs/AERO.md` | Public model documentation: formulas, parameters, sources, fit quality |

## Physics conventions

- **SI internally**: meters, seconds, kilograms, radians, rad/s. Golf units
  (mph, degrees, rpm, yards, feet) appear only in `ShotInput`, `ShotStats`,
  and `units`.
- **RK4 only**, fixed 1 ms timestep for flight and within each smooth ground
  phase (steps end exactly at phase switches). No adaptive or lower-order
  integrators.
- **Frame**: right-handed, origin at the tee. +X downrange, +Y up,
  +Z right of the target line from the golfer's view.
- **Lateral signs: positive = right** everywhere: `launch_azimuth_deg`,
  `sidespin_rpm`, `lateral_m`, `lateral_yards`, `landing_sidespin_rads`, and
  the `z` of every public vector. This matches the Flight Relay Protocol.
- **Spin**: one angular velocity vector, right-hand rule, from launch to
  rest: `ω = backspin · (−sin az, 0, cos az) − sidespin · ŷ`, i.e.
  `(0, -sidespin, backspin)` for a straight launch. Positive backspin lifts;
  positive sidespin (clockwise from above) curves right.
- Any change touching lateral or spin signs needs a test asserting the
  observable direction (see the `LATERAL SIGN CONVENTION` and
  `GROUND LATERAL SIGNS` tests in `trajectory.rs`).

## Scope

Ball-flight physics and shot statistics only. Game logic (target scoring,
range modes), rendering, I/O, and device protocols belong in consuming
applications.

## Rust conventions

- Rust 2024 edition, stable toolchain (`rust-toolchain.toml`).
- One runtime dependency (`glam`). Justify any addition.
- No `unwrap()` in library paths; `expect` only on internal invariants that
  cannot fail (e.g. a flight always has its launch point). `#[must_use]` on
  value-returning fns.
- Never modify the `version` field in `Cargo.toml`; release CI bumps it.
- Comments and docs describe current behaviour only — no history.

## Docs policy

Public docs explain how the software works: formulas, parameters, published
sources, fit quality, limitations. Calibration datasets, fitting tools, and
research notes are not part of this repository.

## Gates

All must pass before a change is done:

```bash
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
```

`cargo test` includes the crate-level and README doctests.

## Releases

Push a `vX.Y.Z` tag. CI runs the gates, publishes to crates.io, creates a
GitHub Release, then commits the next patch version to `main`.
