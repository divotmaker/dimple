//! The ground model against published measurements.
//!
//! - Biber, S. W., Jones, K. M., Champneys, A. R., Green, R. and Szalai, R.
//!   (2023), "Measurements and linearized models for golf ball bounce on a
//!   green", *Sports Engineering* 26:50, doi:10.1007/s12283-023-00442-4,
//!   CC BY 4.0. Bounces on an artificial mat on a rigid base (Campaign A)
//!   and on a natural tee (Campaign B). Used here: the fitted contact-angle
//!   coefficients k_P of Table 3, the "typical driver landing" of Fig. 3,
//!   and the piecewise-affine landing-to-lift-off fits of eqs (8) (A) and
//!   (9) (B), with their error measures, eqs (5)–(6). The paper publishes
//!   its per-bounce measurements as figures and fits, so the fits stand in
//!   for them here; they carry the paper's 6.6–7.8% fit error (Table 2).
//! - Brede, A. D. (1991), "Correction for slope in green speed measurement
//!   of golf course putting greens", *Agron. J.* 83:425–426,
//!   doi:10.2134/agronj1991.00021962008300020032x: up- and downhill
//!   stimpmeter readings a and b combine as 2ab/(a+b), which reads the same
//!   on every slope from 0 to 5.6%.
//! - Whitlark, B. and Pringle, M. (2012), "Can't Stop Your Ball on the
//!   Green?", *USGA Green Section Record* 50(25): landing shots on a 1.5%
//!   bermudagrass green, mowed to about 12 ft. Summary values only: mean
//!   bounce and roll 32.5 ft downhill and 13.5 ft uphill, firmer turf
//!   (0.08 in on the TruFirm) 9–12 ft further uphill, green speeds 11 ft 3 in
//!   and 12 ft 3 in no different, some wedge distances negative.
//! - Cross, R. (2024), "Trajectory of a golf ball on a sloping green", *Eur.
//!   J. Phys.* 45, 035002, doi:10.1088/1361-6404/ad2c2d: under constant rolling
//!   friction a ball rolled across an incline stops at a displacement angle
//!   independent of its launch speed, at a distance proportional to v₀².

use super::*;
use crate::environment::Environment;
use crate::units;

const BALL: BallModel = BallModel::TOUR;
const R: f64 = 0.04267 / 2.0;

/// No air: ballistic hops, gravity only.
const VACUUM: f64 = 0.0;

/// Air density (kg/m³) at sea level, as `simulate_shot` passes it.
fn air() -> f64 {
    Environment::SEA_LEVEL.air_density()
}

/// A plane rising at `grade` (rise over run) toward +X, made of `surface`.
fn incline(
    grade: f64,
    surface: Surface,
) -> GroundPhysics<FnTerrain<impl Fn(f64, f64) -> GroundSample + Copy>> {
    let normal = DVec3::new(-grade, 1.0, 0.0).normalize();
    GroundPhysics::new(FnTerrain(move |x: f64, _z: f64| GroundSample {
        height: grade * x,
        normal,
        surface,
    }))
}

/// A ball at the origin of a plane rising at `grade`, rolling without slip
/// at `speed` in the direction `heading` (horizontal, unit), which the
/// velocity follows up or down the plane.
fn rolling_start(grade: f64, heading: DVec3, speed: f64) -> BallState {
    let normal = DVec3::new(-grade, 1.0, 0.0).normalize();
    let along = DVec3::new(heading.x, grade * heading.x, heading.z).normalize();
    let vel = along * speed;
    BallState {
        time: 0.0,
        pos: DVec3::ZERO,
        vel,
        spin: normal.cross(vel) / R,
    }
}

/// A stimpmeter roll along +X on a plane rising at `grade` (negative:
/// downhill).
fn stimp_roll(surface: Surface, grade: f64) -> GroundResult {
    incline(grade, surface).settle(
        &rolling_start(grade, DVec3::X, STIMPMETER_SPEED),
        &BALL,
        VACUUM,
    )
}

/// Length (ft) of a stimpmeter roll along a plane rising at `grade`,
/// measured along the surface as a tape on the green would.
fn stimp_roll_ft(surface: Surface, grade: f64) -> f64 {
    let g = stimp_roll(surface, grade);
    assert_eq!(g.termination, Termination::Rest, "grade {grade}");
    units::meters_to_feet(g.rest.pos.length())
}

/// Brede's slope-corrected reading 2ab/(a+b).
fn brede(up: f64, down: f64) -> f64 {
    2.0 * up * down / (up + down)
}

/// `surface` with `roll_law`, re-pinned to the stimp it reads now.
fn same_stimp(surface: Surface, roll_law: RollLaw) -> Surface {
    Surface {
        roll_law,
        ..surface
    }
    .with_stimp(surface.stimp_feet())
}

#[test]
fn brede_harmonic_mean_exact_under_constant_law() {
    // Constant rolling resistance decelerates at g(μ cos θ ± sin θ)/(1+k) up
    // and down the slope, so 2ab/(a+b) = (1+k)v₀²/(2gμ cos θ): the level
    // reading over cos θ, the tape measuring along the surface.
    for preset in [Surface::GREEN, Surface::FAIRWAY] {
        let surface = same_stimp(preset, RollLaw::Constant);
        let level = surface.stimp_feet();
        for pct in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            let grade: f64 = pct / 100.0;
            let up = stimp_roll_ft(surface, grade);
            let down = stimp_roll_ft(surface, -grade);
            assert!(up < level && level < down, "{pct}%: {up} {level} {down}");
            let cos = 1.0 / grade.hypot(1.0);
            let corrected = brede(up, down) * cos;
            assert!(
                (corrected / level - 1.0).abs() < 1e-6,
                "{pct}%: 2ab/(a+b)·cos θ {corrected} vs level {level}"
            );
        }
    }
}

#[test]
fn brede_harmonic_mean_under_preset_law() {
    // Under the presets' quadratic law the up- and downhill rolls cover
    // different speeds, so the correction is approximate. GREEN at 1–4%:
    // up/down 9.95/13.41, 8.82/16.26, 7.93/20.73, 7.20/28.94 ft against a
    // level 11.42 ft, i.e. +0.04/+0.17/+0.43/+0.95%; FAIRWAY at most +0.5%
    // to 6%.
    for (preset, max_pct, tolerance) in [
        (Surface::GREEN, 3.0, 0.005),
        (Surface::GREEN, 4.0, 0.015),
        (Surface::FAIRWAY, 6.0, 0.006),
    ] {
        let level = preset.stimp_feet();
        for pct in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            if pct > max_pct {
                break;
            }
            let grade = pct / 100.0;
            let (up, down) = (stimp_roll_ft(preset, grade), stimp_roll_ft(preset, -grade));
            let error = brede(up, down) / level - 1.0;
            assert!(
                error.abs() < tolerance,
                "{pct}%: 2ab/(a+b) {} vs level {level}",
                brede(up, down)
            );
        }
    }
    // The plain average of the two readings, which ignores the slope, is
    // already about 10% long at 2% on GREEN; Brede's correction removes it.
    let (up, down) = (
        stimp_roll_ft(Surface::GREEN, 0.02),
        stimp_roll_ft(Surface::GREEN, -0.02),
    );
    assert!(0.5 * (up + down) > 1.05 * Surface::GREEN.stimp_feet());
}

#[test]
fn green_cannot_hold_downhill_stimp_roll_past_hold_grade() {
    // GREEN holds a ball at rest only to a 5.5% grade (tan θ = mu_roll),
    // where Brede held stimpmeter rolls to 5.6% on his carpet runway and
    // read field greens to 6%.
    let sample = |grade: f64| GroundSample {
        height: 0.0,
        normal: DVec3::new(grade, 1.0, 0.0).normalize(),
        surface: Surface::GREEN,
    };
    assert!(sample(0.054).holds_at_rest());
    assert!(!sample(0.056).holds_at_rest());

    // At 5% the downhill roll stops, about 51 ft out (4.5 × level).
    let down5 = stimp_roll_ft(Surface::GREEN, -0.05);
    assert!((40.0..65.0).contains(&down5), "{down5} ft");

    // At 6% it never does: it settles toward the quadratic law's terminal
    // speed, μ_eff(v) = tan θ at 0.75 m/s, and is still rolling when the
    // time limit ends it.
    let six = GroundPhysics {
        max_time: 20.0,
        ..incline(-0.06, Surface::GREEN)
    }
    .settle(
        &rolling_start(-0.06, DVec3::X, STIMPMETER_SPEED),
        &BALL,
        VACUUM,
    );
    assert_eq!(six.termination, Termination::MaxTime);
    assert!(six.rest.vel.length() > 0.7, "{}", six.rest.vel.length());

    // The constant law at the same stimp has the higher hold grade
    // (mu_roll 0.0687, 6.9%) and stops the same roll, about 90 ft out.
    let constant = same_stimp(Surface::GREEN, RollLaw::Constant);
    assert!(constant.mu_roll > 0.06);
    assert!(stimp_roll_ft(constant, -0.06) > down5);
}

#[test]
fn penner_contact_angle_in_scaled_units() {
    // Penner's θc = 15.4°·(v/18.6 m/s)·(φ/44.4°) is Biber et al.'s
    // β = k_P·s·φ with s = v/R in ball radii per second and φ in radians,
    // k_P = (15.4°/44.4°)·R/(18.6 m/s) for compliance 1.0.
    let k_p = PENNER_THETA_C / PENNER_REF_ANGLE * R / PENNER_REF_SPEED;
    assert!((k_p - 3.978e-4).abs() < 1e-7, "{k_p}");
    for (v_t, v_n) in [(10.0, 15.0), (20.0, 8.0), (2.0, 5.0), (25.0, 25.0)] {
        let s = f64::hypot(v_t, v_n) / R;
        let phi = f64::atan2(v_t, v_n);
        let beta = Surface::GREEN.contact_angle(v_t, v_n);
        assert!((beta / (k_p * s * phi) - 1.0).abs() < 1e-9);
    }

    // Their Table 3 fits are compliance 0.388 (Campaign A, the mat) and
    // 0.747 (Campaign B, the tee): `compliance` scales k_P.
    for (fitted_k_p, compliance) in [(1.543e-4, 0.388), (2.973e-4, 0.747)] {
        assert!((fitted_k_p / k_p - compliance).abs() < 1e-3);
        let surface = Surface {
            compliance,
            ..Surface::GREEN
        };
        let (v_t, v_n) = (12.0, 14.0);
        let s = f64::hypot(v_t, v_n) / R;
        let beta = surface.contact_angle(v_t, v_n);
        let theirs = fitted_k_p * s * f64::atan2(v_t, v_n);
        assert!((beta / theirs - 1.0).abs() < 2e-3, "{beta} vs {theirs}");
    }
}

/// One branch of a piecewise-affine bounce fit, p_F = A·p₀ + b, in Biber et
/// al.'s scaled units: p = (ẋ, ẏ, ω), velocities in ball radii per second,
/// ω in rad/s with backspin positive.
struct Affine {
    a: [[f64; 3]; 3],
    b: [f64; 3],
}

/// A campaign's fit: `roll` where `plane`·p₀ < `switch` (the ball enters
/// rolling during the bounce), `slip` above (it slips throughout).
struct BounceFit {
    plane: [f64; 3],
    switch: f64,
    roll: Affine,
    slip: Affine,
}

impl BounceFit {
    /// The fitted lift-off state for landing `p0`, and whether it is on
    /// the roll branch.
    fn lift_off(&self, p0: [f64; 3]) -> ([f64; 3], bool) {
        let rolls = dot(self.plane, p0) < self.switch;
        let f = if rolls { &self.roll } else { &self.slip };
        let mut out = f.b;
        for (o, row) in out.iter_mut().zip(f.a) {
            *o += dot(row, p0);
        }
        (out, rolls)
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Biber et al. (2023) eq. (8): Campaign A, artificial teeing turf 30 mm
/// thick on 38 mm of wood.
// The published coefficient 0.318 is a fitted value, not 1/π.
#[allow(clippy::approx_constant)]
const MAT: BounceFit = BounceFit {
    plane: [0.891, 0.763, 1.0],
    switch: 600.0,
    roll: Affine {
        a: [
            [0.536, -0.004, -0.205],
            [0.144, -0.407, -0.027],
            [-0.524, 0.091, 0.318],
        ],
        b: [-111.0, 23.7, 39.4],
    },
    slip: Affine {
        a: [
            [0.674, 0.107, -0.005],
            [0.103, -0.390, 0.017],
            [-0.153, 0.277, 1.01],
        ],
        b: [-154.0, 25.0, -267.0],
    },
};

/// Biber et al. (2023) eq. (9): Campaign B, a natural teeing area.
const TEE: BounceFit = BounceFit {
    plane: [0.377, 2.25, 1.0],
    switch: -485.0,
    roll: Affine {
        a: [
            [0.326, 0.222, -0.213],
            [0.185, -0.008, -0.076],
            [-0.333, -0.147, 0.224],
        ],
        b: [175.0, 90.4, -171.0],
    },
    slip: Affine {
        a: [
            [0.606, 0.649, -0.007],
            [0.146, -0.089, -0.006],
            [-0.192, 1.01, 0.887],
        ],
        b: [117.0, 49.3, 35.0],
    },
};

/// A landing at `speed` m/s, `descent_deg` below horizontal, moving +X with
/// `rpm` of backspin (negative: topspin).
fn landing(speed: f64, descent_deg: f64, rpm: f64) -> BallState {
    let a = units::deg_to_rad(descent_deg);
    BallState {
        time: 0.0,
        pos: DVec3::ZERO,
        vel: DVec3::new(speed * a.cos(), -speed * a.sin(), 0.0),
        spin: DVec3::new(0.0, 0.0, units::rpm_to_rads(rpm)),
    }
}

/// The scaled state p = (ẋ, ẏ, ω) of a ball moving in the X–Y plane.
fn scaled(s: &BallState) -> [f64; 3] {
    [s.vel.x / R, s.vel.y / R, s.spin.z]
}

/// Every landing of a grid of speeds (m/s), descent angles (degrees) and
/// spins (rpm), each inside the campaign's measured ranges (Table 1).
fn grid(speeds: &[f64], angles: &[f64], spins: &[f64]) -> Vec<BallState> {
    let mut out = Vec::new();
    for &v in speeds {
        for &a in angles {
            for &rpm in spins {
                out.push(landing(v, a, rpm));
            }
        }
    }
    out
}

/// Campaign A grid: 96 landings at 10–40 m/s, 25–70°, −2000–8000 rpm.
fn mat_grid() -> Vec<BallState> {
    grid(
        &[10.0, 20.0, 30.0, 40.0],
        &[25.0, 40.0, 55.0, 70.0],
        &[-2000.0, 0.0, 2000.0, 4000.0, 6000.0, 8000.0],
    )
}

/// Campaign B grid: 120 landings at 5–30 m/s, 20–80°, −2000–12000 rpm.
fn tee_grid() -> Vec<BallState> {
    grid(
        &[5.0, 10.0, 20.0, 30.0],
        &[20.0, 35.0, 50.0, 65.0, 80.0],
        &[-2000.0, 0.0, 3000.0, 6000.0, 9000.0, 12000.0],
    )
}

/// One landing: its scaled state, the fit's lift-off and branch, and
/// dimple's lift-off.
struct Bounce {
    p0: [f64; 3],
    fit: [f64; 3],
    rolls: bool,
    model: [f64; 3],
}

fn bounces(fit: &BounceFit, landings: &[BallState], surface: Surface) -> Vec<Bounce> {
    landings
        .iter()
        .map(|l| {
            let p0 = scaled(l);
            let (fitted, rolls) = fit.lift_off(p0);
            let out = impact(l, &surface.sample(0.0, 0.0), &BALL);
            Bounce {
                p0,
                fit: fitted,
                rolls,
                model: scaled(&out),
            }
        })
        .collect()
}

/// Biber et al.'s error measures, eqs (5)–(6), of dimple against the fit:
/// δ = mean of |p_F − P_F| / |p₀|, and per component Δq = mean of
/// |q_F − Q_F| / |p₀|. Returns (δ, Δẋ, Δẏ, Δω).
fn errors(b: &[Bounce]) -> (f64, f64, f64, f64) {
    let n = b.len() as f64;
    let mut e = (0.0, 0.0, 0.0, 0.0);
    for x in b {
        let d = [
            x.fit[0] - x.model[0],
            x.fit[1] - x.model[1],
            x.fit[2] - x.model[2],
        ];
        let p = dot(x.p0, x.p0).sqrt();
        e.0 += dot(d, d).sqrt() / p;
        e.1 += d[0].abs() / p;
        e.2 += d[1].abs() / p;
        e.3 += d[2].abs() / p;
    }
    (e.0 / n, e.1 / n, e.2 / n, e.3 / n)
}

/// Mean normal rebound ratio |ẏ_F / ẏ₀| of the fit and of dimple.
fn rebound(b: &[Bounce]) -> (f64, f64) {
    let n = b.len() as f64;
    let fit = b.iter().map(|x| x.fit[1] / -x.p0[1]).sum::<f64>() / n;
    let model = b.iter().map(|x| x.model[1] / -x.p0[1]).sum::<f64>() / n;
    (fit, model)
}

#[test]
fn biber_tee_turf_fit_agreement() {
    // Natural turf (Campaign B). dimple's GREEN against the fit: δ 12.2%
    // (Δẋ 7.4, Δẏ 3.0, Δω 7.7%); the paper's own physics fits reach 19–35%
    // against the data (Table 3), its affine fit 7.8% (Table 2).
    let landings = tee_grid();
    let green = bounces(&TEE, &landings, Surface::GREEN);
    let (delta, _, dy, _) = errors(&green);
    assert!(delta < 0.15, "δ {delta}");
    assert!(dy < 0.05, "Δẏ {dy}");
    let (fairway_delta, ..) = errors(&bounces(&TEE, &landings, Surface::FAIRWAY));
    assert!(fairway_delta < 0.15, "FAIRWAY δ {fairway_delta}");

    // Normal rebound: 0.388 against the fit's 0.408.
    let (fit_rebound, model_rebound) = rebound(&green);
    assert!(
        (fit_rebound - model_rebound).abs() < 0.05,
        "{model_rebound} vs {fit_rebound}"
    );

    // Spin reversal: of the 80 backspin landings, 26 leave with topspin in
    // dimple and 23 in the fit.
    let backspin: Vec<_> = green.iter().filter(|b| b.p0[2] > 0.0).collect();
    let reversed = |q: fn(&Bounce) -> f64| backspin.iter().filter(|b| q(b) < 0.0).count();
    let (fit_rev, model_rev) = (reversed(|b| b.fit[2]), reversed(|b| b.model[2]));
    assert!(fit_rev.abs_diff(model_rev) <= 6, "{model_rev} vs {fit_rev}");

    // Balls that enter rolling lift off slipping with excess topspin,
    // contact-point velocity v_P = ẋ + ω < 0 (their Fig. 6). dimple's
    // tilted contact plane gives the same sign in the ground frame for 56
    // of the 57 roll-branch landings; the exception slips throughout and
    // keeps its backspin.
    let roll: Vec<_> = green.iter().filter(|b| b.rolls).collect();
    let topspin_slip = roll
        .iter()
        .filter(|b| b.model[0] + b.model[2] < 0.0)
        .count();
    assert!(roll.len() > 40);
    assert!(
        topspin_slip * 10 >= roll.len() * 9,
        "{topspin_slip} of {}",
        roll.len()
    );

    // A rigid plane (no tilt, 20.6%) and the fits' μ ≈ 1 (22.2%) both do
    // worse than GREEN's μ 0.4 with Penner's tilt.
    for worse in [
        Surface {
            compliance: 0.0,
            ..Surface::GREEN
        },
        Surface {
            mu_slide: 1.0,
            ..Surface::GREEN
        },
    ] {
        let (d, ..) = errors(&bounces(&TEE, &landings, worse));
        assert!(d > delta + 0.05, "{worse:?}: δ {d} vs {delta}");
    }
}

#[test]
fn biber_artificial_turf_needs_higher_restitution() {
    // A mat on a rigid base (Campaign A) rebounds far more than Penner's
    // turf: the fit's mean normal rebound ratio is 0.58, GREEN's 0.32
    // (Table 3 fits r = 0.42–0.54 for the mat). So no turf preset is a mat.
    let landings = mat_grid();
    let green = bounces(&MAT, &landings, Surface::GREEN);
    let (fit_rebound, model_rebound) = rebound(&green);
    assert!(fit_rebound > 0.55, "{fit_rebound}");
    assert!(model_rebound < 0.4, "{model_rebound}");

    // A constant restitution with a lower contact angle and friction fits
    // it better (δ 30.5% → 17.2%, Δẏ 18.1% → 4.6%).
    let mat = Surface {
        restitution: Restitution::Constant(0.45),
        compliance: 0.4,
        mu_slide: 0.25,
        ..Surface::GREEN
    };
    let (green_delta, _, green_dy, _) = errors(&green);
    let (mat_delta, _, mat_dy, _) = errors(&bounces(&MAT, &landings, mat));
    assert!(mat_dy < 0.5 * green_dy, "Δẏ {mat_dy} vs {green_dy}");
    assert!(
        mat_delta < green_delta - 0.05,
        "δ {mat_delta} vs {green_delta}"
    );
}

#[test]
fn biber_driver_landing_rebound() {
    // Their typical driver landing (Fig. 3): 93.6 ft/s, 37.3°, 34.9 rev/s
    // of backspin. On GREEN the normal rebound is 0.354 against the tee
    // fit's 0.342, the forward ratio 0.214 against 0.277, and the ball
    // leaves with topspin, as in the fit.
    let l = landing(units::feet_to_meters(93.6), 37.3, 34.9 * 60.0);
    let p0 = scaled(&l);
    let (fit, rolls) = TEE.lift_off(p0);
    assert!(rolls);
    let model = scaled(&impact(&l, &Surface::GREEN.sample(0.0, 0.0), &BALL));
    let normal = |p: [f64; 3]| p[1] / -p0[1];
    let forward = |p: [f64; 3]| p[0] / p0[0];
    assert!((normal(model) - normal(fit)).abs() < 0.05);
    assert!((forward(model) - forward(fit)).abs() < 0.10);
    assert!(fit[2] < 0.0 && model[2] < 0.0);
}

/// Whitlark & Pringle's (2012) green: `GREEN` mowed to stimp 12 ft.
fn whitlark_green() -> Surface {
    Surface::GREEN.with_stimp(12.0)
}

/// Their fairway 7-iron landing (Table 1): 55 mph, 46°, 5100 rpm.
const SEVEN_IRON: (f64, f64, f64) = (55.0, 46.0, 5100.0);
/// Their fairway sand wedge landing: 49 mph, 51°, 9100 rpm.
const SAND_WEDGE: (f64, f64, f64) = (49.0, 51.0, 9100.0);

/// Bounce and roll (ft) from landing to rest along a plane rising at
/// `grade` in the direction of play, through sea-level air; negative when
/// the ball finishes behind its landing point.
fn bounce_and_roll_ft(surface: Surface, grade: f64, (mph, descent, rpm): (f64, f64, f64)) -> f64 {
    let g = incline(grade, surface).settle(
        &landing(units::mph_to_mps(mph), descent, rpm),
        &BALL,
        air(),
    );
    assert_eq!(g.termination, Termination::Rest);
    units::meters_to_feet(g.rest.pos.x * grade.hypot(1.0))
}

#[test]
fn whitlark_slope_asymmetry_on_a_green() {
    // 7-iron on a 1.5% grade: uphill 26.3 ft, level 33.3 ft, downhill
    // 42.5 ft. Whitlark & Pringle measured mean totals of 13.5 ft uphill
    // and 32.5 ft downhill over their shot mix, a 19 ft difference.
    let green = whitlark_green();
    let up = bounce_and_roll_ft(green, 0.015, SEVEN_IRON);
    let level = bounce_and_roll_ft(green, 0.0, SEVEN_IRON);
    let down = bounce_and_roll_ft(green, -0.015, SEVEN_IRON);
    assert!(
        0.0 < up && up < level && level < down,
        "{up} {level} {down}"
    );
    assert!((10.0..30.0).contains(&(down - up)), "{}", down - up);

    // The high-spin sand wedge stops by its pitch mark on every slope, a
    // little behind it (1.8 ft uphill, 0.8 ft level, 0.3 ft downhill); "in
    // some cases this distance was negative".
    for grade in [0.015, 0.0, -0.015] {
        let wedge = bounce_and_roll_ft(green, grade, SAND_WEDGE);
        assert!(wedge.abs() < 5.0, "grade {grade}: {wedge} ft");
    }
}

#[test]
fn whitlark_green_speed_barely_moves_bounce_and_roll() {
    // Stimp 11 ft 3 in against 12 ft 3 in made no difference to their
    // bounce and roll; dimple's uphill 7-iron moves 25.5 → 26.6 ft.
    let slow = bounce_and_roll_ft(Surface::GREEN.with_stimp(11.25), 0.015, SEVEN_IRON);
    let fast = bounce_and_roll_ft(Surface::GREEN.with_stimp(12.25), 0.015, SEVEN_IRON);
    assert!((fast - slow).abs() < 2.0, "{slow} vs {fast}");
}

#[test]
fn whitlark_compliance_is_the_firmness_knob() {
    // Firmer turf by 0.08 in on the TruFirm ran 9–12 ft further at their
    // uphill sites. Compliance 0.9 against 1.0 adds 9.0 ft (26.3 → 35.3)
    // to the uphill 7-iron, the same order.
    let green = whitlark_green();
    let soft = bounce_and_roll_ft(green, 0.015, SEVEN_IRON);
    let firm = bounce_and_roll_ft(
        Surface {
            compliance: 0.9,
            ..green
        },
        0.015,
        SEVEN_IRON,
    );
    assert!((5.0..20.0).contains(&(firm - soft)), "{soft} → {firm}");

    // The impact friction is not a firmness knob for this landing: the
    // Coulomb bound does not bind, so FIRM_GREEN's μ 0.25 leaves the
    // 7-iron where it was (it matters for the high-spin wedge).
    let low_mu = bounce_and_roll_ft(Surface::FIRM_GREEN.with_stimp(12.0), 0.015, SEVEN_IRON);
    assert!((low_mu - soft).abs() < 1.0, "{soft} vs {low_mu}");
}

#[test]
fn cross_break_angle_independent_of_speed_under_constant_law() {
    // Launched 10° off straight uphill, a ball under constant rolling
    // resistance stops at the same angle from its launch line at every
    // speed, and at a distance proportional to v₀²: the motion scales
    // exactly with v₀. On GREEN's stimp: 1.125° and 0.8318 m/(m/s)² on a
    // 1° plane, 2.035° and 0.6941 on a 2° plane.
    //
    // The presets' quadratic law breaks the invariance: over 0.5–3 m/s the
    // angle drifts from 1.36° to 1.12° on the 1° plane and from 2.41° to
    // 2.04° on the 2° plane.
    let surface = same_stimp(Surface::GREEN, RollLaw::Constant);
    let off = units::deg_to_rad(10.0);
    let heading = DVec3::new(off.cos(), 0.0, off.sin());
    for slope_deg in [1.0, 2.0] {
        let grade = units::deg_to_rad(slope_deg).tan();
        let ground = incline(grade, surface);
        let mut first: Option<(f64, f64)> = None;
        for v0 in [0.5, 1.0, 2.0, 3.0] {
            let start = rolling_start(grade, heading, v0);
            let g = ground.settle(&start, &BALL, VACUUM);
            assert_eq!(g.termination, Termination::Rest);
            let d = g.rest.pos - start.pos;
            let angle = units::rad_to_deg(d.angle_between(start.vel));
            let scaled_distance = d.length() / (v0 * v0);
            let (a0, s0) = *first.get_or_insert((angle, scaled_distance));
            assert!(
                (angle - a0).abs() < 0.01,
                "{slope_deg}° {v0} m/s: {angle}° vs {a0}°"
            );
            assert!(
                (scaled_distance / s0 - 1.0).abs() < 1e-3,
                "{slope_deg}° {v0} m/s: {scaled_distance} vs {s0}"
            );
        }
    }
}
