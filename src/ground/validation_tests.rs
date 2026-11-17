//! The ground model against published measurements.
//!
//! - Biber, S. W., Jones, K. M., Champneys, A. R., Green, R. and Szalai, R.
//!   (2023), "Measurements and linearized models for golf ball bounce on a
//!   green", *Sports Engineering* 26:50, doi:10.1007/s12283-023-00442-4,
//!   CC BY 4.0. Bounces on an artificial mat on a rigid base (Campaign A)
//!   and on a natural tee (Campaign B); dimple's contact-angle law is fitted
//!   to Campaign B's bounces. Used here: the "typical driver landing" of
//!   Fig. 3 and the piecewise-affine landing-to-lift-off fits of eqs (8) (A)
//!   and (9) (B), with their error measures, eqs (5)–(6). The paper publishes
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
//!   (0.08 in on the USGA firmness gauge) 9–12 ft further uphill, green speeds 11 ft 3 in
//!   and 12 ft 3 in no different, some wedge distances negative.
//! - USGA and R&A Rules Limited (2021), "Proposed Bounce Model for Use in
//!   Evaluating Optimum Overall Distance": the 2019 range test's roll-out
//!   line, 88.3 − 1.79 α yd for landings at 10–36°, 75–90 ft/s and
//!   1,800 rpm.
//! - Johnson, A. D. G., Otto, S. R., Richardson, P. C. W. and Isaac, S. P.
//!   (2008), "An Investigation into Ball-Turf Impact Characteristics", The
//!   R&A: Table 2's driver landing (64 mph, 36°, 2,250 rpm) and its
//!   bounce and roll on five fairways, 4–31 yd.
//! - USGA and R&A (2007), "Second Report on the Study of Spin Generation",
//!   appendix on golf ball impact with a green: 57 landings on one firm
//!   green, Stimpmeter 10–11 ft, and its printed fits for outbound backspin
//!   and first bounce.
//! - Cross, R. (2024), "Trajectory of a golf ball on a sloping green", *Eur.
//!   J. Phys.* 45, 035002, doi:10.1088/1361-6404/ad2c2d: under constant rolling
//!   friction a ball rolled across an incline stops at a displacement angle
//!   independent of its launch speed, at a distance proportional to v₀².

use super::*;
use crate::aero::Tour;
use crate::air::Wind;
use crate::units;

const BALL: Tour = Tour::SEA_LEVEL;
const R: f64 = 0.04267 / 2.0;

static NO_AIR: Tour = Tour::new(0.0, Wind::None);
static STILL_AIR: Tour = Tour::SEA_LEVEL;

/// No air: ballistic hops, gravity only.
const VACUUM: &Tour = &NO_AIR;

/// Still air at sea level.
fn air() -> &'static Tour {
    &STILL_AIR
}

/// A plane rising at `grade` (rise over run) toward +X, made of `surface`.
fn incline(
    grade: f64,
    surface: Surface,
) -> GroundPhysics<FnTerrain<impl Fn(f64, f64) -> GroundSample + Copy>> {
    let normal = Vec3::new(-grade, 1.0, 0.0).normalize();
    GroundPhysics::new(FnTerrain(move |x: f64, _z: f64| GroundSample {
        height: grade * x,
        normal,
        surface,
    }))
}

/// Angle (rad) between `a` and `b`.
fn angle_between(a: Vec3, b: Vec3) -> f64 {
    (a.dot(b) / (a.length_squared() * b.length_squared()).sqrt())
        .clamp(-1.0, 1.0)
        .acos()
}

/// A ball at the origin of a plane rising at `grade`, rolling without slip
/// at `speed` in the direction `heading` (horizontal, unit), which the
/// velocity follows up or down the plane.
fn rolling_start(grade: f64, heading: Vec3, speed: f64) -> BallState {
    let normal = Vec3::new(-grade, 1.0, 0.0).normalize();
    let along = Vec3::new(heading.x, grade * heading.x, heading.z).normalize();
    let vel = along * speed;
    BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel,
        spin: normal.cross(vel) / R,
    }
}

/// A stimpmeter roll along +X on a plane rising at `grade` (negative:
/// downhill).
fn stimp_roll(surface: Surface, grade: f64) -> GroundResult {
    incline(grade, surface).settle(&rolling_start(grade, Vec3::X, STIMPMETER_SPEED), VACUUM)
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

#[test]
fn brede_harmonic_mean_exact_under_constant_law() {
    // Constant rolling resistance decelerates at g(μ cos θ ± sin θ)/(1+k) up
    // and down the slope, so 2ab/(a+b) = (1+k)v₀²/(2gμ cos θ): the level
    // reading over cos θ, the tape measuring along the surface. Exact only
    // without plowing.
    for surface in [
        TEST_TURF,
        Surface {
            mu_roll: 0.12,
            ..TEST_TURF
        },
    ] {
        let level = surface.stimp_feet();
        for pct in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            let grade: f64 = pct / 100.0;
            if grade >= surface.mu_roll {
                break;
            }
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
fn brede_harmonic_mean_under_plowing() {
    // With plowing the up- and downhill rolls cover different speeds, so
    // the correction is approximate. A green (μ_roll 0.055, plow 0.2 /s,
    // level 9.95 ft) at 1–4%: up/down 8.82/11.42, 7.92/13.42, 7.19/16.32,
    // 6.58/20.93 ft, i.e. +0.03/+0.13/+0.32/+0.68%; a fairway (0.12,
    // 0.36 /s) at most +0.32% to 6%.
    let green = Surface {
        plow: 0.2,
        ..TEST_TURF
    };
    let fairway = Surface {
        mu_roll: 0.12,
        plow: 0.36,
        ..TEST_TURF
    };
    for (surface, max_pct, tolerance) in [
        (green, 3.0, 0.004),
        (green, 4.0, 0.008),
        (fairway, 6.0, 0.004),
    ] {
        let level = surface.stimp_feet();
        for pct in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] {
            if pct > max_pct {
                break;
            }
            let grade = pct / 100.0;
            let (up, down) = (
                stimp_roll_ft(surface, grade),
                stimp_roll_ft(surface, -grade),
            );
            let error = brede(up, down) / level - 1.0;
            assert!(
                error.abs() < tolerance,
                "{pct}%: 2ab/(a+b) {} vs level {level}",
                brede(up, down)
            );
        }
    }
    // The plain average of the two readings, which ignores the slope, is
    // already about 7% long at 2% on that green; Brede's correction removes
    // it.
    let (up, down) = (stimp_roll_ft(green, 0.02), stimp_roll_ft(green, -0.02));
    assert!(0.5 * (up + down) > 1.05 * green.stimp_feet());
}

#[test]
fn green_cannot_hold_downhill_stimp_roll_past_hold_grade() {
    // A green with μ_roll 0.055 (stimp 14.3 ft under the constant law)
    // holds a ball at rest only to a 5.5% grade (tan θ = mu_roll), where
    // Brede held stimpmeter rolls to 5.6% on his carpet runway and read
    // field greens to 6%.
    let sample = |grade: f64| GroundSample {
        height: 0.0,
        normal: Vec3::new(grade, 1.0, 0.0).normalize(),
        surface: TEST_TURF,
    };
    assert!(sample(0.054).holds_at_rest());
    assert!(!sample(0.056).holds_at_rest());

    // At 5% the downhill roll stops, about 157 ft out (11 × level).
    let down5 = stimp_roll_ft(TEST_TURF, -0.05);
    assert!((140.0..175.0).contains(&down5), "{down5} ft");

    // At 6% it never does: under the constant law it gathers speed and is
    // still rolling when the time limit ends it.
    let six = GroundPhysics {
        max_time: 20.0,
        ..incline(-0.06, TEST_TURF)
    }
    .settle(&rolling_start(-0.06, Vec3::X, STIMPMETER_SPEED), VACUUM);
    assert_eq!(six.termination, Termination::MaxTime);
    assert!(
        six.rest.vel.length() > STIMPMETER_SPEED,
        "{}",
        six.rest.vel.length()
    );

    // With plowing 0.2 /s it settles toward the terminal speed where
    // μ_roll·g·cos θ + plow·v = g·sin θ, about 0.24 m/s.
    let plowed = Surface {
        plow: 0.2,
        ..TEST_TURF
    };
    let six = GroundPhysics {
        max_time: 60.0,
        ..incline(-0.06, plowed)
    }
    .settle(&rolling_start(-0.06, Vec3::X, STIMPMETER_SPEED), VACUUM);
    assert_eq!(six.termination, Termination::MaxTime);
    let theta = 0.06_f64.atan();
    let terminal = G * (theta.sin() - plowed.mu_roll * theta.cos()) / plowed.plow;
    assert!(
        (six.rest.vel.length() - terminal).abs() < 0.01,
        "{} vs {terminal}",
        six.rest.vel.length()
    );
}

/// One branch of a piecewise-affine bounce fit, `p_F` = A·p₀ + b, in Biber et
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
        pos: Vec3::ZERO,
        vel: Vec3::new(speed * a.cos(), -speed * a.sin(), 0.0),
        spin: Vec3::new(0.0, 0.0, units::rpm_to_rads(rpm)),
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
/// δ = mean of |`p_F` − `P_F`| / |p₀|, and per component Δq = mean of
/// |`q_F` − `Q_F`| / |p₀|. Returns (δ, Δẋ, Δẏ, Δω).
#[allow(clippy::cast_precision_loss, reason = "a few dozen bounces")]
fn errors(b: &[Bounce]) -> (f64, f64, f64, f64) {
    let n = b.len() as f64;
    let mut e = (0.0, 0.0, 0.0, 0.0);
    for x in b {
        let d = [
            x.fit[0] - x.model[0],
            x.fit[1] - x.model[1],
            x.fit[2] - x.model[2],
        ];
        let p0 = dot(x.p0, x.p0).sqrt();
        e.0 += dot(d, d).sqrt() / p0;
        e.1 += d[0].abs() / p0;
        e.2 += d[1].abs() / p0;
        e.3 += d[2].abs() / p0;
    }
    (e.0 / n, e.1 / n, e.2 / n, e.3 / n)
}

/// Mean normal rebound ratio |`ẏ_F` / ẏ₀| of the fit and of dimple.
#[allow(clippy::cast_precision_loss, reason = "a few dozen bounces")]
fn rebound(b: &[Bounce]) -> (f64, f64) {
    let n = b.len() as f64;
    let fit = b.iter().map(|x| x.fit[1] / -x.p0[1]).sum::<f64>() / n;
    let model = b.iter().map(|x| x.model[1] / -x.p0[1]).sum::<f64>() / n;
    (fit, model)
}

/// The tee-turf constants of the contact law: compliance 1.0 and sliding
/// friction 0.35, fitted with it to Campaign B's 693 bounces.
const TEE_TURF: Surface = Surface {
    mu_slide: 0.35,
    ..TEST_TURF
};

#[test]
fn biber_tee_turf_fit_agreement() {
    // Natural turf (Campaign B). The tee constants against the fit: δ 10.5%
    // (Δẋ 5.5, Δẏ 2.7, Δω 7.2%); the paper's own physics fits reach 19–35%
    // against the data (Table 3), its affine fit 7.8% (Table 2). The
    // presets, whose friction is set on greens and fairways, score 17–22%
    // on this turf.
    let landings = tee_grid();
    let tee = bounces(&TEE, &landings, TEE_TURF);
    let (delta, _, dy, _) = errors(&tee);
    assert!(delta < 0.12, "δ {delta}");
    assert!(dy < 0.05, "Δẏ {dy}");
    for preset in [Surface::GREEN, Surface::FAIRWAY, Surface::FIRM_GREEN] {
        let (preset_delta, ..) = errors(&bounces(&TEE, &landings, preset));
        assert!(preset_delta < 0.25, "{preset:?}: δ {preset_delta}");
    }

    // Normal rebound: 0.371 against the fit's 0.408.
    let (fit_rebound, model_rebound) = rebound(&tee);
    assert!(
        (fit_rebound - model_rebound).abs() < 0.05,
        "{model_rebound} vs {fit_rebound}"
    );

    // Spin reversal: of the 80 backspin landings, 19 leave with topspin in
    // dimple and 23 in the fit.
    let backspin: Vec<_> = tee.iter().filter(|b| b.p0[2] > 0.0).collect();
    let reversed = |q: fn(&Bounce) -> f64| backspin.iter().filter(|b| q(b) < 0.0).count();
    let (fit_rev, model_rev) = (reversed(|b| b.fit[2]), reversed(|b| b.model[2]));
    assert!(fit_rev.abs_diff(model_rev) <= 6, "{model_rev} vs {fit_rev}");

    // Balls that enter rolling lift off slipping with excess topspin,
    // contact-point velocity v_P = ẋ + ω < 0 (their Fig. 6). dimple's
    // tilted contact plane gives the same sign in the ground frame for 55
    // of the 57 roll-branch landings.
    let roll: Vec<_> = tee.iter().filter(|b| b.rolls).collect();
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

    // A rigid plane (no tilt, 21.9%) and the fits' μ ≈ 1 (22.0%) both do
    // worse than the tilted contact plane with μ 0.35.
    for worse in [
        Surface {
            compliance: 0.0,
            ..TEE_TURF
        },
        Surface {
            mu_slide: 1.0,
            ..TEE_TURF
        },
    ] {
        let (d, ..) = errors(&bounces(&TEE, &landings, worse));
        assert!(d > delta + 0.05, "{worse:?}: δ {d} vs {delta}");
    }
}

#[test]
fn biber_artificial_turf_needs_higher_restitution() {
    // A mat on a rigid base (Campaign A) rebounds far more than Penner's
    // turf: the fit's mean normal rebound ratio is 0.58, GREEN's 0.34
    // (Table 3 fits r = 0.42–0.54 for the mat). So no turf preset is a mat.
    let landings = mat_grid();
    let green = bounces(&MAT, &landings, Surface::GREEN);
    let (fit_rebound, model_rebound) = rebound(&green);
    assert!(fit_rebound > 0.55, "{fit_rebound}");
    assert!(model_rebound < 0.4, "{model_rebound}");

    // A constant restitution with a lower contact angle and friction fits
    // it better (δ 29.4% → 16.8%, Δẏ 17.1% → 3.7%).
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
    // of backspin. On GREEN the normal rebound is 0.356 against the tee
    // fit's 0.342, the forward ratio 0.236 against 0.277, and the ball
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
    let g = incline(grade, surface).settle(&landing(units::mph_to_mps(mph), descent, rpm), air());
    assert_eq!(g.termination, Termination::Rest);
    units::meters_to_feet(g.rest.pos.x * grade.hypot(1.0))
}

#[test]
fn whitlark_slope_asymmetry_on_a_green() {
    // 7-iron on a 1.5% grade: uphill 23.7 ft, level 31.0 ft, downhill
    // 41.1 ft. Whitlark & Pringle measured mean totals of 13.5 ft uphill
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

    // The high-spin sand wedge stops near its pitch mark on every slope:
    // 1.3 ft behind it uphill, 3.3 ft past it level, 6.4 ft past it
    // downhill; "in some cases this distance was negative".
    let wedges = [0.015, 0.0, -0.015].map(|grade| bounce_and_roll_ft(green, grade, SAND_WEDGE));
    for wedge in wedges {
        assert!(wedge.abs() < 8.0, "{wedge} ft");
    }
    assert!(wedges[0] < 0.0, "uphill {} ft", wedges[0]);
}

#[test]
fn whitlark_green_speed_barely_moves_bounce_and_roll() {
    // Stimp 11 ft 3 in against 12 ft 3 in made no difference to their
    // bounce and roll; dimple's uphill 7-iron moves 23.1 → 23.9 ft.
    let slow = bounce_and_roll_ft(Surface::GREEN.with_stimp(11.25), 0.015, SEVEN_IRON);
    let fast = bounce_and_roll_ft(Surface::GREEN.with_stimp(12.25), 0.015, SEVEN_IRON);
    assert!((fast - slow).abs() < 2.0, "{slow} vs {fast}");
}

#[test]
fn whitlark_compliance_is_the_firmness_knob() {
    // Firmer turf by 0.08 in (2 mm) on the TruFirm ran 9–12 ft further at
    // their uphill sites. On greens compliance falls 0.090 per mm of
    // penetration (R&A 2008), so 2 mm firmer is compliance 0.18 lower: it
    // adds 17.2 ft (23.7 → 40.9) to the uphill 7-iron, the same order.
    let green = whitlark_green();
    let soft = bounce_and_roll_ft(green, 0.015, SEVEN_IRON);
    let firm = bounce_and_roll_ft(
        Surface {
            compliance: green.compliance - 2.0 * 0.090,
            ..green
        },
        0.015,
        SEVEN_IRON,
    );
    assert!((5.0..25.0).contains(&(firm - soft)), "{soft} → {firm}");

    // The impact friction is hardly a firmness knob for this landing:
    // FIRM_GREEN's μ 0.17 against GREEN's 0.22 moves the uphill 7-iron
    // 2.3 ft (it matters more for the high-spin wedge: 5.9 ft downhill).
    let low_mu = bounce_and_roll_ft(
        Surface {
            mu_slide: Surface::FIRM_GREEN.mu_slide,
            ..green
        },
        0.015,
        SEVEN_IRON,
    );
    assert!((low_mu - soft).abs() < 3.0, "{soft} vs {low_mu}");
}

#[test]
fn cross_break_angle_independent_of_speed_under_constant_law() {
    // Launched 10° off straight uphill, a ball under constant rolling
    // resistance stops at the same angle from its launch line at every
    // speed, and at a distance proportional to v₀²: the motion scales
    // exactly with v₀. With μ_roll 0.055: 1.368° and 0.9897 m/(m/s)² on a
    // 1° plane, 2.428° and 0.8008 on a 2° plane.
    //
    // Plowing breaks the invariance: on `GREEN` over 0.5–3 m/s the angle
    // drifts from 1.37° to 1.06° on the 1° plane and from 2.43° to 1.93° on
    // the 2° plane.
    let surface = TEST_TURF;
    let off = units::deg_to_rad(10.0);
    let heading = Vec3::new(off.cos(), 0.0, off.sin());
    for slope_deg in [1.0, 2.0] {
        let grade = units::deg_to_rad(slope_deg).tan();
        let ground = incline(grade, surface);
        let mut first: Option<(f64, f64)> = None;
        for v0 in [0.5, 1.0, 2.0, 3.0] {
            let start = rolling_start(grade, heading, v0);
            let g = ground.settle(&start, VACUUM);
            assert_eq!(g.termination, Termination::Rest);
            let d = g.rest.pos - start.pos;
            let angle = units::rad_to_deg(angle_between(d, start.vel));
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

#[test]
fn usga_2019_runner_line_on_fairway() {
    // The USGA's 2019 range test (USGA/R&A 2021): roll-out from landings at
    // 75–90 ft/s and 1,800 rpm follows 88.3 − 1.79 α yd. FAIRWAY, averaging
    // landings at 75 and 90 ft/s, is within 3 yd of the line at 20–36°
    // (rms 2 yd). A guard on a fit target: FAIRWAY's `plow` is fitted to
    // this line, and it bounds the contact angle's θ0 and power.
    let roll_out = |alpha: f64| {
        let total: f64 = [75.0, 90.0]
            .iter()
            .map(|&fps| {
                let l = landing(units::feet_to_meters(fps), alpha, 1800.0);
                let g = Surface::FAIRWAY.settle(&l, air());
                assert_eq!(g.termination, Termination::Rest);
                units::meters_to_yards(g.rest.pos.x)
            })
            .sum();
        total / 2.0 - (88.3 - 1.79 * alpha)
    };
    let mut sum_sq = 0.0;
    let mut n = 0.0;
    for alpha in (20..=36).step_by(2).map(f64::from) {
        let error = roll_out(alpha);
        assert!(error.abs() < 3.0, "{alpha}°: {error:+.1} yd");
        sum_sq += error * error;
        n += 1.0;
    }
    let rms = (sum_sq / n).sqrt();
    assert!(rms < 2.5, "rms {rms:.2} yd");
}

#[test]
fn rnda_2008_driver_bounce_and_roll_on_fairways() {
    // The R&A's driver landing (2008, Table 2): 64 mph, 36°, 2,250 rpm onto
    // five fairways bounced and rolled 4–31 yd by firmness (median 24 yd).
    // FAIRWAY: 24.8 yd. A guard on a fit target: FAIRWAY's `compliance`,
    // `mu_slide` and `plow` are fitted to these fairways.
    let g = Surface::FAIRWAY.settle(&landing(units::mph_to_mps(64.0), 36.0, 2250.0), air());
    let yards = units::meters_to_yards(g.rest.pos.x);
    assert!((4.0..31.0).contains(&yards), "{yards:.1} yd");
}

#[test]
fn usga_2006_firm_green_printed_fits() {
    // 57 landings on a firm green at 66, 75 and 85 ft/s, 36–59° and
    // 0–10,000 rpm (USGA and R&A 2007). The report's fits: outbound backspin
    // = inbound − 4,690 rpm, first bounce = 72.3 − 1.2 θ ft. FIRM_GREEN over
    // that grid: backspin out within 1,250 rpm of the fit (mean +384 rpm),
    // first bounce within 7 ft (mean −1.8 ft). The measured first bounces
    // span 3–32 ft. A guard on a fit target: FIRM_GREEN's `compliance` and
    // `mu_slide` and the greens' `plow` are fitted to this green.
    let mut spin_sum = 0.0;
    let mut bounce_sum = 0.0;
    let mut n = 0.0;
    for fps in [66.0, 75.0, 85.0] {
        for theta in [35.0, 40.0, 45.0, 50.0, 55.0] {
            for rpm in [2000.0, 5000.0, 8000.0, 10000.0] {
                let l = landing(units::feet_to_meters(fps), theta, rpm);
                let g = Surface::FIRM_GREEN.settle(&l, air());
                let first = g.contacts[0];
                let out_rpm = first.spin_out.z * 30.0 / std::f64::consts::PI;
                let spin_error = out_rpm - (rpm - 4690.0);
                let second = g.contacts.get(1).map_or(g.rest.pos.x, |c| c.pos.x);
                let bounce = units::meters_to_feet(second - first.pos.x);
                let bounce_error = bounce - (72.3 - 1.2 * theta);
                let case = format!("{fps} ft/s, {theta}°, {rpm} rpm");
                assert!(
                    spin_error.abs() < 1250.0,
                    "{case}: spin {spin_error:+.0} rpm"
                );
                assert!(
                    bounce_error.abs() < 7.0,
                    "{case}: bounce {bounce_error:+.1} ft"
                );
                spin_sum += spin_error;
                bounce_sum += bounce_error;
                n += 1.0;
            }
        }
    }
    assert!(
        (spin_sum / n).abs() < 500.0,
        "mean spin {:+.0}",
        spin_sum / n
    );
    assert!(
        (bounce_sum / n).abs() < 3.0,
        "mean bounce {:+.1}",
        bounce_sum / n
    );
}
