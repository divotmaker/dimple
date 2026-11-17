//! Generates `docs/FIT.md`: how dimple's ground model scores against each
//! dataset its laws are fitted to or checked against.
//!
//! ```text
//! cargo run --example fit_doc   # fixture rows + docs/fit-summary.csv → docs/FIT.md
//! ```
//!
//! **Fixture** rows use only numbers the sources print (fitted equations,
//! table rows, lines, Stimpmeter readings), held as constants below, so they
//! regenerate anywhere. **Summary** rows score per-impact datasets that
//! dimple does not distribute; the maintainers score them, and only their
//! size, ranges and the resulting metrics are published, in
//! `docs/fit-summary.csv`, from which they are read here. The test at the
//! bottom recomputes the fixture rows and compares them with that file.
//!
//! Everything here goes through the public API: [`impact`],
//! [`GroundModel::settle`] and the [`Surface`] presets.

#![allow(
    clippy::many_single_char_names,
    reason = "the impact formulas keep the papers' symbols"
)]

use std::f64::consts::PI;
use std::fmt::Write as _;
use std::path::PathBuf;

use dimple::ground::{GroundModel, GroundSample, impact};
use dimple::{BallState, Surface, Tour, Vec3};

const FPS: f64 = 0.3048;
const YD: f64 = 0.9144;
const RPM: f64 = 2.0 * PI / 60.0;
/// Ball radius used for the scaled state of Biber et al.'s error measure
/// (the model's own radius).
const R: f64 = 0.042_67 / 2.0;
/// Cross's ball radius and inertia factor, for his ey, COF and D.
const RC: f64 = 0.0214;
const AC: f64 = 0.4;
const AIR: Tour = Tour::SEA_LEVEL;

// ---------------------------------------------------------------- model

/// Outbound (horizontal m/s, vertical m/s up, backspin rad/s) of one impact
/// on level ground: landing `speed` m/s at `deg` below horizontal with
/// `back` rad/s of backspin.
fn bounce(s: Surface, speed: f64, deg: f64, back: f64) -> [f64; 3] {
    let a = deg.to_radians();
    let state = landing(speed, a, back);
    let at = GroundSample {
        height: 0.0,
        normal: Vec3::new(0.0, 1.0, 0.0),
        surface: s,
    };
    let out = impact(&state, &at, &AIR);
    [out.vel.x, out.vel.y, out.spin.z]
}

fn landing(speed: f64, a: f64, back: f64) -> BallState {
    BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel: Vec3::new(speed * a.cos(), -speed * a.sin(), 0.0),
        spin: Vec3::new(0.0, 0.0, back),
    }
}

/// Settles a landing: (first bounce length, total bounce and roll), meters.
fn settle(s: Surface, speed: f64, deg: f64, rpm: f64) -> (f64, f64) {
    let g = s.settle(&landing(speed, deg.to_radians(), rpm * RPM), &AIR);
    let first = g.contacts.get(1).map_or(g.rest.pos.x, |c| c.pos.x);
    (first, g.rest.pos.x)
}

/// Cross's metrics from an inbound and outbound vector: ey, COF, D (mm).
fn cross_metrics(v1: f64, th1: f64, back1: f64, out: [f64; 3]) -> [f64; 3] {
    let a = th1.to_radians();
    let (vx1, vy1) = (v1 * a.cos(), v1 * a.sin());
    let (w1, w2) = (-back1, -out[2]);
    let ey = out[1] / vy1;
    let cof = (vx1 - out[0]) / (vy1 + out[1]);
    let d = (RC * (vx1 - out[0]) - AC * RC * RC * (w2 - w1)) / (vy1 + out[1]) * 1000.0;
    [ey, cof, d]
}

// ---------------------------------------------------------------- scoring

/// One measured impact: inbound speed, angle, backspin (rad/s); outbound
/// horizontal and vertical velocity and backspin.
#[derive(Clone, Copy)]
struct Shot {
    v: f64,
    deg: f64,
    back: f64,
    out: [f64; 3],
}

impl Shot {
    fn from_polar(v: f64, deg: f64, back_rpm: f64, vo: f64, ango: f64, out_back_rpm: f64) -> Self {
        let ao = ango.to_radians();
        Self {
            v,
            deg,
            back: back_rpm * RPM,
            out: [vo * ao.cos(), vo * ao.sin(), out_back_rpm * RPM],
        }
    }
    fn inbound(&self) -> (f64, f64) {
        let a = self.deg.to_radians();
        (self.v * a.cos(), -self.v * a.sin())
    }
    /// Error of a predicted outbound vector, Biber et al.'s δ: % of the
    /// inbound norm in scaled units (radii/s, rad/s).
    fn delta(&self, p: [f64; 3]) -> f64 {
        let (vx, vy) = self.inbound();
        let n0 = ((vx / R).powi(2) + (vy / R).powi(2) + self.back.powi(2)).sqrt();
        let e = [
            (p[0] - self.out[0]) / R,
            (p[1] - self.out[1]) / R,
            p[2] - self.out[2],
        ];
        100.0 * (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]).sqrt() / n0
    }
}

#[allow(clippy::cast_precision_loss, reason = "dataset sizes")]
fn len<T>(v: &[T]) -> f64 {
    v.len() as f64
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / len(v)
}

fn rms(v: &[f64]) -> f64 {
    (v.iter().map(|x| x * x).sum::<f64>() / len(v)).sqrt()
}

/// Least-squares `y = c0 + c1·(v − 23) + c2·(θ − 45) + c3·(ω + 270)`, Cross's
/// form (ω in rad/s, topspin positive). Returns c0, c1, c2.
#[allow(
    clippy::needless_range_loop,
    reason = "Gauss-Jordan elimination reads clearest with indices"
)]
fn linfit(y: &[f64], pts: &[(f64, f64, f64)]) -> [f64; 3] {
    let mut m = [[0.0; 5]; 4];
    for (&yi, &(v, t, w)) in y.iter().zip(pts) {
        let x = [1.0, v - 23.0, t - 45.0, w + 270.0];
        for i in 0..4 {
            for j in 0..4 {
                m[i][j] += x[i] * x[j];
            }
            m[i][4] += x[i] * yi;
        }
    }
    for i in 0..4 {
        let p = m[i][i];
        for k in i..5 {
            m[i][k] /= p;
        }
        for r in 0..4 {
            if r != i {
                let f = m[r][i];
                for k in i..5 {
                    m[r][k] -= f * m[i][k];
                }
            }
        }
    }
    [m[0][4], m[1][4], m[2][4]]
}

// ---------------------------------------------------------------- rows

/// One row of the fit table.
#[derive(Clone, Debug, PartialEq)]
struct Row {
    id: String,
    status: String,
    law: String,
    dataset: String,
    range: String,
    metrics: String,
}

fn row(id: &str, status: &str, law: &str, dataset: &str, range: &str, metrics: String) -> Row {
    Row {
        id: id.into(),
        status: status.into(),
        law: law.into(),
        dataset: dataset.into(),
        range: range.into(),
        metrics,
    }
}

/// Signed with a typographic minus.
fn sg(x: f64, dp: usize) -> String {
    let s = format!("{x:+.dp$}");
    s.replacen('-', "−", 1)
}

fn un(x: f64, dp: usize) -> String {
    format!("{x:.dp$}").replacen('-', "−", 1)
}

// Cross 2025, eqs. 14–16 (v m/s, θ deg, ω rad/s with topspin positive).
fn eq14(v: f64, t: f64, w: f64) -> f64 {
    0.301 - 0.0047 * (v - 23.0) - 0.0083 * (t - 44.9) + 6.553e-5 * (w + 275.0)
}
fn eq15(v: f64, t: f64, w: f64) -> f64 {
    0.574 + 0.0095 * (v - 23.0) - 0.009 * (t - 45.0) - 0.0002 * (w + 270.0)
}
fn eq16(v: f64, t: f64, w: f64) -> f64 {
    8.23 + 0.30 * (v - 23.1) - 0.14 * (t - 45.5) - 0.0009 * (w + 255.0)
}

/// Cross 2025, Table 1: inbound speed, angle, spin (rpm, backspin negative);
/// outbound the same; printed ey, COF, D (mm).
const TABLE1: [[f64; 9]; 4] = [
    [10.12, 31.1, -3344.0, 4.67, 31.0, -358.0, 0.460, 0.611, 5.46],
    [9.70, 36.2, -1847.0, 4.47, 33.1, 967.0, 0.426, 0.501, 4.11],
    [9.10, 45.1, -789.0, 3.88, 40.3, 1432.0, 0.386, 0.387, 3.52],
    [3.42, 45.2, 1379.0, 2.71, 24.4, 1061.0, 0.459, -0.017, 1.37],
];

const CROSS25: &str = "Cross (2025)";
const GREEN_RANGE: &str =
    "17–27 m/s, 35–55°, 0–7000 rpm backspin (the range of Cross's fits); 135-point grid";

/// The rows built from printed numbers only.
#[allow(clippy::too_many_lines, reason = "one block per fixture row")]
fn fixture_rows() -> Vec<Row> {
    let mut rows = Vec::new();

    // Cross 2025 eqs. 14–16 on the grid of his fits.
    let mut grid = Vec::new();
    for v in [17.0, 22.0, 27.0] {
        for t in [35.0, 45.0, 55.0] {
            for i in 0..15 {
                grid.push((v, t, f64::from(i) * 500.0 * RPM));
            }
        }
    }
    let pts: Vec<_> = grid.iter().map(|&(v, t, b)| (v, t, -b)).collect();
    let metric = |s: Surface, k: usize| -> Vec<f64> {
        grid.iter()
            .map(|&(v, t, b)| cross_metrics(v, t, b, bounce(s, v, t, b))[k])
            .collect()
    };
    let printed: [fn(f64, f64, f64) -> f64; 3] =
        [|v, t, w| eq14(v, t, w).clamp(0.1, 0.5), eq15, eq16];
    let minus = |s: Surface, k: usize| -> (f64, f64) {
        let d: Vec<f64> = metric(s, k)
            .iter()
            .zip(&pts)
            .map(|(m, &(v, t, w))| m - printed[k](v, t, w))
            .collect();
        (mean(&d), rms(&d))
    };
    let g = Surface::GREEN;
    for (id, law, k, name, cross, dp, unit) in [
        (
            "cross2025-eq16",
            "Contact angle",
            2,
            "D (normal-force offset)",
            [8.23, 0.30, -0.14],
            2,
            " mm",
        ),
        (
            "cross2025-eq14",
            "Restitution",
            0,
            "ey",
            [0.301, -0.0047, -0.0083],
            3,
            "",
        ),
        (
            "cross2025-eq15",
            "Friction (`mu_slide`)",
            1,
            "COF",
            [0.574, 0.0095, -0.0090],
            3,
            "",
        ),
    ] {
        let c = linfit(&metric(g, k), &pts);
        let (b, r) = minus(g, k);
        let sdp = if k == 2 { 2 } else { 4 };
        let eq = match k {
            0 => "14",
            1 => "15",
            _ => "16",
        };
        rows.push(row(
            id,
            "fixture, check",
            law,
            &format!("{CROSS25}, eq. {eq} (`GREEN`)"),
            GREEN_RANGE,
            format!(
                "{name}: mean {}{unit} (eq. {eq}: {}); per m/s {} ({}); per degree {} ({}); model minus eq. {eq}: {} / {}{unit} (mean / rms)",
                un(c[0], dp),
                un(cross[0], dp),
                sg(c[1], sdp),
                sg(cross[1], sdp),
                sg(c[2], sdp),
                sg(cross[2], sdp),
                sg(b, dp),
                un(r, dp),
            ),
        ));
    }

    // Sliding COF at 35°.
    let cof = |s: Surface| -> String {
        [17.0, 22.0, 27.0]
            .iter()
            .map(|&v| {
                let b = 7000.0 * RPM;
                un(cross_metrics(v, 35.0, b, bounce(s, v, 35.0, b))[1], 2)
            })
            .collect::<Vec<_>>()
            .join(" / ")
    };
    rows.push(row(
        "cross2025-slide",
        "fixture, check",
        "Friction (`mu_slide`)",
        &format!("{CROSS25}, μ from sliding bounces"),
        "35°, 7000 rpm backspin, 17 / 22 / 27 m/s",
        format!(
            "COF on the level surface: `GREEN` {}, `FIRM_GREEN` {} (Cross: 0.7 ± 0.1; 0.6 on a firm green, 0.8 on a soft one)",
            cof(Surface::GREEN),
            cof(Surface::FIRM_GREEN)
        ),
    ));

    // Cross 2025 Table 1.
    let shots: Vec<Shot> = TABLE1
        .iter()
        .map(|r| Shot::from_polar(r[0], r[1], -r[2], r[3], r[4], -r[5]))
        .collect();
    let (mut parts, mut deltas) = (Vec::new(), Vec::new());
    for (r, x) in TABLE1.iter().zip(&shots) {
        let p = bounce(g, x.v, x.deg, x.back);
        deltas.push(x.delta(p));
        let m = cross_metrics(x.v, x.deg, x.back, p);
        parts.push(format!(
            "{} m/s at {}°: out {} m/s, {}°, {} rpm (measured {}, {}, {}); ey {} ({}), COF {} ({}), δ {}",
            un(r[0], 2),
            un(r[1], 1),
            un(p[0].hypot(p[1]), 2),
            un(p[1].atan2(p[0]).to_degrees(), 1),
            un(-p[2] / RPM, 0),
            un(r[3], 2),
            un(r[4], 1),
            un(r[5], 0),
            un(m[0], 3),
            un(r[6], 3),
            un(m[1], 3),
            un(r[7], 3),
            un(x.delta(p), 1),
        ));
    }
    rows.push(row(
        "cross2025-table1",
        "fixture, check",
        "Friction (`mu_slide`), whole impact",
        &format!("{CROSS25}, Table 1 (`GREEN`; each row the mean of at least five shots)"),
        "3.4–10.1 m/s, 31–45° (spin: rpm, topspin positive)",
        format!("δ {}; {}", un(mean(&deltas), 1), parts.join("; ")),
    ));

    // USGA and R&A 2007: the firm green's printed fits.
    let (mut spin, mut first) = (Vec::new(), Vec::new());
    for fps in [66.0, 75.0, 85.0] {
        for theta in [35.0, 40.0, 45.0, 50.0, 55.0] {
            for rpm in [2000.0, 5000.0, 8000.0, 10000.0] {
                let p = bounce(Surface::FIRM_GREEN, fps * FPS, theta, rpm * RPM);
                spin.push(p[2] / RPM - (rpm - 4690.0));
                let (fb, _) = settle(Surface::FIRM_GREEN, fps * FPS, theta, rpm);
                first.push(fb / FPS - (72.3 - 1.2 * theta));
            }
        }
    }
    rows.push(row(
        "usga2007-fits",
        "fixture, fit",
        "Whole impact, `FIRM_GREEN`",
        "USGA and The R&A (2007), the report's two fits",
        "66–85 ft/s, 35–55°, 2000–10000 rpm; 60-point grid",
        format!(
            "outbound backspin minus (inbound − 4690 rpm): {} / {} rpm (mean / rms); first bounce minus (72.3 − 1.2 θ ft): {} / {} ft",
            sg(mean(&spin), 0),
            un(rms(&spin), 0),
            sg(mean(&first), 1),
            un(rms(&first), 1)
        ),
    ));

    // Stimpmeter readings.
    rows.push(row(
        "stimp",
        "fixture, fit",
        "Roll (`mu_roll`)",
        "USGA and The R&A (2020), fairway reading; USGA and The R&A (2007), green reading",
        "Stimpmeter, 1.83 m/s exit speed",
        format!(
            "`FAIRWAY` reads {} ft (measured 5 ft 8 in = 5.67 ft); `GREEN` and `FIRM_GREEN` read {} ft (measured 10–11 ft)",
            un(Surface::FAIRWAY.stimp_feet(), 2),
            un(Surface::GREEN.stimp_feet(), 2)
        ),
    ));

    // USGA 2019 range test: 88.3 − 1.79 α yd.
    let total = |s: Surface, fps: f64, a: f64, rpm: f64| settle(s, fps * FPS, a, rpm).1 / YD;
    let (mut hi, mut lo) = (Vec::new(), Vec::new());
    for a in (12..=36).step_by(2).map(f64::from) {
        let t = 0.5
            * (total(Surface::FAIRWAY, 75.0, a, 1800.0) + total(Surface::FAIRWAY, 90.0, a, 1800.0));
        let e = t - (88.3 - 1.79 * a);
        if a < 20.0 { lo.push(e) } else { hi.push(e) }
    }
    rows.push(row(
        "usga2019-line",
        "fixture, fit",
        "Run-out, end to end (`FAIRWAY`)",
        "USGA and R&A Rules Limited (2021): 2019 range test, roll-out 88.3 − 1.79 α yd",
        "10–36° (compared at 12–36°, every 2°), 75–90 ft/s, 1800 rpm (model: mean of 75 and 90 ft/s)",
        format!(
            "20–36°: {} / {} yd (bias / rms), largest {} yd; 12–18°: {} yd bias (no impact is measured below 16°)",
            sg(mean(&hi), 1),
            un(rms(&hi), 1),
            un(hi.iter().fold(0.0, |m: f64, e| m.max(e.abs())), 1),
            sg(mean(&lo), 1)
        ),
    ));

    // Philadelphia Cricket Club lines and the PGA TOUR line.
    let angles = [20.0, 30.0, 39.0, 50.0];
    let model: Vec<f64> = angles
        .iter()
        .map(|&a| total(Surface::FAIRWAY, 100.0, a, 1800.0))
        .collect();
    let vs = |i0: f64, s: f64| -> String {
        let e: Vec<f64> = angles
            .iter()
            .zip(&model)
            .map(|(&a, m)| m - (i0 + s * a))
            .collect();
        format!("{} / {}", sg(mean(&e), 1), un(rms(&e), 1))
    };
    rows.push(row(
        "philadelphia-lines",
        "fixture, check",
        "Run-out, end to end (`FAIRWAY`)",
        "USGA and R&A Rules Limited (2021): Philadelphia Cricket Club tests of 2011, firmer hole 110.8 − 1.91 α, softer hole 100.8 − 2.01 α yd",
        "20–50°, 90–125 ft/s (model: 100 ft/s, 1800 rpm)",
        format!(
            "model {} yd at 20 / 30 / 39 / 50°; minus the firmer line {} yd, the softer {} yd (bias / rms)",
            model.iter().map(|x| un(*x, 1)).collect::<Vec<_>>().join(" / "),
            vs(110.8, -1.91),
            vs(100.8, -2.01)
        ),
    ));
    let tour = total(Surface::FAIRWAY, 100.0, 39.0, 2500.0);
    rows.push(row(
        "tour-line",
        "fixture, check",
        "Run-out, end to end (`FAIRWAY`)",
        "USGA and R&A Rules Limited (2021): PGA TOUR 2016–18 mean roll-out 81 − 1.52 α yd",
        "39° (model: 100 ft/s, 2500 rpm)",
        format!(
            "{} yd against the line's {} yd",
            un(tour, 1),
            un(81.0 - 1.52 * 39.0, 1)
        ),
    ));
    rows
}

// ---------------------------------------------------------------- csv and render

const COLUMNS: [&str; 6] = ["id", "status", "law", "dataset", "range", "metrics"];

fn quote(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "\"\""))
}

fn write_csv(rows: &[Row]) -> String {
    let mut out = COLUMNS.join(",") + "\n";
    for r in rows {
        let f = [&r.id, &r.status, &r.law, &r.dataset, &r.range, &r.metrics];
        out += &f.iter().map(|s| quote(s)).collect::<Vec<_>>().join(",");
        out.push('\n');
    }
    out
}

fn parse_csv(text: &str) -> Vec<Row> {
    let mut records = Vec::new();
    let (mut rec, mut field, mut quoted) = (Vec::new(), String::new(), false);
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        match (c, quoted) {
            ('"', true) if chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            ('"', _) => quoted = !quoted,
            (',', false) => rec.push(std::mem::take(&mut field)),
            ('\n', false) => {
                rec.push(std::mem::take(&mut field));
                records.push(std::mem::take(&mut rec));
            }
            _ => field.push(c),
        }
    }
    records
        .into_iter()
        .skip(1)
        .filter(|r| r.len() == COLUMNS.len())
        .map(|r| Row {
            id: r[0].clone(),
            status: r[1].clone(),
            law: r[2].clone(),
            dataset: r[3].clone(),
            range: r[4].clone(),
            metrics: r[5].clone(),
        })
        .collect()
}

/// Display order of the rows (laws first, then whole impacts, roll, run-out).
const ORDER: [&str; 22] = [
    "biber2023-tilt",
    "cross2025-eq16",
    "biber2023-e",
    "cross2025-eq14",
    "cross2021-drops",
    "cross2025-eq15",
    "cross2025-slide",
    "cross2025-table1",
    "randa2008-sw-mu",
    "randa2008-5i",
    "pinehurst2005-greens",
    "usga2007-green",
    "usga2007-fits",
    "pinehurst2005-fairways",
    "randa2008-driver",
    "cross2021-fairway",
    "stimp",
    "usga2007-roll",
    "randa2008-roll",
    "usga2019-line",
    "philadelphia-lines",
    "tour-line",
];

fn sort_rows(rows: &mut [Row]) {
    rows.sort_by_key(|r| {
        ORDER
            .iter()
            .position(|id| *id == r.id)
            .unwrap_or(ORDER.len())
    });
}

fn render(rows: &[Row]) -> String {
    let mut s = String::new();
    s += "# Fit quality of the ground model\n\n";
    s += "<!-- Generated by `cargo run --example fit_doc`; do not edit by hand. -->\n\n";
    s += "How dimple's ground model scores against every dataset its laws are fitted to or checked\n";
    s += "against, one row per law and dataset, with the presets as they ship. The model and the\n";
    s += "constants are described in [AERO.md](AERO.md#ground-model); the sources are listed under\n";
    s += "References in the [README](../README.md#references).\n\n";
    s += "**Regenerating.** `cargo run --example fit_doc` recomputes the **fixture** rows and rewrites\n";
    s += "this file and `docs/fit-summary.csv`. Fixture rows use only numbers the sources print (fitted\n";
    s += "equations, table rows, lines, Stimpmeter readings), so they regenerate anywhere; a test fails\n";
    s += "when a change to the model moves one of them until this file is regenerated. **Summary** rows\n";
    s += "are scored by the maintainers on per-impact datasets that dimple does not include or\n";
    s += "redistribute: only their size, ranges and the resulting metrics are published, read from\n";
    s += "`docs/fit-summary.csv`.\n\n";
    s += "**Status.** Each row is `fixture` or `summary` as above, and `fit` or `check`. `fit` marks a\n";
    s += "dataset the constants were fitted to (in-sample: the row shows how closely the fit\n";
    s += "reproduces it); `check` marks a dataset no constant was fitted to (held out).\n\n";
    s += "**Measures.** δ is the mean error of the outbound vector (both velocity components and the\n";
    s += "spin, in scaled units) as a percentage of the inbound one, the measure of Biber et al. (2023).\n";
    s += "ey is the vertical rebound ratio, COF the horizontal over the vertical impulse and D the\n";
    s += "offset of the normal force ahead of the ball's center, all on the level surface (Cross 2025).\n";
    s += "Biases are model minus measured.\n\n";
    s += "| Law | Dataset | Measured range | Metrics | Status |\n|---|---|---|---|---|\n";
    for r in rows {
        let _ = writeln!(
            s,
            "| {} | {} | {} | {} | {} |",
            r.law, r.dataset, r.range, r.metrics, r.status
        );
    }
    s += "\n## Where the model runs outside the data\n\n";
    s += "- **Landings flatter than 16°.** No impact is measured below 16°. On the USGA's 2019 range\n";
    s += "  test the roll-out is within 3 yd at 20–36° (rms 2 yd) and runs about 11 yd long at\n";
    s += "  12–18°.\n";
    s += "- **Shallow rebounds.** For landings flatter than about 35° the rebound is 0.07–0.15 livelier\n";
    s += "  than measured (Biber et al. 2023; Cross 2021); see the ey bands above.\n";
    s += "- **`GREEN` speed range.** `GREEN`'s friction is fitted to landings at 17–27 m/s. One green\n";
    s += "  measured at 9–10 m/s (Cross 2025, Table 1) slides at about 0.29; a chip with 3,300 rpm of\n";
    s += "  backspin leaves that green 1.3 m/s slower than this preset gives.\n";
    s += "- **Rough and sand** have no measured impacts or rolls; see the README.\n";
    s += "- **Aerodynamics** are not covered here: their fit quality is the \"Fit quality\" table of\n";
    s += "  [AERO.md](AERO.md), scored on launch-monitor and reference sets that are not distributed.\n";
    s
}

fn main() {
    let docs = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("docs");
    let csv_path = docs.join("fit-summary.csv");
    let mut rows = fixture_rows();
    let old = std::fs::read_to_string(&csv_path).unwrap_or_default();
    rows.extend(
        parse_csv(&old)
            .into_iter()
            .filter(|r| r.status.starts_with("summary")),
    );
    sort_rows(&mut rows);
    std::fs::write(&csv_path, write_csv(&rows)).expect("write docs/fit-summary.csv");
    std::fs::write(docs.join("FIT.md"), render(&rows)).expect("write docs/FIT.md");
    println!(
        "wrote {} rows to docs/FIT.md and docs/fit-summary.csv",
        rows.len()
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The numbers in a metrics string, each with its number of decimals.
    fn numbers(s: &str) -> Vec<(f64, i32)> {
        let s = s.replace('−', "-");
        let b: Vec<char> = s.chars().collect();
        let mut out = Vec::new();
        let mut i = 0;
        while i < b.len() {
            let signed =
                (b[i] == '-' || b[i] == '+') && b.get(i + 1).is_some_and(char::is_ascii_digit);
            if b[i].is_ascii_digit() || signed {
                let start = i;
                i += 1;
                while i < b.len() && (b[i].is_ascii_digit() || b[i] == '.') {
                    i += 1;
                }
                let tok: String = b[start..i].iter().collect();
                let dp = tok
                    .split_once('.')
                    .map_or(0, |(_, d)| i32::try_from(d.len()).unwrap_or(0));
                out.push((tok.parse().unwrap_or(f64::NAN), dp));
            } else {
                i += 1;
            }
        }
        out
    }

    #[test]
    fn fixture_rows_match_the_published_summary() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/fit-summary.csv");
        let published = parse_csv(&std::fs::read_to_string(path).expect("docs/fit-summary.csv"));
        let fresh = fixture_rows();
        let stored: Vec<_> = published
            .iter()
            .filter(|r| r.status.starts_with("fixture"))
            .collect();
        assert_eq!(
            stored.len(),
            fresh.len(),
            "fixture rows: regenerate with `cargo run --example fit_doc`"
        );
        for r in &fresh {
            let p = stored
                .iter()
                .find(|p| p.id == r.id)
                .unwrap_or_else(|| panic!("{} missing: regenerate docs/FIT.md", r.id));
            assert_eq!(r.status, p.status, "{}: regenerate docs/FIT.md", r.id);
            let (a, b) = (numbers(&r.metrics), numbers(&p.metrics));
            assert_eq!(a.len(), b.len(), "{}: {} vs {}", r.id, r.metrics, p.metrics);
            for ((x, _), (y, dp)) in a.iter().zip(&b) {
                let unit = 10f64.powi(-dp);
                assert!(
                    (x - y).abs() <= unit * 1.000_001,
                    "{}: {x} vs published {y}; regenerate with `cargo run --example fit_doc`",
                    r.id
                );
            }
        }
    }
}
