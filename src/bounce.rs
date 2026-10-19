use glam::DVec3;

use crate::units::G;

/// Golf ball radius (m). USGA min diameter = 42.67mm.
const BALL_RADIUS: f64 = 0.021_335;

/// Hop height (m) below which the ball stops bouncing and rolls: 5 mm,
/// Penner's (2002) bounce-to-roll transition. A ball keeps bouncing while
/// its rebound would rise at least this high, `v_y² / (2g) ≥ 5 mm`
/// (`v_y` ≥ ~0.313 m/s). Penner tested transition heights from 1 to 10 mm
/// and found only a minor effect on the run.
pub(crate) const MIN_HOP_HEIGHT: f64 = 0.005;

/// Safety bound on impacts per bounce sequence. Every hop rises at least
/// [`MIN_HOP_HEIGHT`] and loses energy, so real shots end by the hop-height
/// rule far below this; the bound only guards against degenerate surfaces
/// (e.g. `Restitution::Constant(1.0)` with no friction).
pub(crate) const MAX_IMPACTS: usize = 200;

/// Penner (2002) eq. 8 reference contact angle θc (rad): 15.4°, his fit to
/// Haake's green impact.
const PENNER_THETA_C: f64 = 15.4 * std::f64::consts::PI / 180.0;
/// Impact speed (m/s) of Haake's measured green impact.
const PENNER_REF_SPEED: f64 = 18.6;
/// Incidence angle from vertical (rad) of Haake's measured green impact: 44.4°.
const PENNER_REF_ANGLE: f64 = 44.4 * std::f64::consts::PI / 180.0;

/// Normal coefficient of restitution between ball and turf.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Restitution {
    /// Constant normal coefficient of restitution (0..1).
    Constant(f64),
    /// Penner (2002) eq. 5, a speed-dependent fit to measured ball–turf
    /// restitution:
    /// e = 0.510 − 0.0375·v + 0.000903·v² for v ≤ 20 m/s, 0.120 above
    /// (v = normal impact speed against the effective contact plane, m/s).
    Penner2002,
}

impl Restitution {
    /// Coefficient of restitution at the given normal impact speed (m/s).
    ///
    /// The speed is taken as a magnitude; its sign is ignored.
    #[must_use]
    pub fn coefficient_at(&self, normal_speed: f64) -> f64 {
        match *self {
            Self::Constant(e) => e,
            Self::Penner2002 => {
                let v = normal_speed.abs();
                if v <= 20.0 {
                    0.510 - 0.0375 * v + 0.000_903 * v * v
                } else {
                    0.120
                }
            }
        }
    }
}

/// Surface parameters needed for a bounce computation.
#[derive(Debug, Clone, Copy)]
pub struct BounceSurface {
    /// Normal coefficient of restitution.
    pub restitution: Restitution,
    /// Sliding friction coefficient.
    pub mu: f64,
    /// Turf compliance: dimensionless multiplier on Penner's (2002) eq. 8
    /// effective contact angle, θc = 15.4° · (v / 18.6 m/s) · (φ / 44.4°),
    /// with v the impact speed and φ the incidence angle from vertical.
    /// `1.0` is Penner's green; `0.0` is a rigid flat surface.
    pub compliance: f64,
}

impl BounceSurface {
    /// Standard fairway: Penner's restitution and friction with compliance
    /// 0.9, slightly below the green's so rebounds come off slightly firmer.
    /// The value is chosen, not measured; the run is sensitive to it (see
    /// `docs/AERO.md`).
    pub const FAIRWAY: Self = Self {
        restitution: Restitution::Penner2002,
        mu: 0.40,
        compliance: 0.9,
    };

    /// Standard green: Penner's (2002) green model.
    pub const GREEN: Self = Self {
        restitution: Restitution::Penner2002,
        mu: 0.40,
        compliance: 1.0,
    };

    /// Effective contact angle β (rad) for an impact with horizontal speed
    /// `v_h` and downward speed `v_n` (both ≥ 0).
    ///
    /// `compliance` times Penner's eq. 8 θc, capped at the incidence angle
    /// from vertical: at the cap the ball strikes the effective plane
    /// head-on, and a larger tilt would turn its tangential velocity
    /// backward. Zero when there is no horizontal motion (the tilt direction
    /// is undefined).
    fn contact_angle(&self, v_h: f64, v_n: f64) -> f64 {
        if v_h <= 0.01 {
            return 0.0;
        }
        let speed = v_h.hypot(v_n);
        let incidence = v_h.atan2(v_n);
        let theta_c = PENNER_THETA_C * (speed / PENNER_REF_SPEED) * (incidence / PENNER_REF_ANGLE);
        (self.compliance * theta_c).min(incidence)
    }
}

/// Result of a single bounce computation.
#[derive(Debug, Clone, Copy)]
pub struct BounceResult {
    /// Rebound velocity (m/s) in world space.
    pub vel: DVec3,
    /// Angular velocity after bounce (rad/s), about the horizontal axis
    /// across the plane of travel, signed relative to downrange (+X)
    /// travel: positive is the forward-rolling sense for a ball moving
    /// downrange (topspin), negative opposes it (backspin). A ball rolling
    /// back toward the tee has negative `omega`. This is the convention
    /// [`simulate_rollout`](crate::rollout::simulate_rollout) takes; it
    /// equals minus the `backspin` argument of [`compute_bounce`].
    pub omega: f64,
    /// Sidespin after bounce (rad/s, positive = curves right).
    ///
    /// Bookkeeping only: the bounce decays it by a constant factor and it
    /// never affects [`vel`](Self::vel). Rollout does not take sidespin, so
    /// it has no lateral effect after landing; the ball keeps the lateral
    /// velocity it had at impact.
    pub sidespin: f64,
    /// Whether the ball keeps bouncing: its rebound rises at least 5 mm
    /// (Penner's 2002 bounce-to-roll transition height). Otherwise it rolls.
    pub still_bouncing: bool,
}

/// Compute the result of a golf ball bouncing on a surface.
///
/// Uses the oblique impact model from Penner (2002) and Cross (2002):
/// - Decompose velocity into normal and tangential components
/// - Apply Penner's effective contact angle (turf compliance)
/// - Normal: COR-based rebound
/// - Tangential: slip vs grip depending on friction and spin
///
/// Coordinate frame: +X downrange, +Y up, +Z right of the target line
/// (right-handed).
/// Ground plane is Y=0, surface normal is +Y.
///
/// `vel`: incoming velocity (m/s), must have negative Y component (falling).
/// `backspin`: incoming backspin (rad/s), signed relative to downrange (+X)
/// travel: positive opposes forward rolling for a ball moving downrange,
/// which is the flight's `landing_backspin_rads`. For a ball moving back
/// toward the tee, positive `backspin` is topspin relative to its travel.
/// `sidespin`: incoming sidespin (rad/s, positive = curves right).
/// `surface`: surface parameters.
#[must_use]
pub fn compute_bounce(
    vel: DVec3,
    backspin: f64,
    sidespin: f64,
    surface: &BounceSurface,
) -> BounceResult {
    // Spin is signed relative to downrange travel, but the impact is solved
    // along the ball's own direction of horizontal travel. `track_sign` maps
    // between the two: -1 when the ball moves back toward the tee.
    let track_sign = if vel.x < 0.0 { -1.0 } else { 1.0 };

    // Spin along the direction of travel: positive = topspin (contact point
    // moves backward, promotes rolling), negative = backspin.
    let omega_in = -backspin * track_sign;

    let v_n_abs = vel.y.abs(); // normal (vertical) impact speed
    let v_x = vel.x; // downrange
    let v_z = vel.z; // lateral
    let v_h = (v_x * v_x + v_z * v_z).sqrt(); // horizontal speed

    // Penner's effective contact angle: the ball deforms the turf and meets
    // the rising lip of its own crater, so it bounces as if off a rigid
    // plane inclined by beta toward the incoming ball. In the vertical plane
    // of travel (along-track horizontal, vertical) the effective outward
    // normal is n' = (-sin β, cos β) and the effective tangent is
    // t' = (cos β, sin β); this steepens the rebound.
    let beta = surface.contact_angle(v_h, v_n_abs);

    let cos_b = beta.cos();
    let sin_b = beta.sin();

    // Incoming speed along -n' and along t' (Penner 2002, eq. 6).
    // The clamp only absorbs rounding at the beta cap.
    let v_n_eff = v_n_abs * cos_b + v_h * sin_b;
    let v_h_eff = (v_h * cos_b - v_n_abs * sin_b).max(0.0);

    // Normal rebound. Restitution depends on the normal speed against the
    // effective plane (Penner 2002, eq. 5 uses |v_iy'|).
    let e_n = surface.restitution.coefficient_at(v_n_eff);
    let v_n_out = e_n * v_n_eff;

    // Tangential: contact point velocity determines slip vs grip.
    //
    // Contact point velocity: v_c = v_tangential - R * omega
    // With backspin (omega < 0): v_c = v_h_eff + R*|backspin| (slides forward)
    // With topspin (omega > 0): v_c = v_h_eff - R*topspin (closer to rolling)
    let v_c = v_h_eff - BALL_RADIUS * omega_in;

    // Normal impulse (mass-normalized): J_n/m = (1 + e_n) * v_n_eff
    let j_n = (1.0 + e_n) * v_n_eff;

    // Maximum friction impulse: J_t_max/m = mu * J_n/m
    let j_t_max = surface.mu * j_n;

    // Impulse needed to reach rolling condition (v_c_out = 0):
    // j_t_roll = |v_c| * (2/7) for a solid sphere (alpha = 2/5)
    let j_t_roll = (2.0 / 7.0) * v_c.abs();

    let (v_h_out, omega_out) = if j_t_roll <= j_t_max {
        // GRIP: friction is sufficient to reach rolling condition.
        // Ball leaves with v_c = 0 (rolling): v_out = R * omega_out.
        //
        // v_out = (5/7)*v_h_eff + (2/7)*R*omega_in
        // With backspin (omega_in < 0): v_out is reduced.
        // With enough backspin: v_out goes negative → spin-back.
        let v_out = (5.0 / 7.0) * v_h_eff + (2.0 / 7.0) * BALL_RADIUS * omega_in;
        let omega = v_out / BALL_RADIUS;
        (v_out, omega)
    } else {
        // SLIP: ball slides throughout bounce. Apply maximum friction impulse.
        let sign_vc = v_c.signum();
        let v_out = v_h_eff - j_t_max * sign_vc;
        let omega = omega_in + (5.0 / (2.0 * BALL_RADIUS)) * j_t_max * sign_vc;
        (v_out, omega)
    };

    // Reconstruct world-space velocity: v_out = v_n_out·n' + v_h_out·t'
    // (Penner 2002, eq. 7). v_h_out is signed (negative = spin-back), so the
    // along-track result is signed too: positive continues forward.
    let v_h_world = v_h_out * cos_b - v_n_out * sin_b;
    // A strong spin-back can point the rebound into the ground. The ball
    // cannot leave downward, so it stays on the surface and rolls with its
    // along-track velocity.
    let v_n_world = (v_n_out * cos_b + v_h_out * sin_b).max(0.0);

    // Distribute horizontal velocity back into X and Z components
    // preserving the original direction ratio.
    let (out_x, out_z) = if v_h > 1e-10 {
        (v_h_world * v_x / v_h, v_h_world * v_z / v_h)
    } else {
        // Vertical drop: no horizontal component
        (0.0, 0.0)
    };

    // Sidespin: largely unaffected by bounce. Slight decay.
    let sidespin_out = sidespin * 0.9;

    let still_bouncing = v_n_world * v_n_world >= 2.0 * G * MIN_HOP_HEIGHT;

    BounceResult {
        vel: DVec3::new(out_x, v_n_world, out_z),
        // Back from the direction of travel to downrange-signed spin.
        omega: omega_out * track_sign,
        sidespin: sidespin_out,
        still_bouncing,
    }
}

/// Compute all bounces until the ball transitions to rolling.
///
/// Returns the sequence of bounce results and the final state for rollout.
/// The sequence stops at the first rebound that would rise less than 5 mm
/// (Penner 2002); a safety bound of 200 impacts, never reached by real
/// shots, guards degenerate surfaces. Typical full shots make tens of small
/// hops before rolling.
/// Between bounces, the ball follows a parabolic arc that ignores aero
/// forces. Drag is not negligible on a drive's first hops (rebounds of up
/// to ~8 m/s with over a second in the air); ignoring it lengthens the run
/// somewhat.
///
/// `vel`: velocity at first ground contact (from flight simulation).
/// `backspin`: spin at landing (rad/s), signed relative to downrange travel
/// as in [`compute_bounce`]: positive = backspin for a ball moving downrange.
/// `sidespin`: sidespin at landing (rad/s, positive = curves right).
/// `surface`: surface parameters (assumed constant across bounces).
#[must_use]
pub fn bounce_sequence(
    vel: DVec3,
    backspin: f64,
    sidespin: f64,
    surface: &BounceSurface,
) -> Vec<BounceResult> {
    let mut results = Vec::new();
    let mut v = vel;
    let mut bs = backspin;
    let mut ss = sidespin;

    for _ in 0..MAX_IMPACTS {
        let result = compute_bounce(v, bs, ss, surface);
        let done = !result.still_bouncing;
        v = result.vel;
        bs = -result.omega;
        ss = result.sidespin;
        results.push(result);
        if done {
            break;
        }
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units;

    fn rpm_to_rads(rpm: f64) -> f64 {
        units::rpm_to_rads(rpm)
    }

    // ══════════════════════════════════════════════════════════════════
    //  STRUCTURAL TESTS — invariants for any surface
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn rebound_goes_up() {
        let vel = DVec3::new(0.0, -10.0, 0.0);
        let r = compute_bounce(vel, 0.0, 0.0, &BounceSurface::FAIRWAY);
        assert!(r.vel.y > 0.0, "rebound should go up: vy={:.3}", r.vel.y);
    }

    #[test]
    fn energy_decreases() {
        let vel = DVec3::new(20.0, -15.0, 0.0);
        let r = compute_bounce(vel, rpm_to_rads(5000.0), 0.0, &BounceSurface::FAIRWAY);
        assert!(
            r.vel.length() < vel.length(),
            "should lose energy: in={:.1}, out={:.1}",
            vel.length(),
            r.vel.length()
        );
    }

    #[test]
    fn no_lateral_without_sidespin() {
        let vel = DVec3::new(20.0, -15.0, 0.0);
        let r = compute_bounce(vel, rpm_to_rads(3000.0), 0.0, &BounceSurface::GREEN);
        assert!(
            r.vel.z.abs() < 1e-10,
            "no lateral without sidespin: vz={:.6}",
            r.vel.z
        );
    }

    #[test]
    fn vertical_drop_no_horizontal() {
        // Ball dropping straight down with no spin should bounce straight up.
        let vel = DVec3::new(0.0, -5.0, 0.0);
        let r = compute_bounce(vel, 0.0, 0.0, &BounceSurface::GREEN);
        assert!(r.vel.y > 0.0);
        assert!(r.vel.x.abs() < 1e-6, "no horizontal: vx={:.6}", r.vel.x);
    }

    // ══════════════════════════════════════════════════════════════════
    //  COR TESTS — normal rebound
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn green_rebound_steeper_and_slower_than_fairway() {
        // Both presets share Penner's restitution and friction; the green's
        // larger compliance tilts the effective contact plane further, which
        // steepens the rebound and takes more forward speed out of it.
        let vel = DVec3::new(15.0, -10.0, 0.0);
        let spin = rpm_to_rads(3000.0);
        let rf = compute_bounce(vel, spin, 0.0, &BounceSurface::FAIRWAY);
        let rg = compute_bounce(vel, spin, 0.0, &BounceSurface::GREEN);
        assert!(
            rg.vel.x < rf.vel.x,
            "green rebound slower: green vx={:.3}, fairway vx={:.3}",
            rg.vel.x,
            rf.vel.x
        );
        assert!(
            rg.vel.y.atan2(rg.vel.x) > rf.vel.y.atan2(rf.vel.x),
            "green rebound steeper: green=({:.3}, {:.3}), fairway=({:.3}, {:.3})",
            rg.vel.x,
            rg.vel.y,
            rf.vel.x,
            rf.vel.y
        );
    }

    #[test]
    fn penner2002_restitution_matches_eq5() {
        // Penner (2002) eq. 5: e = 0.510 − 0.0375·v + 0.000903·v² for
        // v ≤ 20 m/s, 0.120 above.
        let e = Restitution::Penner2002;
        for (v, expected) in [
            (0.0, 0.510),
            (5.0, 0.345_075),
            (10.0, 0.225_3),
            (15.0, 0.150_675),
            (20.0, 0.121_2),
        ] {
            assert!(
                (e.coefficient_at(v) - expected).abs() < 1e-12,
                "e({v}) = {:.6}, expected {expected}",
                e.coefficient_at(v)
            );
        }
        assert!((e.coefficient_at(-10.0) - e.coefficient_at(10.0)).abs() < f64::EPSILON);
    }

    #[test]
    fn penner2002_restitution_constant_above_20_mps() {
        let e = Restitution::Penner2002;
        for v in [20.001, 25.0, 40.0] {
            assert!(
                (e.coefficient_at(v) - 0.120).abs() < f64::EPSILON,
                "e({v}) = {}",
                e.coefficient_at(v)
            );
        }
        // The two branches meet within 0.0012 at 20 m/s.
        assert!((e.coefficient_at(20.0) - e.coefficient_at(20.001)).abs() < 0.002);
    }

    #[test]
    fn constant_restitution_ignores_speed() {
        let e = Restitution::Constant(0.3);
        assert!((e.coefficient_at(1.0) - 0.3).abs() < f64::EPSILON);
        assert!((e.coefficient_at(30.0) - 0.3).abs() < f64::EPSILON);
    }

    // ══════════════════════════════════════════════════════════════════
    //  SPIN TESTS — backspin effects
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn backspin_reduces_forward_speed() {
        let vel = DVec3::new(15.0, -10.0, 0.0);
        let r_low = compute_bounce(vel, rpm_to_rads(2000.0), 0.0, &BounceSurface::GREEN);
        let r_high = compute_bounce(vel, rpm_to_rads(8000.0), 0.0, &BounceSurface::GREEN);
        assert!(
            r_high.vel.x < r_low.vel.x,
            "high spin should reduce forward speed: high={:.3}, low={:.3}",
            r_high.vel.x,
            r_low.vel.x
        );
    }

    #[test]
    fn very_high_backspin_can_reverse() {
        // Wedge shot: steep angle, very high backspin on green → spin back.
        let speed = 25.0;
        let angle = units::deg_to_rad(50.0);
        let vel = DVec3::new(speed * angle.cos(), -speed * angle.sin(), 0.0);
        let backspin = rpm_to_rads(10000.0);
        let r = compute_bounce(vel, backspin, 0.0, &BounceSurface::GREEN);
        assert!(
            r.vel.x < 0.0,
            "very high backspin on green should spin back: vx={:.3}",
            r.vel.x
        );
    }

    #[test]
    fn driver_bounces_forward() {
        // Driver: shallow angle, low spin → bounces forward.
        let speed = 30.0;
        let angle = units::deg_to_rad(30.0);
        let vel = DVec3::new(speed * angle.cos(), -speed * angle.sin(), 0.0);
        let backspin = rpm_to_rads(2500.0);
        let r = compute_bounce(vel, backspin, 0.0, &BounceSurface::FAIRWAY);
        assert!(
            r.vel.x > 0.0,
            "driver should bounce forward: vx={:.3}",
            r.vel.x
        );
    }

    #[test]
    fn seven_iron_checks_on_green() {
        let speed = 25.0;
        let angle = units::deg_to_rad(45.0);
        let vel = DVec3::new(speed * angle.cos(), -speed * angle.sin(), 0.0);
        let backspin = rpm_to_rads(7000.0);
        let r = compute_bounce(vel, backspin, 0.0, &BounceSurface::GREEN);
        let r_nospin = compute_bounce(vel, 0.0, 0.0, &BounceSurface::GREEN);
        assert!(
            r.vel.x < r_nospin.vel.x,
            "7-iron with spin should check: spin vx={:.3}, nospin vx={:.3}",
            r.vel.x,
            r_nospin.vel.x
        );
    }

    // ══════════════════════════════════════════════════════════════════
    //  PENNER COMPLIANCE — steeper rebound on soft surfaces
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn compliance_steepens_rebound() {
        // Isolate compliance: same restitution and mu, only compliance differs.
        let vel = DVec3::new(20.0, -15.0, 0.0);
        let spin = rpm_to_rads(5000.0);

        let firm = BounceSurface {
            restitution: Restitution::Constant(0.35),
            mu: 0.40,
            compliance: 0.2,
        };
        let soft = BounceSurface {
            restitution: Restitution::Constant(0.35),
            mu: 0.40,
            compliance: 1.0,
        };

        let r_firm = compute_bounce(vel, spin, 0.0, &firm);
        let r_soft = compute_bounce(vel, spin, 0.0, &soft);

        // Higher compliance → steeper rebound (higher vy/vx ratio)
        let angle_firm = r_firm.vel.y.atan2(r_firm.vel.x.abs().max(0.001));
        let angle_soft = r_soft.vel.y.atan2(r_soft.vel.x.abs().max(0.001));
        assert!(
            angle_soft > angle_firm,
            "soft surface steeper: firm={:.1}°, soft={:.1}°",
            units::rad_to_deg(angle_firm),
            units::rad_to_deg(angle_soft)
        );
    }

    #[test]
    fn compliance_rotates_mirror_rebound_by_two_beta() {
        // Elastic, frictionless impact at 45°: the effective plane is tilted
        // by Penner's θc = 15.4° · (v / 18.6 m/s) · (φ / 44.4°) toward the
        // ball, so the specular rebound off it is the flat-ground mirror
        // direction rotated up by 2β, at full speed.
        let vel = DVec3::new(10.0, -10.0, 0.0);
        let surface = BounceSurface {
            restitution: Restitution::Constant(1.0),
            mu: 0.0,
            compliance: 1.0,
        };
        let r = compute_bounce(vel, 0.0, 0.0, &surface);

        let speed = vel.length();
        let incidence = std::f64::consts::FRAC_PI_4;
        let beta = units::deg_to_rad(15.4) * (speed / 18.6) * (incidence / units::deg_to_rad(44.4));
        let mirror_angle = std::f64::consts::FRAC_PI_4;
        let rebound_angle = r.vel.y.atan2(r.vel.x);
        assert!(
            (rebound_angle - (mirror_angle + 2.0 * beta)).abs() < 1e-12,
            "rebound {:.4}° should be mirror {:.4}° + 2β {:.4}°",
            units::rad_to_deg(rebound_angle),
            units::rad_to_deg(mirror_angle),
            units::rad_to_deg(2.0 * beta)
        );
        assert!(
            (r.vel.length() - vel.length()).abs() < 1e-12,
            "elastic frictionless impact keeps speed: in={:.6}, out={:.6}",
            vel.length(),
            r.vel.length()
        );
        assert!(r.vel.x > 0.0 && r.vel.z.abs() < 1e-12);
    }

    #[test]
    fn no_compliance_is_exact_mirror() {
        let vel = DVec3::new(10.0, -10.0, 0.0);
        let surface = BounceSurface {
            restitution: Restitution::Constant(1.0),
            mu: 0.0,
            compliance: 0.0,
        };
        let r = compute_bounce(vel, 0.0, 0.0, &surface);
        assert!(
            (r.vel - DVec3::new(10.0, 10.0, 0.0)).length() < 1e-12,
            "compliance=0 should mirror: got ({:.6}, {:.6}, {:.6})",
            r.vel.x,
            r.vel.y,
            r.vel.z
        );
    }

    #[test]
    fn compliance_tilt_capped_at_incidence_angle() {
        // 10 × θc ≈ 119° would overshoot; β is capped at the 45° incidence
        // angle from vertical, so the ball meets the effective plane head-on
        // and an elastic frictionless impact sends it straight back.
        let vel = DVec3::new(10.0, -10.0, 0.0);
        let surface = BounceSurface {
            restitution: Restitution::Constant(1.0),
            mu: 0.0,
            compliance: 10.0,
        };
        let r = compute_bounce(vel, 0.0, 0.0, &surface);
        assert!(
            (r.vel + vel).length() < 1e-9,
            "head-on impact should reverse: got ({:.6}, {:.6}, {:.6})",
            r.vel.x,
            r.vel.y,
            r.vel.z
        );
    }

    #[test]
    fn spin_back_into_ground_rolls_instead_of_rebounding_down() {
        // Dead turf (e=0.05), strong tilt (2 × θc ≈ 0.41 rad) and extreme
        // backspin: the ball grips and leaves the effective plane moving
        // backward fast enough that the reconstructed rebound points below
        // horizontal. It must stay on the ground and roll back, never
        // rebound downward.
        let vel = DVec3::new(10.0, -10.0, 0.0);
        let surface = BounceSurface {
            restitution: Restitution::Constant(0.05),
            mu: 1.0,
            compliance: 2.0,
        };
        let r = compute_bounce(vel, 1500.0, 0.0, &surface);
        assert!(
            r.vel.y.abs() < f64::EPSILON,
            "no downward rebound: vy={:.6}",
            r.vel.y
        );
        assert!(!r.still_bouncing);
        assert!(r.vel.x < 0.0, "spins back: vx={:.3}", r.vel.x);
    }

    // ══════════════════════════════════════════════════════════════════
    //  PENNER (2002) VALIDATION — GREEN reproduces the published model
    // ══════════════════════════════════════════════════════════════════

    /// Impact at `speed` m/s, `from_vertical_deg` from the vertical, moving +X.
    fn impact(speed: f64, from_vertical_deg: f64) -> DVec3 {
        let phi = units::deg_to_rad(from_vertical_deg);
        DVec3::new(speed * phi.sin(), -speed * phi.cos(), 0.0)
    }

    #[test]
    fn penner_drive_first_bounce_heights_on_green() {
        // Penner (2002) §3.1: drives land at 26.8 m/s, 68.0° from vertical
        // with 224 rad/s backspin (35 m/s club-head speed) and at 26.2 m/s,
        // 46.2° with 315 rad/s (55 m/s); the first bounce rises 2.82 m and
        // 1.17 m respectively, and the ball leaves with topspin.
        for (speed, angle, backspin, height) in
            [(26.8, 68.0, 224.0, 2.82), (26.2, 46.2, 315.0, 1.17)]
        {
            let r = compute_bounce(impact(speed, angle), backspin, 0.0, &BounceSurface::GREEN);
            let peak = r.vel.y * r.vel.y / (2.0 * units::G);
            assert!(
                (peak - height).abs() < 0.03,
                "{speed} m/s at {angle}°: first bounce peak {peak:.3} m, Penner {height} m"
            );
            assert!(
                r.vel.x > 0.0 && r.omega > 0.0,
                "leaves forward with topspin: vx={:.3}, ω={:.1}",
                r.vel.x,
                r.omega
            );
        }
    }

    #[test]
    fn penner_nine_iron_bounces_backward_on_green() {
        // Penner (2002) §3.2: a 9-iron landing at 23.8 m/s, 34.0° from
        // vertical with 809 rad/s backspin (below the 903 rad/s needed to
        // reverse along the effective plane) still bounces backward once the
        // rebound is rotated back to the ground frame.
        let r = compute_bounce(impact(23.8, 34.0), 809.0, 0.0, &BounceSurface::GREEN);
        assert!(r.vel.x < 0.0, "bounces backward: vx={:.3}", r.vel.x);
        assert!(r.vel.y > 0.0 && r.still_bouncing);
    }

    // ══════════════════════════════════════════════════════════════════
    //  BOUNCE SEQUENCE
    // ══════════════════════════════════════════════════════════════════

    fn hop_height(b: &BounceResult) -> f64 {
        b.vel.y * b.vel.y / (2.0 * units::G)
    }

    /// Asserts that a sequence ended by the 5 mm hop-height rule: every
    /// rebound but the last rises at least 5 mm, the last rises less, and
    /// the safety bound was not reached.
    fn assert_ends_by_hop_height(bounces: &[BounceResult]) {
        assert!(
            bounces.len() < MAX_IMPACTS,
            "reached the {MAX_IMPACTS}-impact safety bound"
        );
        let (last, rest) = bounces.split_last().unwrap();
        for b in rest {
            assert!(b.still_bouncing && hop_height(b) >= MIN_HOP_HEIGHT);
        }
        assert!(
            !last.still_bouncing && hop_height(last) < MIN_HOP_HEIGHT,
            "last hop {:.2} mm should be below 5 mm",
            hop_height(last) * 1000.0
        );
    }

    #[test]
    fn bounce_sequence_terminates() {
        let vel = DVec3::new(20.0, -15.0, 0.0);
        let bounces = bounce_sequence(vel, rpm_to_rads(3000.0), 0.0, &BounceSurface::FAIRWAY);
        assert_ends_by_hop_height(&bounces);
    }

    #[test]
    fn bounce_sequence_stops_below_threshold() {
        // Straight drop at 15 m/s: no tilt, so each rebound is e(v)·v with
        // Penner's e rising toward 0.51 as the ball slows. Rebounds of
        // 2.26, 0.97, 0.46 and 0.23 m/s rise 260, 48, 10.8 and 2.6 mm: the
        // fourth is below the 5 mm transition height.
        let bounces = bounce_sequence(
            DVec3::new(0.0, -15.0, 0.0),
            0.0,
            0.0,
            &BounceSurface::FAIRWAY,
        );
        assert_eq!(bounces.len(), 4);
        assert_ends_by_hop_height(&bounces);
    }

    #[test]
    fn bounce_sequence_safety_bound() {
        // Elastic, frictionless, rigid: every rebound returns at full speed,
        // so only the safety bound ends the sequence.
        let surface = BounceSurface {
            restitution: Restitution::Constant(1.0),
            mu: 0.0,
            compliance: 0.0,
        };
        let bounces = bounce_sequence(DVec3::new(5.0, -5.0, 0.0), 0.0, 0.0, &surface);
        assert_eq!(bounces.len(), MAX_IMPACTS);
        assert!(bounces.iter().all(|b| b.still_bouncing));
    }

    #[test]
    fn bounce_sequence_height_decreases() {
        let vel = DVec3::new(15.0, -10.0, 0.0);
        let bounces = bounce_sequence(vel, rpm_to_rads(2000.0), 0.0, &BounceSurface::FAIRWAY);
        for w in bounces.windows(2) {
            assert!(
                w[1].vel.y <= w[0].vel.y + 0.01,
                "bounce height should decrease: {:.3} then {:.3}",
                w[0].vel.y,
                w[1].vel.y
            );
        }
    }

    #[test]
    fn soft_drop_ends_by_hop_height() {
        // Gentle, shallow landing on green: the tilt keeps turning forward
        // speed into small hops that shrink slowly, and the sequence ends
        // once a hop would rise less than 5 mm.
        let vel = DVec3::new(5.0, -1.0, 0.0);
        let bounces = bounce_sequence(vel, rpm_to_rads(500.0), 0.0, &BounceSurface::GREEN);
        assert_ends_by_hop_height(&bounces);
    }

    // ══════════════════════════════════════════════════════════════════
    //  SPIN FRAME — spin is signed relative to downrange travel
    // ══════════════════════════════════════════════════════════════════

    #[test]
    fn bounce_mirrors_for_backward_travel() {
        // Mirroring the impact in x (ball travelling back toward the tee)
        // mirrors the rebound velocity in x and flips the downrange-signed
        // spin, both on input and output.
        for (vel, backspin) in [
            (DVec3::new(5.0, -3.0, 0.0), -5.0 / BALL_RADIUS),
            (DVec3::new(15.0, -10.0, 2.0), 300.0),
            (DVec3::new(23.8, -16.0, 0.0), 900.0),
        ] {
            for surface in [BounceSurface::FAIRWAY, BounceSurface::GREEN] {
                let fwd = compute_bounce(vel, backspin, 0.0, &surface);
                let mirrored = DVec3::new(-vel.x, vel.y, vel.z);
                let bwd = compute_bounce(mirrored, -backspin, 0.0, &surface);
                let expected = DVec3::new(-fwd.vel.x, fwd.vel.y, fwd.vel.z);
                assert!(
                    (bwd.vel - expected).length() < 1e-9 && (bwd.omega + fwd.omega).abs() < 1e-9,
                    "mirror of ({:.3}, {:.3}, {:.3}) ω {:.1}: got ({:.3}, {:.3}, {:.3}) ω {:.1}",
                    fwd.vel.x,
                    fwd.vel.y,
                    fwd.vel.z,
                    fwd.omega,
                    bwd.vel.x,
                    bwd.vel.y,
                    bwd.vel.z,
                    bwd.omega
                );
            }
        }
    }

    #[test]
    fn ball_rolling_backward_keeps_rolling_backward() {
        // A ball moving back toward the tee with matching rolling spin
        // (negative downrange-signed omega, i.e. positive `backspin`) grips
        // and leaves still moving backward with spin in the
        // backward-rolling sense.
        let r = compute_bounce(
            DVec3::new(-5.0, -3.0, 0.0),
            5.0 / BALL_RADIUS,
            0.0,
            &BounceSurface::GREEN,
        );
        assert!(
            r.vel.x < -4.0 && r.vel.y > 0.0 && r.omega < -150.0,
            "leaves rolling backward: vx={:.3}, vy={:.3}, ω={:.1}",
            r.vel.x,
            r.vel.y,
            r.omega
        );
    }

    /// Net run (m) along x from the landing point: hops, then rollout.
    fn run_out(vel: DVec3, backspin: f64, surface: &BounceSurface) -> f64 {
        let bounces = bounce_sequence(vel, backspin, 0.0, surface);
        assert_ends_by_hop_height(&bounces);
        let hops: f64 = bounces
            .iter()
            .filter(|b| b.still_bouncing)
            .map(|b| b.vel.x * 2.0 * b.vel.y / units::G)
            .sum();
        let last = bounces.last().unwrap();
        let roll = crate::rollout::simulate_rollout(
            DVec3::new(hops, 0.0, 0.0),
            last.vel,
            last.omega,
            &crate::rollout::RolloutSurface::GREEN,
        );
        roll.final_pos.x
    }

    #[test]
    fn penner_nine_iron_spin_back_grows_with_backspin() {
        // Penner (2002) Fig. 7: the 9-iron landing (23.8 m/s, 34.0° from
        // vertical) on his green with 809, 1000 and 1200 rad/s backspin.
        // More backspin carries the ball further back from where it landed.
        let runs: Vec<f64> = [809.0, 1000.0, 1200.0]
            .iter()
            .map(|&w| run_out(impact(23.8, 34.0), w, &BounceSurface::GREEN))
            .collect();
        assert!(
            runs[0] <= 0.0 && runs[1] < runs[0] - 1.0 && runs[2] < runs[1] - 1.0,
            "backward run should grow with backspin: {runs:.2?} m"
        );
    }

    #[test]
    fn penner_nine_iron_firm_green_bounces_forward_then_spins_back() {
        // Penner (2002) Fig. 8: with μ = 0.25 the friction cannot check the
        // backspin on the first impact, so the ball bounces forward still
        // carrying backspin, then runs back; with 1000 and 1200 rad/s it
        // finishes behind its landing point.
        let firm = BounceSurface {
            mu: 0.25,
            ..BounceSurface::GREEN
        };
        for w in [809.0, 1000.0, 1200.0] {
            let first = compute_bounce(impact(23.8, 34.0), w, 0.0, &firm);
            assert!(
                first.vel.x > 0.0 && first.omega < 0.0,
                "{w} rad/s: first bounce forward with backspin: vx={:.3}, ω={:.1}",
                first.vel.x,
                first.omega
            );
        }
        let low = run_out(impact(23.8, 34.0), 1000.0, &firm);
        let high = run_out(impact(23.8, 34.0), 1200.0, &firm);
        assert!(
            low < 0.0 && high < low,
            "firm green spin-back: 1000 rad/s {low:.2} m, 1200 rad/s {high:.2} m"
        );
    }
}
