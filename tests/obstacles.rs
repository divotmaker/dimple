//! Bounces off obstacles (walls, cart paths, a tree trunk) with
//! `ground::impact`, and a shot that continues after one with
//! `simulate_shot_from`.
//!
//! The surfaces here are illustrative values, not measured ones.

use dimple::ground::{
    Contact, GroundModel, GroundPoint, GroundResult, GroundSample, Phase, Restitution, impact,
};
use dimple::{
    AeroModel, BallState, ShotInput, Surface, Tour, Vec3, simulate_shot, simulate_shot_from, units,
};

/// A hard surface (a cart path, a wall): rigid and flat (`compliance` 0),
/// a constant restitution. Illustrative values.
const HARD: Surface = Surface {
    restitution: Restitution::Constant(0.6),
    mu_slide: 0.5,
    compliance: 0.0,
    ..Surface::FAIRWAY
};

/// Bark. Illustrative values.
const BARK: Surface = Surface {
    restitution: Restitution::Constant(0.4),
    mu_slide: 0.6,
    compliance: 0.0,
    ..Surface::FAIRWAY
};

const BALL: &Tour = &Tour::SEA_LEVEL;

fn hit(surface: Surface, normal: Vec3) -> GroundSample {
    GroundSample {
        height: 0.0,
        normal,
        surface,
    }
}

fn energy(s: &BallState) -> f64 {
    let r = BALL.radius_m();
    0.5 * s.vel.length_squared() + 0.5 * BALL.inertia_factor() * r * r * s.spin.length_squared()
}

/// A ball arriving at `mph`, `descent` degrees below horizontal, moving
/// downrange, with backspin `rpm`.
fn arriving(mph: f64, descent: f64, rpm: f64) -> BallState {
    let v = units::mph_to_mps(mph);
    let d = units::deg_to_rad(descent);
    BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel: Vec3::new(v * d.cos(), -v * d.sin(), 0.0),
        spin: Vec3::new(0.0, 0.0, units::rpm_to_rads(rpm)),
    }
}

#[test]
fn backspin_checks_a_ball_off_a_cart_path() {
    // A full wedge landing, with and without its backspin.
    let up = Vec3::new(0.0, 1.0, 0.0);
    let spun = arriving(49.0, 51.0, 9100.0);
    let plain = arriving(49.0, 51.0, 0.0);
    let spun_out = impact(&spun, &hit(HARD, up), BALL);
    let plain_out = impact(&plain, &hit(HARD, up), BALL);
    assert!(
        spun_out.vel.x < plain_out.vel.x,
        "{} vs {}",
        spun_out.vel.x,
        plain_out.vel.x
    );
    // The friction takes the backspin out as it slows the ball: the ball
    // grips and leaves rolling forward, more slowly, and so with less spin,
    // than the ball that arrived without spin.
    assert!(spun_out.spin.z < 0.0, "{}", spun_out.spin.z);
    assert!(
        spun_out.spin.length() < plain_out.spin.length(),
        "{:?} vs {:?}",
        spun_out.spin,
        plain_out.spin
    );
    // The bounce height does not depend on the spin on a rigid flat surface.
    assert!((spun_out.vel.y - plain_out.vel.y).abs() < 1e-12);
    assert!(energy(&spun_out) <= energy(&spun));
}

#[test]
fn sidespin_kicks_a_ball_sideways_off_a_wall() {
    // A ball flying into a wall that faces the tee, with sidespin about the
    // vertical. Positive sidespin (a fade's, curving right) is about −Y.
    let wall = hit(HARD, Vec3::new(-1.0, 0.0, 0.0));
    let ball = |sidespin_rpm: f64| BallState {
        time: 0.0,
        pos: Vec3::ZERO,
        vel: Vec3::new(20.0, -3.0, 0.0),
        spin: Vec3::new(0.0, -units::rpm_to_rads(sidespin_rpm), 0.0),
    };
    let fade = impact(&ball(2000.0), &wall, BALL);
    let draw = impact(&ball(-2000.0), &wall, BALL);
    let none = impact(&ball(0.0), &wall, BALL);
    assert!(
        fade.vel.x < 0.0 && draw.vel.x < 0.0,
        "it rebounds off the wall"
    );
    // The spinning ball rolls along the wall: a fade's spin kicks it left
    // (−Z), a draw's right, by the same amount; no spin, no kick.
    assert!(fade.vel.z < -0.5, "{}", fade.vel.z);
    assert!((fade.vel.z + draw.vel.z).abs() < 1e-12);
    assert!(none.vel.z.abs() < 1e-12);
    assert!(energy(&fade) <= energy(&ball(2000.0)));
}

#[test]
fn without_spin_a_wall_returns_the_normal_velocity_times_the_restitution() {
    let e = 0.6;
    let normals = [
        Vec3::new(-1.0, 0.0, 0.0),             // a wall facing the tee
        Vec3::new(-0.6, 0.0, -0.8),            // a trunk, hit off center
        Vec3::new(0.3, -0.9, 0.3).normalize(), // the underside of an overhang
        Vec3::new(0.0, 1.0, 0.0),              // a cart path
    ];
    for n in normals {
        // Head on: straight back, slowed by the restitution.
        let head_on = BallState {
            time: 0.0,
            pos: Vec3::ZERO,
            vel: n * -15.0,
            spin: Vec3::ZERO,
        };
        let out = impact(&head_on, &hit(HARD, n), BALL);
        assert!(
            (out.vel - n * (15.0 * e)).length() < 1e-12,
            "{n:?}: {:?}",
            out.vel
        );
        assert_eq!(out.spin, Vec3::ZERO);

        // At an angle: the normal component is reversed and scaled by e;
        // friction acts only in the surface's plane.
        let oblique = BallState {
            vel: n * -15.0 + Vec3::new(1.0, 2.0, 3.0).cross(n),
            ..head_on
        };
        let out = impact(&oblique, &hit(HARD, n), BALL);
        let (v_in, v_out) = (oblique.vel.dot(n), out.vel.dot(n));
        assert!((v_out + e * v_in).abs() < 1e-12, "{n:?}: {v_in} → {v_out}");
    }
}

/// Flat [`Surface::FAIRWAY`] with one tree trunk: a vertical cylinder of
/// bark rising from the ground.
struct Tree {
    x: f64,
    z: f64,
    radius: f64,
}

impl Tree {
    /// Horizontal distance from the bark: positive outside the trunk.
    fn gap(&self, p: Vec3) -> f64 {
        (p.x - self.x).hypot(p.z - self.z) - self.radius
    }

    /// The bark's outward normal nearest `p`.
    fn outward(&self, p: Vec3) -> Vec3 {
        Vec3::new(p.x - self.x, 0.0, p.z - self.z).normalize()
    }
}

impl GroundModel for Tree {
    fn height(&self, _x: f64, _z: f64) -> f64 {
        0.0
    }

    /// The first of the ground crossing and the trunk crossing in the step
    /// from `a` to `b`.
    fn landing(&self, a: &BallState, b: &BallState) -> Option<BallState> {
        let ground = (a.pos.y >= 0.0 && b.pos.y <= 0.0).then(|| {
            if a.pos.y > b.pos.y {
                a.pos.y / (a.pos.y - b.pos.y)
            } else {
                0.0
            }
        });
        let trunk = (self.gap(a.pos) >= 0.0 && self.gap(b.pos) < 0.0).then(|| {
            // Bisect for the crossing, keeping the end outside the trunk.
            let (mut lo, mut hi) = (0.0, 1.0);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if self.gap(a.pos.lerp(b.pos, mid)) >= 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            lo
        });
        let (frac, on_ground) = match (ground, trunk) {
            (Some(g), Some(t)) if t < g => (t, false),
            (Some(g), _) => (g, true),
            (None, Some(t)) => (t, false),
            (None, None) => return None,
        };
        let mut pos = a.pos.lerp(b.pos, frac);
        if on_ground {
            pos.y = 0.0;
        }
        Some(BallState {
            time: a.time + (b.time - a.time) * frac,
            pos,
            vel: a.vel.lerp(b.vel, frac),
            spin: a.spin.lerp(b.spin, frac),
        })
    }

    /// On the ground: the fairway's physics. On the trunk: bounce off the
    /// bark and fly on with `simulate_shot_from`, which comes back here for
    /// the next landing.
    fn settle(&self, landing: &BallState, model: &dyn AeroModel) -> GroundResult {
        let on_trunk = landing.pos.y > 0.0 && self.gap(landing.pos).abs() < 1e-9;
        if !on_trunk {
            return Surface::FAIRWAY.settle(landing, model);
        }
        let bark = GroundSample {
            height: 0.0,
            normal: self.outward(landing.pos),
            surface: BARK,
        };
        let rebound = impact(landing, &bark, model);
        let after = simulate_shot_from(&rebound, model, self).ground;

        // A GroundResult holds contacts and sampled ground points; the
        // flight from the trunk to the next landing is kept only as its
        // first point (the rebound) and the next contact.
        let mut contacts = vec![Contact::new(landing, &bark, &rebound)];
        contacts.extend(after.contacts);
        let mut points = vec![GroundPoint::new(&rebound, Phase::Hop)];
        points.extend(after.points);
        GroundResult::new(contacts, points, after.rest, after.termination)
    }
}

#[test]
#[allow(
    clippy::float_cmp,
    reason = "a ground landing is placed exactly on the ground"
)]
fn a_shot_into_a_tree_trunk_bounces_back_and_comes_to_rest_short_of_it() {
    let tree = Tree {
        x: 40.0,
        z: 0.0,
        radius: 0.2,
    };
    let input = ShotInput::driver();
    let open = simulate_shot(&input, BALL, &Surface::FAIRWAY);
    let shot = simulate_shot(&input, BALL, &tree);

    // The flight ends at the bark, in the air, on the face toward the tee.
    let trunk_hit = shot.ground.contacts[0];
    assert!(trunk_hit.pos.y > 1.0, "{:?}", trunk_hit.pos);
    assert!((trunk_hit.pos.x - 39.8).abs() < 1e-6);
    assert!((trunk_hit.normal - Vec3::new(-1.0, 0.0, 0.0)).length() < 1e-9);
    assert!((shot.flight.carry_m - 39.8).abs() < 1e-6);
    // It rebounds toward the tee, lands again and comes to rest short of
    // the tree (this hard a drive runs back past the tee).
    assert!(trunk_hit.vel_out.x < 0.0);
    assert!(shot.ground.contacts.len() > 1);
    assert_eq!(shot.ground.contacts[1].pos.y, 0.0);
    let rest = shot.ground.rest.pos;
    assert!(rest.x < 38.0, "{rest:?}");
    assert!(shot.total_m < open.total_m / 4.0);
    // Time runs on through the bounce.
    assert!(shot.ground.contacts[1].time > trunk_hit.time);

    // A shot started wide of the trunk never touches it.
    let wide = simulate_shot(
        &ShotInput {
            launch_azimuth_deg: 3.0,
            ..input
        },
        BALL,
        &tree,
    );
    assert_eq!(wide.ground.contacts[0].pos.y, 0.0);
    assert!(wide.total_m > open.total_m * 0.95);
}
