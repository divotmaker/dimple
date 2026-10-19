use glam::DVec3;

/// State vector for RK4 integration: position and velocity.
#[derive(Debug, Clone, Copy)]
pub struct State {
    pub pos: DVec3,
    pub vel: DVec3,
}

/// Advance state by one RK4 step.
///
/// `accel_fn(position, velocity, time) -> acceleration`
#[must_use]
pub fn step<F>(state: State, t: f64, dt: f64, accel_fn: F) -> State
where
    F: Fn(DVec3, DVec3, f64) -> DVec3,
{
    let k1_v = accel_fn(state.pos, state.vel, t);
    let k1_x = state.vel;

    let k2_v = accel_fn(
        state.pos + k1_x * (dt / 2.0),
        state.vel + k1_v * (dt / 2.0),
        t + dt / 2.0,
    );
    let k2_x = state.vel + k1_v * (dt / 2.0);

    let k3_v = accel_fn(
        state.pos + k2_x * (dt / 2.0),
        state.vel + k2_v * (dt / 2.0),
        t + dt / 2.0,
    );
    let k3_x = state.vel + k2_v * (dt / 2.0);

    let k4_v = accel_fn(state.pos + k3_x * dt, state.vel + k3_v * dt, t + dt);
    let k4_x = state.vel + k3_v * dt;

    State {
        pos: state.pos + (k1_x + 2.0 * k2_x + 2.0 * k3_x + k4_x) * (dt / 6.0),
        vel: state.vel + (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v) * (dt / 6.0),
    }
}

/// A state that [`step_state`] can integrate: a point in a vector space.
///
/// The time derivative of a state has the same type as the state itself.
pub trait Integrable: Copy {
    /// `self + rate · h`, component by component.
    #[must_use]
    fn add_scaled(&self, rate: &Self, h: f64) -> Self;
}

/// Advance an arbitrary state by one RK4 step.
///
/// `deriv(state, time) -> d(state)/dt`. [`step`] is the position/velocity
/// special case used for flight.
#[must_use]
pub fn step_state<S, F>(state: S, t: f64, dt: f64, deriv: F) -> S
where
    S: Integrable,
    F: Fn(&S, f64) -> S,
{
    let half = dt / 2.0;
    let k1 = deriv(&state, t);
    let k2 = deriv(&state.add_scaled(&k1, half), t + half);
    let k3 = deriv(&state.add_scaled(&k2, half), t + half);
    let k4 = deriv(&state.add_scaled(&k3, dt), t + dt);
    state
        .add_scaled(&k1, dt / 6.0)
        .add_scaled(&k2, dt / 3.0)
        .add_scaled(&k3, dt / 3.0)
        .add_scaled(&k4, dt / 6.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::G;

    #[test]
    fn vacuum_parabola() {
        // Launch at 45 degrees, 50 m/s, no drag
        let speed = 50.0;
        let angle = std::f64::consts::FRAC_PI_4;
        let state = State {
            pos: DVec3::ZERO,
            vel: DVec3::new(speed * angle.cos(), speed * angle.sin(), 0.0),
        };

        let dt = 0.001;
        let mut s = state;
        let mut t = 0.0;

        // Gravity only
        let gravity = |_pos: DVec3, _vel: DVec3, _t: f64| DVec3::new(0.0, -G, 0.0);

        while s.pos.y >= 0.0 || t < 0.1 {
            s = step(s, t, dt, gravity);
            t += dt;
        }

        // Analytical range: v^2 * sin(2*theta) / g
        let expected_range = speed * speed * (2.0 * angle).sin() / G;
        assert!(
            (s.pos.x - expected_range).abs() < 0.5,
            "range {:.1} vs expected {:.1}",
            s.pos.x,
            expected_range
        );
    }

    #[test]
    fn constant_acceleration() {
        let state = State {
            pos: DVec3::ZERO,
            vel: DVec3::ZERO,
        };
        let dt = 0.001;
        let a = DVec3::new(1.0, 0.0, 0.0);

        let mut s = state;
        for i in 0..1000 {
            s = step(s, i as f64 * dt, dt, |_, _, _| a);
        }

        let t = 1.0;
        let expected_pos = 0.5 * a * t * t;
        let expected_vel = a * t;

        assert!((s.pos - expected_pos).length() < 1e-6);
        assert!((s.vel - expected_vel).length() < 1e-6);
    }

    #[derive(Debug, Clone, Copy)]
    struct Scalar(f64);

    impl Integrable for Scalar {
        fn add_scaled(&self, rate: &Self, h: f64) -> Self {
            Self(self.0 + rate.0 * h)
        }
    }

    #[test]
    fn step_state_exponential_decay() {
        // dy/dt = -y: RK4 is fourth-order accurate, so 1000 steps of 1 ms
        // land within ~1e-13 of exp(-1).
        let mut y = Scalar(1.0);
        for i in 0..1000 {
            y = step_state(y, f64::from(i) * 0.001, 0.001, |s, _| Scalar(-s.0));
        }
        assert!((y.0 - (-1.0_f64).exp()).abs() < 1e-12, "y(1) = {}", y.0);
    }

    #[test]
    fn step_state_matches_step_for_constant_acceleration() {
        #[derive(Debug, Clone, Copy)]
        struct PosVel(DVec3, DVec3);
        impl Integrable for PosVel {
            fn add_scaled(&self, rate: &Self, h: f64) -> Self {
                Self(self.0 + rate.0 * h, self.1 + rate.1 * h)
            }
        }
        let g = DVec3::new(0.0, -G, 0.0);
        let a = step(
            State {
                pos: DVec3::ZERO,
                vel: DVec3::new(3.0, 4.0, 1.0),
            },
            0.0,
            0.25,
            |_, _, _| g,
        );
        let b = step_state(
            PosVel(DVec3::ZERO, DVec3::new(3.0, 4.0, 1.0)),
            0.0,
            0.25,
            |s, _| PosVel(s.1, g),
        );
        assert!((a.pos - b.0).length() < 1e-14 && (a.vel - b.1).length() < 1e-14);
    }
}
