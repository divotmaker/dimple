use crate::trajectory::{FlightResult, ShotResult};
use crate::units;

/// Statistics for a single shot, in display-friendly units.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct ShotStats {
    /// Carry distance (yards).
    pub carry_yards: f64,
    /// Carry + rollout distance (yards). Equal to carry for flight-only stats.
    pub total_yards: f64,
    /// Total minus carry (yards). Zero for flight-only stats.
    pub rollout_yards: f64,
    /// Maximum height (feet).
    pub apex_feet: f64,
    /// Time from launch to landing (s).
    pub flight_time_s: f64,
    /// Lateral offset of the landing point from the target line (yards).
    /// Positive = right, negative = left, from the golfer's view.
    pub lateral_yards: f64,
    /// Angle below horizontal at landing (degrees).
    pub descent_angle_deg: f64,
}

impl ShotStats {
    /// Compute shot statistics from a flight result (no ground phase).
    #[must_use]
    pub fn from_flight(result: &FlightResult) -> Self {
        Self {
            carry_yards: result.carry_yards,
            total_yards: result.carry_yards,
            rollout_yards: 0.0,
            apex_feet: units::meters_to_feet(result.apex_m),
            flight_time_s: result.flight_time,
            lateral_yards: units::meters_to_yards(result.lateral_m),
            descent_angle_deg: Self::descent_angle(result),
        }
    }

    /// Compute shot statistics from a full shot result (flight + ground).
    #[must_use]
    pub fn from_shot(result: &ShotResult) -> Self {
        let rollout_yards = result.total_yards - result.flight.carry_yards;
        Self {
            carry_yards: result.flight.carry_yards,
            total_yards: result.total_yards,
            rollout_yards,
            apex_feet: units::meters_to_feet(result.flight.apex_m),
            flight_time_s: result.flight.flight_time,
            lateral_yards: units::meters_to_yards(result.flight.lateral_m),
            descent_angle_deg: Self::descent_angle(&result.flight),
        }
    }

    fn descent_angle(result: &FlightResult) -> f64 {
        if let Some(last) = result.points.last() {
            let vy = last.vel.y;
            let vx = last.vel.x;
            let vz = last.vel.z;
            let horizontal_speed = (vx * vx + vz * vz).sqrt();
            if horizontal_speed > 1e-10 {
                units::rad_to_deg((-vy).atan2(horizontal_speed))
            } else {
                90.0
            }
        } else {
            0.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aero::BallModel;
    use crate::environment::Environment;
    use crate::trajectory::{ShotInput, simulate_flight};

    #[test]
    fn shot_stats_from_seven_iron() {
        let result = simulate_flight(
            &ShotInput::seven_iron(),
            &BallModel::TOUR,
            &Environment::SEA_LEVEL,
        );
        let stats = ShotStats::from_flight(&result);
        assert!(stats.carry_yards > 100.0);
        assert!(stats.apex_feet > 0.0);
        assert!(stats.descent_angle_deg > 0.0);
    }

    #[test]
    fn shot_stats_lateral_positive_right() {
        let slice = ShotInput {
            sidespin_rpm: 1500.0,
            ..ShotInput::seven_iron()
        };
        let result = simulate_flight(&slice, &BallModel::TOUR, &Environment::SEA_LEVEL);
        let stats = ShotStats::from_flight(&result);
        assert!(
            stats.lateral_yards > 0.0,
            "slice lateral {:.1}",
            stats.lateral_yards
        );
        assert!((units::yards_to_meters(stats.lateral_yards) - result.lateral_m).abs() < 1e-9);
    }
}
