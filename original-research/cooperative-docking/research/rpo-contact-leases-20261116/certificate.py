"""Scalar bounded-disturbance certificates and the clock-rate lower bound."""

from dataclasses import dataclass
from math import sqrt


@dataclass(frozen=True)
class Approach:
    gap_min: float
    gap_max: float
    speed_min: float
    speed_max: float
    lag_min: float
    lag_max: float
    disturbance: float
    accel_max: float
    contact_speed_min: float
    contact_speed_max: float


@dataclass(frozen=True)
class Certificate:
    acceleration: float
    acceleration_max: float
    earliest: float
    latest: float
    contact_speed_min: float
    contact_speed_max: float


def guaranteed_remaining(lease, elapsed, receiver_drift, chaser_drift, receiver_processing=0.0):
    """Lower bound in physical seconds, without either clock's epoch."""
    age_upper = max(0.0, elapsed / (1 - chaser_drift) - receiver_processing / (1 + receiver_drift))
    return lease / (1 + receiver_drift) - age_upper


def stopping_distance(speed, latency, brake, disturbance):
    speed = max(0.0, speed)
    return (speed * latency + 0.5 * disturbance * latency**2
            + (speed + disturbance * latency)**2 / (2 * (brake - disturbance)))


def displacement(speed, disturbance, acceleration, lag, time):
    return speed * time + 0.5 * (disturbance * time**2 + acceleration * max(time - lag, 0)**2)


def _crossing(gap, speed, disturbance, acceleration, lag, horizon):
    """Enclose the first root of this increasing displacement on [0, horizon]."""
    low, high = 0.0, horizon
    for _ in range(48):
        middle = (low + high) / 2
        if displacement(speed, disturbance, acceleration, lag, middle) < gap:
            low = middle
        else:
            high = middle
    return low, high


def certify(approach, horizon, *, times=True):
    """Find the minimum constant inward acceleration certified over horizon.

    Initial gap and speed lie in a box, additive acceleration disturbance is
    bounded pointwise, and only command onset has an uncertain lag.
    """
    c = approach
    d = c.disturbance
    if horizon <= 0 or c.gap_min <= 0 or c.speed_min <= d * horizon:
        return None
    lag_min, lag_max = min(c.lag_min, horizon), min(c.lag_max, horizon)
    slow_before_onset = max(0.0, displacement(c.speed_min, -d, 0, 0, lag_min))
    fast_before_onset = displacement(c.speed_max, d, 0, 0, lag_max)
    remaining_max = max(0.0, c.gap_max - slow_before_onset)

    deficit = c.gap_max - displacement(c.speed_min, -d, 0, 0, horizon)
    thrust_time = horizon - lag_max
    lower = 0.0
    if deficit > 0:
        if thrust_time == 0:
            return None
        lower = 2 * deficit / thrust_time**2
    # Minimize the lower energy envelope over the unknown initial gap.
    # It is piecewise affine with a single kink at fast_before_onset.
    energy_lines = []
    kink = min(c.gap_max, max(c.gap_min, fast_before_onset))
    for gap in (c.gap_min, kink, c.gap_max):
        intercept = c.speed_min**2 - 2 * d * gap
        slope = 2 * max(0.0, gap - fast_before_onset)
        energy_lines.append((intercept, slope))
        speed_deficit = c.contact_speed_min**2 - intercept
        if speed_deficit > 0:
            if slope == 0:
                return None
            lower = max(lower, speed_deficit / slope)
    energy_max = c.speed_max**2 + 2 * d * c.gap_max
    headroom = c.contact_speed_max**2 - energy_max
    if headroom < 0:
        return None
    upper = c.accel_max
    if remaining_max > 0:
        upper = min(upper, headroom / (2 * remaining_max))
    if lower > upper:
        return None

    earliest, latest = 0.0, horizon
    if times:
        earliest = _crossing(c.gap_min, c.speed_max, d, lower, lag_min, horizon)[0]
        latest = _crossing(c.gap_max, c.speed_min, -d, lower, lag_max, horizon)[1]
    energy_min = min(intercept + lower * slope for intercept, slope in energy_lines)
    return Certificate(lower, upper, earliest, latest, sqrt(max(0.0, energy_min)),
                       sqrt(energy_max + 2 * lower * remaining_max))
