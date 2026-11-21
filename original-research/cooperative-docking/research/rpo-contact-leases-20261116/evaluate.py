"""Deterministic analytical workload and independent piecewise-motion check."""

import csv
import itertools
import json
import sys
from collections import Counter
from math import sqrt
from pathlib import Path
from random import Random

from certificate import Approach, certify, guaranteed_remaining

DRIFT = 0.001
HANDOFF_S = 0.04
MAX_REVERSE_S = 0.28125
MAX_ROUND_TRIP_S = 0.59375
UPLINK_S = 0.03125


def approach(gap, speed, lag_min=0.03125, lag_max=0.0625):
    return Approach(gap - 0.0006, gap + 0.0006, speed - 0.0012, speed + 0.0012,
                    lag_min, lag_max, 0.008, 0.05, 0.020, 0.070)


def assess(gap, speed, duration, forward, reverse):
    state = approach(gap, speed)
    shortest_hold = duration / (1 + DRIFT)
    # Use fastest requesting clock and fastest receiver for this benchmark.
    elapsed = (1 + DRIFT) * (forward + reverse)
    turnaround = (1 + DRIFT) * (forward - UPLINK_S)
    measured = guaranteed_remaining(duration, elapsed, DRIFT, DRIFT)
    corrected = guaranteed_remaining(duration, elapsed, DRIFT, DRIFT, turnaround)
    remaining = {
        "fixed_round_trip": shortest_hold - MAX_ROUND_TRIP_S,
        "fixed_reverse": shortest_hold - MAX_REVERSE_S,
        "measured_round_trip": measured,
        "turnaround_corrected": corrected,
        "combined": max(shortest_hold - MAX_REVERSE_S, corrected),
        "oracle": shortest_hold - reverse,
    }
    result = {name: certify(state, available - HANDOFF_S, times=False)
              for name, available in remaining.items()}
    return {
        "gap_m": gap, "speed_m_s": speed, "lease_s": duration, "forward_s": forward, "reverse_s": reverse,
        **{f"{name}_admits": solution is not None for name, solution in result.items()},
        **{f"{name}_accel_m_s2": solution.acceleration if solution else ""
           for name, solution in result.items()},
    }


def first_contact(gap, speed, lag, command, switch_times, disturbances):
    """Exact event solution for a specified piecewise constant scalar profile."""
    boundaries = sorted(set(switch_times) | {lag})
    elapsed = 0.0
    for stop in boundaries:
        active = command if elapsed >= lag else 0.0
        disturbance_index = sum(elapsed >= switch for switch in switch_times[:-1])
        acceleration = active + disturbances[disturbance_index]
        dt = stop - elapsed
        closure = speed * dt + acceleration * dt**2 / 2
        if closure >= gap:
            final_speed = sqrt(speed**2 + 2 * acceleration * gap)
            return elapsed + 2 * gap / (speed + final_speed), final_speed
        gap -= closure
        speed += acceleration * dt
        elapsed = stop
    return None


def independent_check(cases):
    rng = Random(20261116)
    count = 0
    max_time_violation = 0.0
    max_speed_violation = 0.0
    for row in cases:
        if not row["turnaround_corrected_admits"]:
            continue
        c = approach(row["gap_m"], row["speed_m_s"])
        local_elapsed = (1 + DRIFT) * (row["forward_s"] + row["reverse_s"])
        turnaround = (1 + DRIFT) * (row["forward_s"] - UPLINK_S)
        horizon = guaranteed_remaining(row["lease_s"], local_elapsed, DRIFT, DRIFT, turnaround) - HANDOFF_S
        cert = certify(c, horizon)
        corners = itertools.product((c.gap_min, c.gap_max), (c.speed_min, c.speed_max),
                                    (c.lag_min, c.lag_max), ((-1, -1), (-1, 1), (1, -1), (1, 1)))
        profiles = [(g, v, lag, [0.4 * horizon, horizon], [s * c.disturbance for s in signs])
                    for g, v, lag, signs in corners]
        for _ in range(8):
            switches = sorted(rng.uniform(0, horizon) for _ in range(4)) + [horizon]
            profiles.append((rng.uniform(c.gap_min, c.gap_max), rng.uniform(c.speed_min, c.speed_max),
                             rng.uniform(c.lag_min, c.lag_max), switches,
                             [rng.uniform(-c.disturbance, c.disturbance) for _ in switches]))
        for g, v, lag, switches, w in profiles:
            event = first_contact(g, v, lag, cert.acceleration, switches, w)
            if event is None:
                # At the attained endpoint, tolerate only floating-point roundoff.
                rounded_horizon = switches[-1] + 1e-12
                event = first_contact(g, v, lag, cert.acceleration, switches[:-1] + [rounded_horizon], w)
            if event is None:
                raise AssertionError((row, cert, g, v, lag, w))
            time, speed = event
            max_time_violation = max(max_time_violation, time - cert.latest, cert.earliest - time)
            max_speed_violation = max(max_speed_violation, speed - cert.contact_speed_max,
                                      cert.contact_speed_min - speed)
            count += 1
    assert max_time_violation < 2e-12 and max_speed_violation < 2e-12
    return {"piecewise_profiles_checked": count, "max_time_violation_s": max_time_violation,
            "max_speed_violation_m_s": max_speed_violation}


def clock_check():
    rng = Random(29031)
    smallest_slack = float("inf")
    for _ in range(20_000):
        lease = rng.uniform(0.1, 3.0)
        upstream = rng.uniform(0, 0.5)  # Physical request delay and receiver processing.
        processing = rng.uniform(0, upstream)
        downstream = rng.uniform(0, 4.0)
        c_rate = 1 + rng.uniform(-DRIFT, DRIFT)
        r_rate = 1 + rng.uniform(-DRIFT, DRIFT)
        c_epoch = rng.uniform(-10_000, 10_000)
        sent = c_epoch
        observed = c_epoch + c_rate * (upstream + downstream)
        bound = guaranteed_remaining(lease, observed - sent, DRIFT, DRIFT, processing * r_rate)
        observed_actual = lease / r_rate - downstream
        smallest_slack = min(smallest_slack, observed_actual - bound)
    # All unknown forward time can vanish with the clock rate extrema.
    hold = 1.2
    elapsed = 0.35
    tight = guaranteed_remaining(hold, elapsed, DRIFT, DRIFT)
    assert abs(tight - (hold / (1 + DRIFT) - elapsed / (1 - DRIFT))) < 1e-15
    assert smallest_slack >= -1e-10
    return {"clock_rate_and_latency_trials": 20_000, "smallest_remaining_time_slack_s": smallest_slack}


def main():
    output = Path(sys.argv[1])
    output.mkdir(parents=True, exist_ok=True)
    grid = [
        assess(*point) for point in itertools.product(
            (0.024, 0.034, 0.046, 0.060, 0.075, 0.090), (0.028, 0.034, 0.040, 0.046),
            (0.9, 1.2, 1.5, 1.8, 2.1), (0.0625, 0.1875, 0.3125), (0.03125, 0.15625, 0.28125))
    ]
    totals = Counter()
    for row in grid:
        for name, value in row.items():
            if name.endswith("_admits") and value:
                totals[name] += 1
        if row["fixed_round_trip_admits"]:
            assert row["oracle_admits"]
        if row["fixed_reverse_admits"]:
            assert row["combined_admits"]
        if row["combined_admits"]:
            assert row["oracle_admits"]
        if row["turnaround_corrected_admits"]:
            assert row["combined_admits"]
        if row["measured_round_trip_admits"]:
            assert row["turnaround_corrected_admits"]
    summary = {"grid_points": len(grid), "admissions": dict(totals), **independent_check(grid), **clock_check()}
    with (output / "decision_grid.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(grid[0]))
        writer.writeheader()
        writer.writerows(grid)
    (output / "analytical_verification.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
