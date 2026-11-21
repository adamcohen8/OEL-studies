"""Check the abstract reservation contract and scalar dynamics independently."""
from itertools import product
import json
from pathlib import Path
from random import Random
from math import pi, sin, sqrt

from scipy.integrate import solve_ivp

from certificate import certify
from evaluate import HANDOFF_S, approach


def contract_check():
    """Enumerate traces of an abstract receiver with a retained grant ledger."""
    actions = ("new_a", "new_b", "duplicate_a", "abort", "reset", "loss", "tick")
    for sequence in product(actions, repeat=6):
        ledger = {}
        commitments = []
        operational = True
        for t, action in enumerate(sequence):
            busy = any(expiry > t for expiry in ledger.values())
            if action in ("abort", "reset") and not busy:
                operational = False
            if operational and action.startswith("new"):
                identity = action[-1]
                if identity not in ledger:
                    ledger[identity] = t + (3 if identity == "a" else 2)
                    commitments.append((t, ledger[identity]))
            if action == "duplicate_a" and "a" in ledger and ledger["a"] > t:
                assert ledger["a"] == next(e for start, e in commitments if e - start == 3)
            assert not any(start <= t < expiry for start, expiry in commitments) or operational
            assert busy or not any(start < t < expiry for start, expiry in commitments)
    return len(actions) ** 6


def scalar_check():
    rng = Random(20261118)
    accepted = 0
    violations = 0
    min_time_slack = min_speed_slack = float("inf")
    for _ in range(512):
        dstate = approach(rng.uniform(0.024, 0.09), rng.uniform(0.028, 0.046))
        horizon = rng.uniform(0.70, 2.30)
        solution = certify(dstate, horizon)
        if solution is None:
            continue
        accepted += 1
        initial = (rng.uniform(dstate.gap_min, dstate.gap_max),
                   rng.uniform(dstate.speed_min, dstate.speed_max))
        lag = rng.uniform(dstate.lag_min, dstate.lag_max)
        frequency = rng.uniform(0.2, 12)
        phase = rng.uniform(0, 2 * pi)

        def integrate(start, stop, state, command):
            def dynamics(t, y):
                return (-y[1], command + dstate.disturbance * sin(2 * pi * frequency * t + phase))

            def contact(t, y):
                return y[0]

            contact.terminal = True
            contact.direction = -1
            return solve_ivp(dynamics, (start, stop), state, events=contact,
                             rtol=2e-10, atol=1e-12, max_step=0.015625)

        result = integrate(0, min(lag, horizon), initial, 0)
        if not result.t_events[0].size and lag < horizon:
            result = integrate(lag, horizon, result.y[:, -1], solution.acceleration)
        if not result.success or not result.t_events[0].size:
            violations += 1
            continue
        time = result.t_events[0][0]
        speed = result.y_events[0][0, 1]
        slack_t = min(time - solution.earliest, solution.latest - time)
        slack_v = min(speed - solution.contact_speed_min, solution.contact_speed_max - speed)
        min_time_slack = min(min_time_slack, slack_t)
        min_speed_slack = min(min_speed_slack, slack_v)
        violations += int(slack_t < -1e-8 or slack_v < -1e-8)
    assert accepted > 0
    return {"sampled_states": 512, "admitted_and_integrated": accepted, "violations": violations,
            "min_time_slack_s": min_time_slack, "min_speed_slack_m_s": min_speed_slack}


def nominal_reference():
    """Count nominally feasible refusal points with exact time and state."""
    import csv
    grid = Path(__file__).parents[1] / "rpo-contact-leases-20261116/evidence/decision_grid.csv"
    with grid.open() as stream:
        rows = list(csv.DictReader(stream))
    exact_refusals = [row for row in rows if row["oracle_admits"] == "False"]
    witnesses = 0
    for row in exact_refusals:
        gap, speed, duration, age = [float(row[name]) for name in
                                    ("gap_m", "speed_m_s", "lease_s", "reverse_s")]
        horizon = duration / 1.001 - age - HANDOFF_S
        lag = 0.046875
        remaining = gap - speed * lag
        for i in range(201):
            command = 0.05 * i / 200
            impact_speed = sqrt(speed**2 + 2 * command * remaining)
            impact_time = lag + 2 * remaining / (speed + impact_speed)
            if impact_time <= horizon and 0.02 <= impact_speed <= 0.07:
                witnesses += 1
                break
    assert rows and exact_refusals
    return {"grid_points": len(rows), "exact_expiry_certificate_refusals": len(exact_refusals),
            "nominal_witnesses_among_refusals": witnesses}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    path = parser.parse_args().output
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    result = {"abstract_receiver_traces": contract_check(), "scalar": scalar_check(),
              "nominal_reference": nominal_reference()}
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(result)
