"""Bracketing test with real process scheduling and emulated rate offsets."""

import json
import multiprocessing as mp
import sys
from pathlib import Path
from random import Random

from evaluate import DRIFT, HANDOFF_S, first_contact
from run_ipc import NS, REQUEST_STATE, chaser, receiver


def threshold():
    # Independently derive the fast feasible command for this fixture. The
    # lower speed requirement already holds for u=0 across the initial box.
    c = REQUEST_STATE
    assert c.speed_min**2 - 2 * c.disturbance * c.gap_max > c.contact_speed_min**2
    remaining = c.gap_max - (c.speed_min * c.lag_min - 0.5 * c.disturbance * c.lag_min**2)
    cap = min(c.accel_max, (c.contact_speed_max**2 - c.speed_max**2 - 2 * c.disturbance * c.gap_max)
              / (2 * remaining))
    contact, _ = first_contact(c.gap_max, c.speed_min, c.lag_max, cap, [2.0], [-c.disturbance])
    return contact + HANDOFF_S


def trial_specs(seed):
    rng = Random(seed)
    edge = threshold()
    trials = []
    for index, delta in enumerate((0.006, 0.012, 0.015, 0.019, 0.023, 0.029), start=1):
        lease = rng.uniform(0.78, 0.89)
        center_delay = lease / (1 + DRIFT) - 0.020 - edge
        for side, sign in (("inside", -1), ("outside", 1)):
            trials.append({"name": f"ipc_{index:02d}_{side}", "side": side, "lease_s": lease,
                           "reverse_s": center_delay + sign * delta, "required_remaining_s": edge,
                           "mode": "reordered" if (index + sign) % 3 == 0 else "normal"})
    rng.shuffle(trials)
    return trials


def trial(context, spec):
    end1, end2 = context.Pipe()
    result = context.Queue()
    owner = context.Process(target=receiver, args=(end1, result, spec["lease_s"], spec["reverse_s"], spec["mode"]))
    client = context.Process(target=chaser, args=(end2, result, 1.6))
    owner.start()
    client.start()
    end1.close()
    end2.close()
    try:
        rows = dict(result.get(timeout=5) for _ in range(2))
    finally:
        for participant in (owner, client):
            participant.join(timeout=2)
            if participant.is_alive():
                participant.terminate()
                participant.join()
        result.close()
    c, r = rows["chaser"], rows["receiver"]
    failures = []
    if owner.exitcode != 0 or client.exitcode != 0 or c["decision"] == "no_grant":
        failures.append("did not finish the exchange")
    if spec["mode"] == "reordered" and c["ignored"] != 1:
        failures.append("old reply was not ignored")
    # The fixture runs the chaser at its slow rate bound and the receiver at
    # its fast bound. The host's independent scoring timestamps recover
    # elapsed and processing durations directly, without the peer timestamps.
    physical_elapsed = (c["decision_ns"] - c["tx_ns"]) / NS
    physical_processing = (r["issued_ns"] - r["rx_ns"]) / NS
    hold = round(spec["lease_s"] * NS) / NS / (1 + DRIFT)
    reference = hold - (physical_elapsed - physical_processing)
    expected = "commit" if reference >= spec["required_remaining_s"] else "reject"
    if c["decision"] != expected:
        failures.append("measured-time classification")
    if abs(c.get("certified_remaining_s", 0) - reference) > 5e-9:
        failures.append("measured clock accounting")
    slack = None
    if c["decision"] == "commit":
        slack = (r["expiry_ns"] - c["decision_ns"]) / NS - c["latest_contact_delay_s"] - HANDOFF_S
        if slack < 0:
            failures.append("accepted deadline exceeds measured expiry")
    print(spec["name"], c["decision"], f"threshold residual {(reference-spec['required_remaining_s'])*1000:.2f} ms",
          flush=True)
    return {**spec, **rows, "measured_expected": expected, "reference_remaining_s": reference,
            "margin_to_expiry_after_handoff_s": slack, "score_failures": failures}


if __name__ == "__main__":
    plan_path = Path(sys.argv[1]) / "frozen_plan.json"
    out = plan_path.with_name("ipc_results.json")
    if out.exists():
        raise RuntimeError(f"Existing timing trial results at {out}.")
    plan = json.loads(plan_path.read_text())
    ctx = mp.get_context("spawn")
    results = [trial(ctx, case) for case in plan["ipc_trials"]]
    out.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
