"""Verify the equal-observation example, fixed-minimum control and IPC margins."""
import csv
import json
from math import fsum, isclose
from pathlib import Path
import sys
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[2]
FLOOR = ROOT / "outputs/referee-pass-20261118"
STUDY = ROOT / "research/rpo-contact-leases-20261116"
sys.path.insert(0, str(ROOT / "fsw"))
from certificate import certify, displacement, guaranteed_remaining, stopping_distance
from evaluate import DRIFT, HANDOFF_S, MAX_REVERSE_S, UPLINK_S, approach


def grid_minimum_admits(row):
    """Use the physical processing minimum in the analytical grid."""
    before = float(row["forward_s"])
    reverse = float(row["reverse_s"])
    assert before - UPLINK_S >= .03125
    elapsed_upper = (before + reverse) * (1 + DRIFT) / (1 - DRIFT)
    age_upper = min(MAX_REVERSE_S, max(0., elapsed_upper - .03125))
    remaining = float(row["lease_s"]) / (1 + DRIFT) - age_upper
    state = approach(float(row["gap_m"]), float(row["speed_m_s"]))
    return certify(state, remaining - HANDOFF_S, times=False) is not None


def example():
    state = approach(.034, .034)
    hold = .9
    with (STUDY / "evidence/decision_grid.csv").open() as stream:
        grid = list(csv.DictReader(stream))
        assert sum(map(grid_minimum_admits, grid)) == 518
        rows = [r for r in grid if
                (float(r["gap_m"]), float(r["speed_m_s"]), float(r["lease_s"])) ==
                (.034, .034, hold) and float(r["forward_s"]) + float(r["reverse_s"]) == .34375]
    selected = {float(r["reverse_s"]): r for r in rows}
    assert state.gap_min > .002 + stopping_distance(state.speed_max, state.lag_max, .08, state.disturbance)
    for reverse in (.03125, .28125):
        row = selected[reverse]
        before = float(row["forward_s"])
        elapsed = (before + reverse) * (1 + DRIFT)
        processing = (before - UPLINK_S) * (1 + DRIFT)
        remaining = guaranteed_remaining(hold, elapsed, DRIFT, DRIFT, processing)
        cert = certify(state, remaining - HANDOFF_S)
        assert (cert is not None) == (row["turnaround_corrected_admits"] == "True")
        if reverse == .03125:
            assert cert.acceleration < .042 < cert.acceleration_max
            assert isclose(cert.acceleration, float(row["turnaround_corrected_accel_m_s2"]))
        else:
            actual_horizon = hold / (1 + DRIFT) - reverse - HANDOFF_S
            closure = displacement(state.speed_min, -state.disturbance, state.accel_max,
                                   state.lag_max, actual_horizon)
            assert cert is None and round((state.gap_max - closure)*1000, 2) == 10.34
            print("Adverse realization: horizon, maximum closure, required clearance",
                  actual_horizon, closure, state.gap_max)
        print("Elapsed, receiver turnaround, remaining (s)", elapsed, processing, remaining)


def cohorts():
    corrected = {r["case"]: r for r in json.loads(
        (ROOT / "outputs/implementation-review/matched-v2/runs.json").read_text())
        if r["policy"] == "corrected_adaptive" and not r["case"].startswith("onset_")}
    groups = {"matched": [], "sensitivity": [], "broader": []}
    seen = set()
    for group in ("matched", "other-cohorts"):
        member = f"{group}/outputs/review-scores.json"
        path = FLOOR / member
        if path.is_file():
            rows = json.loads(path.read_bytes())
        else:
            with ZipFile(FLOOR / "reviewer-runs.zip") as full:
                rows = json.loads(full.read(member))
        for r in rows:
            case = r["case"]
            assert case not in seen and not r["failures"] and r["policy"] == "floor_adaptive"
            assert r["outcome"] == corrected[case]["outcome"]
            seen.add(case)
            cohort = "sensitivity" if case.startswith("s_") else "broader" if case.startswith("b_") else "matched"
            groups[cohort].append(r)
    assert seen == set(corrected) and len(seen) == 57
    published = json.loads((FLOOR / "comparison.json").read_text())
    for cohort, rows in groups.items():
        summary = published[cohort]
        assert (len(rows), sum(r["outcome"] == "capture" for r in rows), sum(r["requests"]-1 for r in rows)) == (
            summary["profiles"], summary["captures"], summary["renewals"])
        assert isclose(fsum(r["occupancy_s"] for r in rows), summary["occupied_s"], abs_tol=1e-12)
        assert sum(r["first_contact_s"] != corrected[r["case"]]["first_contact_s"] for r in rows) == summary["different_contact_times"]
        print(cohort, summary)


def ipc():
    for path, accepted in [(STUDY / "runs/mac-holdout-v1/ipc_results.json", 2),
                           (ROOT / "outputs/peer-revision/ipc-repeat/ipc_results.json", 6)]:
        rows = [r for r in json.loads(path.read_text()) if r["chaser"]["decision"] == "commit"]
        assert len(rows) == accepted
        errors = [abs(r["margin_to_expiry_after_handoff_s"] -
                      (r["receiver"]["rx_ns"] - r["chaser"]["tx_ns"]) / 1e9) for r in rows]
        assert max(errors) < 3e-9
        print("Accepted IPC trials, minimum margin (s), largest outbound discrepancy (ns)",
              accepted, min(r["margin_to_expiry_after_handoff_s"] for r in rows), max(errors)*1e9)


if __name__ == "__main__":
    example()
    cohorts()
    ipc()
