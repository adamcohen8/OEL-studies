"""Compare archived OEL contact events against an independent axial root."""
import argparse
import json
from math import sin
from pathlib import Path

from scipy.optimize import brentq
from sim.review import ReviewWorkspace
from run_comparison import QUERY, NS, diagnostic_frames, digest, write_json

SAMPLES = """
SELECT time_s, gap_m, lateral_m, angle_rad, normal_speed_m_s, rate_rad_s
FROM docking_samples WHERE phase = 'approach' AND time_s >= ?
ORDER BY time_s
"""


def audit(folder, row, geometry):
    commit = row["commit_s"]
    review = ReviewWorkspace.open(folder / "result")
    sampled = review.query(SAMPLES, parameters=(commit,), max_rows=2000)
    diagnostics = review.query(QUERY, max_rows=5000, max_vm_steps=750_000)
    assert not sampled.truncated and not diagnostics.truncated
    states = sampled.rows
    initial = states[0]
    assert abs(initial["time_s"] - commit) < 1e-9
    telemetry = diagnostic_frames(diagnostics.rows, "chaser")
    decision = next(fields for t, fields in telemetry if t == round(commit * NS))
    cmd = decision["selected_acceleration_m_s2"]
    lag = geometry["lag_s"]
    applied = next((t / NS for t, f in telemetry if t >= round(commit * NS) and
                    f["commanded_axial_acceleration_m_s2"] > 0), None)
    start = commit + lag
    assert applied is None or abs(applied - start) < 1e-9
    center = initial["gap_m"]
    gap = center - 0.2 * sin(initial["angle_rad"])
    speed = -initial["normal_speed_m_s"]
    target = geometry["receiver_m_s2"]

    def residual(t):
        return gap - speed * t - target * t**2 / 2 - cmd * max(0, t - lag)**2 / 2

    delay = row["first_contact_s"] - commit
    ref = brentq(residual, 0, delay + 0.025, xtol=1e-13)
    margins = []
    for a, b in zip(states, states[1:]):
        t1, t2 = a["time_s"], b["time_s"]
        input_dv = cmd * (max(0, t2 - start) - max(0, t1 - start))
        observed_dv = a["normal_speed_m_s"] - b["normal_speed_m_s"]
        margins.append(abs(observed_dv - input_dv) / (t2 - t1))
    return {
        "folder": folder.name,
        "reference_error_s": delay - ref,
        "max_sampled_disturbance_m_s2": max(margins, default=0),
        "min_sampled_disk_slack_m": min(0.25 - 0.2 - s["lateral_m"] for s in states),
        "max_sampled_lead_speed_m_s": max(0.2 * s["rate_rad_s"] for s in states),
        "min_sampled_effective_gap_m": min(s["gap_m"] - 0.2 * sin(s["angle_rad"]) for s in states),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", type=Path)
    batch = parser.parse_args().batch
    plan = json.loads((batch / "plan.json").read_text())
    geoms = {case["name"]: case["geometry"] for case, _, _ in plan["trials"]}
    rows = json.loads((batch / "runs.json").read_text())
    summary = json.loads((batch / "summary.json").read_text())
    assert summary["plan_sha256"] == digest(batch / "plan.json")
    assert summary["runs_sha256"] == digest(batch / "runs.json")
    for row in rows:
        folder = batch / row["folder"]
        assert digest(folder / "config.yaml") == row["config_sha256"]
        assert digest(folder / "result/review/run.sqlite") == row["review_sha256"]
    results = [audit(batch / r["folder"], r, geoms[r["case"]])
               for r in rows if r["outcome"] == "capture"]
    summary = {
        "captures": len(results),
        "max_abs_reference_error_s": max(abs(r["reference_error_s"]) for r in results),
        "max_sampled_disturbance_m_s2": max(r["max_sampled_disturbance_m_s2"] for r in results),
        "min_sampled_disk_slack_m": min(r["min_sampled_disk_slack_m"] for r in results),
        "max_sampled_lead_speed_m_s": max(r["max_sampled_lead_speed_m_s"] for r in results),
        "rows": results,
    }
    write_json(batch / "contact-audit.json", summary)
    print({k: v for k, v in summary.items() if k != "rows"})


if __name__ == "__main__":
    main()
