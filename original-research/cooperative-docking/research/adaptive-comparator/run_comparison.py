"""Run matched OEL comparisons of the original and adaptive readiness policies."""

import argparse
import csv
import itertools
import json
import os
import statistics
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from math import cos, pi, sin
from pathlib import Path

import numpy as np
import yaml
from comparator_stack import NS, TOPIC
from evaluate import HANDOFF_S, MAX_REVERSE_S
from lease_stacks import TOPIC as WINDOW_TOPIC
from replay import ROOT, SCENARIO, digest, verify_pinned_source

from sim.review import ReviewWorkspace
from sim.utils.quaternion import quaternion_to_dcm_bn

POLICIES = ("fixed_one_shot", "fixed_adaptive", "fixed_renew_coast", "corrected_adaptive")
R_DFT = -740
C_DFT = 850
QUERY = f"""
SELECT object_id, topic, generated_time_ns, field_name, value_text, value_integer, value_real
FROM fsw_diagnostic_fields WHERE topic IN ('{TOPIC}', '{WINDOW_TOPIC}')
ORDER BY generated_time_ns, object_id, field_name
"""
GAP = "SELECT MIN(gap_m) AS minimum_gap_m FROM docking_samples"
SPEED = "SELECT normal_speed_m_s FROM docking_samples ORDER BY time_s DESC LIMIT 1"


@dataclass(frozen=True)
class Geometry:
    gap_m: float
    closing_m_s: float
    receiver_m_s2: float
    tilt_rad: float
    lag_s: float
    initial_lease_s: float
    gap_bias_m: float
    speed_bias_m_s: float


@dataclass(frozen=True)
class Case:
    name: str
    geometry: Geometry
    reverse_s: tuple
    second_cap_s: float
    turnaround_s: float = 0.0625
    drop: tuple = (False, False)


GEOMETRIES = {
    "level": Geometry(0.047, 0.040, 0.001, 0, 0.0625, 1.05, 0.00035, -0.00075),
    "moving": Geometry(0.052, 0.044, -0.0070, 0, 0.03125, 1.18, -0.00045, 0.0010),
    "tilted": Geometry(0.052, 0.036, 0.0068, 0.025, 0.0625, 1.24, 0.00035, 0.00065),
}
CASES = [
    Case(f"{label}_{speed}_{resource}", g, delays, cap)
    for (label, g), (speed, delays), (resource, cap) in itertools.product(
        GEOMETRIES.items(),
        (("fast", (0.03125, 0.0625)), ("middle", (0.15625, 0.09375)),
         ("slow", (0.28125, 0.21875))),
        (("renewable", 1.90), ("limited", 0.85)),
    )
] + [
    Case("first_reply_lost", GEOMETRIES["moving"], (0.09375, 0.125), 1.90, drop=(True, False)),
    Case("renewal_reply_lost", GEOMETRIES["level"], (0.09375, 0.125), 1.90, drop=(False, True)),
    Case("slow_processing", GEOMETRIES["tilted"], (0.03125, 0.09375), 1.90, turnaround_s=0.125),
]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def make_config(folder, case, policy):
    config = yaml.safe_load(SCENARIO.read_text())
    g = case.geometry
    config["scenario_name"] = f"{case.name}_{policy}"
    config["scenario_description"] = "Matched docking with speed management and causal renewal."
    config["outputs"]["output_dir"] = str(folder / "result")
    config["simulator"].update(dt_s=0.03125, duration_s=3.125)
    config["simulator"]["docking"]["max_step_s"] = 0.00390625
    profile = {
        "initial_request_ns": round(g.initial_lease_s * NS),
        "renewal_request_ns": round(1.9 * NS), "turn_ns": round(case.turnaround_s * NS),
        "caps_ns": [round(g.initial_lease_s * NS), round(case.second_cap_s * NS)],
        "reverse_ns": [round(delay * NS) for delay in case.reverse_s],
        "drop": list(case.drop),
    }
    for oid, obj in config["objects"].items():
        target = oid == "target"
        obj["flight_software"].update(module="comparator_stack", class_name="ComparatorStack")
        obj["flight_software"]["params"] = {
            "role": "receiver" if target else "chaser",
            "policy": policy, "session": config["scenario_name"], "link_profile": profile,
            "lag_ns": round(g.lag_s * NS), "epoch_ns": -209 * NS if target else 1301 * NS,
            "rate_ppm": R_DFT if target else C_DFT, "receiver_accel": g.receiver_m_s2,
            "measurement_bias": [g.gap_bias_m, g.speed_bias_m_s],
            "tilt_bound_rad": g.tilt_rad + 0.0015 if g.tilt_rad else 0,
            "resolve_wrench": bool(g.tilt_rad),
            "navigation_initialization": "ideal",
        }
    half = (pi - g.tilt_rad) / 2
    attitude = [cos(half), sin(half), 0, 0] if g.tilt_rad else [0, 1, 0, 0]
    offset = quaternion_to_dcm_bn(attitude).T @ np.array([0.0, 0.0, 0.5])
    position = (np.array([7000., 0., 0.]) + (np.array([0.006 if g.tilt_rad else 0.0, 0.,
                                                        0.5 + g.gap_m]) - offset) / 1000).tolist()
    config["objects"]["chaser"]["initial_state"].update(
        position_eci_km=position, velocity_eci_km_s=[0, 7.5, -g.closing_m_s / 1000],
        attitude_quat_bn=attitude,
    )
    file = folder / "config.yaml"
    file.write_text(yaml.safe_dump(config, sort_keys=False))
    return file


def run_oel(engine, folder, config):
    env = os.environ.copy()
    for key in ("OEL_ALLOW_EXTERNAL_CONFIG_PATHS", "OEL_CLI_ALLOW_EXTERNAL_CONFIG_PATHS",
                "OEL_ENGINE_EDITION", "OEL_RELEASE_MANIFEST_SHA256", "OEL_INSTALLATION_TRANSACTION_ID"):
        env.pop(key, None)
    env.update(PYTHONPATH=os.pathsep.join((str(engine), str(ROOT / "fsw"))),
               OEL_WORKSPACE_ROOT=str(ROOT), OEL_OUTPUT_ROOT=str(ROOT / "outputs"),
               OEL_ENGINE_VERSION="1.0.0", OEL_INSTALLATION_DISPOSITION="developer")
    cmd = [sys.executable, "-B", str(engine / "run_simulation.py"), "--config", str(config)]
    for name, flags in (("inspect", ["--safe-validate"]), ("validate", ["--validate-only"]), ("run", [])):
        log = folder / f"{name}.log"
        with log.open("w") as stream:
            process = subprocess.run(cmd + flags, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
        if process.returncode:
            raise RuntimeError(f"OEL {name} failed. See {log}")


def diagnostic_frames(rows, role):
    frames = {}
    for row in rows:
        if row["object_id"] != role:
            continue
        frame = frames.setdefault(row["generated_time_ns"], {})
        frame[row["field_name"]] = next(row[key] for key in ("value_text", "value_integer", "value_real")
                                        if row[key] is not None)
    return sorted(frames.items())


def sampled_retreat_clearance(review):
    samples = review.query("SELECT gap_m, angle_rad FROM docking_samples ORDER BY time_s",
                           max_rows=10000)
    if samples.truncated or not samples.rows:
        raise ValueError("Incomplete docking clearance samples.")
    return min(row["gap_m"] - 0.20 * sin(row["angle_rad"]) for row in samples.rows)


def read_result(folder, config, case, policy):
    output = folder / "result"
    review = ReviewWorkspace.open(output)
    diagnostic = review.query(QUERY, max_rows=5000, max_vm_steps=750_000)
    if diagnostic.truncated:
        raise RuntimeError("Incomplete diagnostic query.")
    chaser = diagnostic_frames(diagnostic.rows, "chaser")
    receiver = diagnostic_frames(diagnostic.rows, "target")
    issued, durations, received = {}, {}, {}
    for time_ns, fields in receiver:
        number = fields["issued_grant"]
        if number:
            issued.setdefault(number, time_ns)
            durations.setdefault(number, fields["issued_duration_ns"])
    for time_ns, fields in chaser:
        if fields["received_grant"]:
            received.setdefault(fields["received_grant"], time_ns)
    if any(number not in issued or not 0 < (t - issued[number]) / NS <= MAX_REVERSE_S
           for number, t in received.items()):
        raise ValueError("Simulated reply took longer than the shared fixed-delay guarantee.")
    minimum_gap = review.query(GAP).rows[0]["minimum_gap_m"]
    final_speed = review.query(SPEED).rows[0]["normal_speed_m_s"]
    truth = json.loads((output / "docking.json").read_text())
    contact = next((e for e in truth["events"] if e["event"] == "first_contact"), None)
    soft_capture = next((e for e in truth["events"] if e["event"] == "soft_capture"), None)
    commit = next(((t, fields) for t, fields in chaser if fields["mode"] == "commit"), None)
    abort = next(((t, fields) for t, fields in chaser if fields["mode"].startswith("abort")), None)
    failures = []
    outcome = "retreat" if abort and not contact else "unresolved"
    minimum_edge = None
    headroom = None
    chosen = None
    if contact and commit and soft_capture and contact["capture_armed"]:
        outcome = "capture"
        time_ns, f = commit
        chosen = f["chosen_grant"]
        physical_expiry = issued[chosen] / NS + durations[chosen] / NS / (1 + R_DFT / 1_000_000)
        earliest, latest = [time_ns / NS + f[key] for key in ("earliest_s", "latest_s")]
        headroom = physical_expiry - latest - HANDOFF_S
        speed = -contact["normal_speed_m_s"]
        if not earliest - 1e-8 <= contact["time_s"] <= latest + 1e-8 or headroom < -1e-8:
            failures.append("contact outside certified readiness window")
        if not f["contact_speed_min_m_s"] - 1e-8 <= speed <= f["contact_speed_max_m_s"] + 1e-8:
            failures.append("contact outside certified speed window")
        if any(e["event"] in {"latch_release", "rebound"} and e["time_s"] <= contact["time_s"] + HANDOFF_S
               for e in truth["events"]):
            failures.append("lost capture during handoff")
    elif contact:
        outcome = "unsafe_contact"
        failures.append("contact without certified soft capture")
    elif outcome == "retreat":
        minimum_edge = sampled_retreat_clearance(review)
        if minimum_edge < 0.002 - 1e-7 or final_speed <= 0:
            failures.append("retreat reserve or settled separation not observed")
    else:
        failures.append("no terminal result")
    actual = {
        "case": case.name, "policy": policy, "outcome": outcome, "failures": failures,
        "requests": max(f["requests"] for _, f in chaser), "received": len(received),
        "slowdown": any(f["entered_slowdown"] for _, f in chaser),
        "first_declined": any(f["first_declined_before_renewal"] for _, f in chaser),
        "commit_grant": chosen, "commit_s": commit[0] / NS if commit else None,
        "abort_s": abort[0] / NS if abort else None,
        "first_contact_s": contact["time_s"] if contact else None,
        "contact_speed_m_s": -contact["normal_speed_m_s"] if contact else None,
        "min_gap_m": minimum_gap, "min_sampled_edge_gap_m": minimum_edge, "deadline_margin_beyond_handoff_s": headroom,
        "granted_total_s": sum(durations.values()) / NS,
        "config_sha256": digest(config), "review_sha256": digest(output / "review" / "run.sqlite"),
    }
    return actual


def run_case(engine, output, trial):
    case, policy = trial
    folder = output / f"{case.name}_{policy}"
    folder.mkdir()
    try:
        config = make_config(folder, case, policy)
        run_oel(engine, folder, config)
        result = read_result(folder, config, case, policy)
    except (RuntimeError, ValueError) as exc:
        result = {"case": case.name, "policy": policy, "outcome": "incomplete", "failures": [str(exc)]}
    print(case.name, policy, result["outcome"], result.get("requests"), result.get("first_contact_s"),
          result["failures"], flush=True)
    return result


def summarize(cases, rows, trial_policies):
    if any(row["failures"] for row in rows):
        return {"scorable": False, "discrepancies": [row for row in rows if row["failures"]]}
    by_policy = {}
    for policy in trial_policies:
        group = [r for r in rows if r["policy"] == policy]
        captured = [r for r in group if r["outcome"] == "capture"]
        by_policy[policy] = {
            "captures": len(captured), "retreats": len(group) - len(captured),
            "first_grant_captures": sum(r["commit_grant"] == 1 for r in captured),
            "renewals_requested": sum(r["requests"] == 2 for r in group),
            "slowed": sum(r["slowdown"] for r in group),
            "total_granted_s": sum(r["granted_total_s"] for r in group),
            "retreat_minimum_gap_mm": min((r["min_gap_m"] * 1000 for r in group if r["outcome"] == "retreat"),
                                          default=None),
        }
    paired = []
    by_key = {(r["case"], r["policy"]): r for r in rows}
    for case in cases:
        fixed = by_key[case.name, "fixed_adaptive"]
        improved = by_key[case.name, "corrected_adaptive"]
        paired.append({
            "case": case.name, "fixed": fixed["outcome"], "corrected": improved["outcome"],
            "fixed_requests": fixed["requests"], "corrected_requests": improved["requests"],
            "fixed_contact_s": fixed["first_contact_s"], "corrected_contact_s": improved["first_contact_s"],
            "fixed_minus_corrected_s": (fixed["first_contact_s"] - improved["first_contact_s"]
                                        if fixed["outcome"] == improved["outcome"] == "capture" else None),
        })
    deltas = [r["fixed_minus_corrected_s"] for r in paired if r["fixed_minus_corrected_s"] is not None]
    avoided = [r for r in paired if r["fixed_requests"] > r["corrected_requests"]
               and r["fixed_minus_corrected_s"] is not None]
    fixed_union = []
    for case in cases:
        fixed_results = [by_key[case.name, policy] for policy in ("fixed_adaptive", "fixed_renew_coast")]
        matches = [row["first_contact_s"] for row in fixed_results if row["outcome"] == "capture"]
        fixed_union.append({"case": case.name, "best_fixed_contact_s": min(matches) if matches else None,
                            "corrected_contact_s": by_key[case.name, "corrected_adaptive"]["first_contact_s"]})
    strict_duels = [
        r["best_fixed_contact_s"] - r["corrected_contact_s"] for r in fixed_union
        if r["best_fixed_contact_s"] is not None and r["corrected_contact_s"] is not None
    ]
    return {
        "scorable": True, "profiles": len(cases), "runs": len(rows), "policies": by_policy,
        "corrected_only_capture": [r["case"] for r in paired if r["corrected"] == "capture" and r["fixed"] != "capture"],
        "fixed_only_capture": [r["case"] for r in paired if r["fixed"] == "capture" and r["corrected"] != "capture"],
        "paired_capture_median_time_saved_s": statistics.median(deltas) if deltas else None,
        "paired_capture_time_saved_s": deltas,
        "both_capture_fixed_retries_avoided": avoided,
        "fixed_either_recovery_capture_count": sum(r["best_fixed_contact_s"] is not None for r in fixed_union),
        "corrected_only_vs_either_fixed": [r["case"] for r in fixed_union
                                          if r["corrected_contact_s"] is not None and r["best_fixed_contact_s"] is None],
        "either_fixed_only_vs_corrected": [r["case"] for r in fixed_union
                                           if r["corrected_contact_s"] is None and r["best_fixed_contact_s"] is not None],
        "time_saved_vs_fastest_fixed_s": strict_duels,
        "paired_results": paired,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, default=ROOT.parents[2] / "build" / "worktrees" / "v1.1.0")
    parser.add_argument("batch", help="New directory name under outputs/adaptive-comparator.")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if not args.batch.replace("-", "").replace("_", "").isalnum():
        raise ValueError("Batch name must be a single slug.")
    engine = args.engine.resolve()
    verify_pinned_source(engine)
    output = ROOT / "outputs" / "adaptive-comparator" / args.batch
    output.mkdir(parents=True)
    cases = [CASES[0], CASES[-3], CASES[-1]] if args.smoke else CASES
    trial_policies = POLICIES
    files = [Path(__file__).resolve(), ROOT / "fsw/comparator_stack.py", SCENARIO]
    frozen = {"cases": [asdict(case) for case in cases], "policies": trial_policies,
              "source": {str(path.relative_to(ROOT)): digest(path) for path in files}}
    plan = output / "plan.json"
    write_json(plan, frozen)
    plan_sha = digest(plan)
    (output / "plan.sha256").write_text(plan_sha + "\n")
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(lambda trial: run_case(engine, output, trial),
                             itertools.product(cases, trial_policies)))
    assert digest(plan) == plan_sha
    assert all(digest(ROOT / file) == expected for file, expected in frozen["source"].items())
    verify_pinned_source(engine)
    write_json(output / "runs.json", rows)
    write_csv(output / "runs.csv", rows)
    summary = summarize(cases, rows, trial_policies)
    write_json(output / "summary.json", summary)
    if summary["scorable"]:
        write_csv(output / "pairs.csv", summary["paired_results"])
    print(json.dumps({k: val for k, val in summary.items() if k != "paired_results"}, indent=2))
    return 0 if summary["scorable"] else 1


if __name__ == "__main__":
    sys.exit(main())
