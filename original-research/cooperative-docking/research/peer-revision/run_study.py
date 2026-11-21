"""Matched reservation sensitivity, unfiltered cases and OEL step refinement."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace
from itertools import product
import json
from pathlib import Path
from random import Random
import subprocess
import sys

from sim.review import ReviewWorkspace
from run_comparison import (Case, Geometry, GEOMETRIES, ROOT, SCENARIO, QUERY, NS,
                            R_DFT, diagnostic_frames, digest, make_config, read_result,
                            run_oel, verify_pinned_source, write_csv, write_json)
import yaml

POLICIES = ("fixed_adaptive", "corrected_adaptive")
SEED = 20261117
TOL_S = 0.0001


def design():
    sensitivity = [
        Case(f"s_{label}_{first:.2f}_{cap:.2f}_{link}".replace(".", "p"),
             replace(g, initial_lease_s=first), (0.041, delay), cap)
        for (label, g), first, cap, (link, delay) in product(
            GEOMETRIES.items(), (1.10, 1.40), (0.85, 1.90),
            (("early", 0.052), ("late", 0.201)))
    ]
    rng = Random(SEED)
    broader = [
        Case(f"b_{i:02}", Geometry(rng.uniform(0.041, 0.071), rng.uniform(0.029, 0.046),
                                   rng.uniform(-0.007, 0.007), rng.uniform(0, 0.032),
                                   rng.choice((0.03125, 0.0625)), rng.uniform(0.90, 1.65),
                                   rng.uniform(-0.0005, 0.0005), rng.uniform(-0.001, 0.001)),
             (rng.uniform(0.020, 0.240), rng.uniform(0.020, 0.240)),
             rng.uniform(0.80, 1.90), rng.choice((0.0625, 0.125)),
             (False, rng.random() < 0.2))
        for i in range(12)
    ]
    selected = [sensitivity[i] for i in (1, 9, 17)]
    return ([(c, p, 1) for c, p in product(sensitivity + broader, POLICIES)]
            + [(c, POLICIES[1], n) for c, n in product(selected, (2, 4))])


def score_telemetry(folder):
    review = ReviewWorkspace.open(folder / "result")
    diagnostics = review.query(QUERY, max_rows=5000, max_vm_steps=750_000)
    if diagnostics.truncated:
        raise ValueError("Incomplete diagnostics")
    grants = {}
    for time_ns, fields in diagnostic_frames(diagnostics.rows, "target"):
        if fields["issued_grant"]:
            grants.setdefault(fields["issued_grant"], (
                time_ns / NS, fields["issued_duration_ns"] / NS / (1 + R_DFT / 1e6)))
    end = occupied = 0.0
    for start, duration in sorted(grants.values()):
        occupied += max(0, start + duration - max(start, end))
        end = max(end, start + duration)
    return occupied


def execute(engine, output, trial):
    case, policy, refinement = trial
    folder = output / f"{case.name}_{policy}_r{refinement}"
    folder.mkdir()
    try:
        config = make_config(folder, case, policy)
        if refinement > 1:
            raw = yaml.safe_load(config.read_text())
            raw["simulator"]["dt_s"] /= refinement
            raw["simulator"]["docking"]["max_step_s"] /= refinement
            config.write_text(yaml.safe_dump(raw, sort_keys=False))
        run_oel(engine, folder, config)
        row = read_result(folder, config, case, policy)
        row["receiver_occupied_s"] = score_telemetry(folder)
    except (RuntimeError, ValueError) as exc:
        row = {"case": case.name, "policy": policy, "outcome": "incomplete", "failures": [str(exc)]}
    row.update(refinement=refinement, folder=folder.name)
    write_json(folder / "score.json", row)
    print(case.name, policy, refinement, row["outcome"], row["failures"], flush=True)
    return row


def group(rows):
    return {policy: {
        "captures": sum(r["outcome"] == "capture" for r in rows if r["policy"] == policy),
        "retreats": sum(r["outcome"] == "retreat" for r in rows if r["policy"] == policy),
        "renewal_requests": sum(r["requests"] - 1 for r in rows if r["policy"] == policy),
        "receiver_occupied_s": sum(r["receiver_occupied_s"] for r in rows if r["policy"] == policy),
    } for policy in POLICIES}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch")
    parser.add_argument("--engine", type=Path, default=ROOT.parents[2] / "build/worktrees/v1.1.0")
    args = parser.parse_args()
    if not args.batch.replace("-", "").isalnum():
        parser.error("batch must be a slug")
    engine = args.engine.resolve()
    verify_pinned_source(engine)
    work = design()
    output = ROOT / "outputs/peer-revision" / args.batch
    output.mkdir(parents=True)
    files = [Path(__file__), Path(__file__).parents[1] / "adaptive-comparator/run_comparison.py",
             SCENARIO, ROOT / "fsw/comparator_stack.py", ROOT / "fsw/certificate.py",
             ROOT / "fsw/lease_stacks.py", ROOT / "fsw/evaluate.py"]
    source = {str(p.resolve().relative_to(ROOT)): digest(p) for p in files}
    plan = output / "plan.json"
    write_json(plan, {"seed": SEED, "convergence_tolerance_s": TOL_S,
                      "trials": [[asdict(c), p, r] for c, p, r in work], "source": source,
                      "engine_head": subprocess.check_output(["git", "-C", engine, "rev-parse", "HEAD"],
                                                              text=True).strip(),
                      "python": sys.version})
    plan_hash = digest(plan)
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(lambda trial: execute(engine, output, trial), work))
    assert digest(plan) == plan_hash
    assert all(digest(ROOT / name) == sha for name, sha in source.items())
    verify_pinned_source(engine)
    write_json(output / "runs.json", rows)
    write_csv(output / "runs.csv", rows)
    base = [r for r in rows if r["refinement"] == 1]
    bad = [r for r in rows if r["failures"]]
    summary = {"runs": len(rows), "failures": bad, "plan_sha256": plan_hash,
               "runs_sha256": digest(output / "runs.json")}
    if not bad:
        for prefix in ("s_", "b_"):
            subset = [r for r in base if r["case"].startswith(prefix)]
            summary[prefix] = group(subset)
        summary["strata"] = {
            f"{first:.2f}/{cap:.2f}": group([
                r for r in base if r["case"].startswith("s_") and
                f"_{first:.2f}_{cap:.2f}_".replace(".", "p") in r["case"]])
            for first, cap in product((1.10, 1.40), (0.85, 1.90))
        }
        summary["convergence"] = []
        for refined in rows[len(base):]:
            coarse = next(r for r in base if (r["case"], r["policy"]) ==
                          (refined["case"], refined["policy"]))
            delta = (abs(coarse["first_contact_s"] - refined["first_contact_s"])
                     if coarse["outcome"] == refined["outcome"] == "capture" else None)
            summary["convergence"].append({
                "case": refined["case"], "refinement": refined["refinement"], "delta_s": delta,
                "matches": (coarse["outcome"], coarse["commit_grant"]) ==
                           (refined["outcome"], refined["commit_grant"]) and (delta is None or delta <= TOL_S)
            })
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    return bool(bad)


if __name__ == "__main__":
    sys.exit(main())
