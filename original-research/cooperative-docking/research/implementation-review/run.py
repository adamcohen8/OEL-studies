"""Frozen paired attribution experiment and off-tick actuator tests."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path

from run_comparison import (CASES, ROOT, SCENARIO, digest, make_config, read_result,
                            run_oel, verify_pinned_source, write_json)
from run_study import design, score_telemetry

POLICIES = ("fixed_adaptive", "elapsed_adaptive", "corrected_adaptive")


def trials():
    cases = CASES + [case for case, policy, refinement in design()
                     if policy == "fixed_adaptive" and refinement == 1]
    rows = [(case, policy) for case in cases for policy in POLICIES]
    level = next(c for c in CASES if c.name == "level_fast_limited")
    for lag_ns in (31250000, 31250001, 46875000, 62499999, 62500000):
        case = replace(level, name=f"onset_{lag_ns}",
                       geometry=replace(level.geometry, lag_s=lag_ns / 1e9))
        rows.extend((case, policy) for policy in ("corrected_adaptive", "fixed_one_shot"))
    return rows


def run(engine, output, case, policy):
    folder = output / f"{case.name}_{policy}"
    folder.mkdir()
    config = make_config(folder, case, policy)
    run_oel(engine, folder, config)
    row = read_result(folder, config, case, policy)
    row["receiver_occupied_s"] = score_telemetry(folder)
    write_json(folder / "score.json", row)
    print(case.name, policy, row["outcome"], row["failures"], flush=True)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch")
    args = parser.parse_args()
    if not args.batch.replace("-", "").isalnum():
        parser.error("use a single batch name")
    engine = ROOT.parents[2] / "build/worktrees/v1.1.0"
    verify_pinned_source(engine)
    output = ROOT / "outputs/implementation-review" / args.batch
    output.mkdir(parents=True)
    work = trials()
    sources = [Path(__file__), ROOT / "research/adaptive-comparator/run_comparison.py",
               ROOT / "research/peer-revision/run_study.py", ROOT / "replay.py", SCENARIO,
               *(ROOT / "fsw").glob("*.py")]
    def identity():
        return {str(p.relative_to(ROOT)): digest(p) for p in sources}
    source = identity()
    plan = output / "plan.json"
    write_json(plan, {"trials": [[asdict(c), p] for c, p in work], "source": source,
                      "onset_expected": {"corrected_adaptive": "capture", "fixed_one_shot": "retreat"},
                      "engine_record_sha256": digest(ROOT / ".oel/reproduction/engine-reproduction.json")})
    plan_hash = digest(plan)
    with ThreadPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(lambda pair: run(engine, output, *pair), work))
    verify_pinned_source(engine)
    assert source == identity() and plan_hash == digest(plan)
    write_json(output / "runs.json", rows)
    write_json(output / "identity.json", {"plan_sha256": plan_hash,
                                        "runs_sha256": digest(output / "runs.json")})
    assert not any(r["failures"] for r in rows)


if __name__ == "__main__":
    main()
