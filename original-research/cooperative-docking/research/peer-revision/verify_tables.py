"""Recount reported revision outcomes from the per-run evidence."""
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REV = ROOT / "outputs/peer-revision"


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    batch = REV / "sensitivity-v1"
    summary = read(batch / "summary.json")
    assert sha(batch / "plan.json") == summary["plan_sha256"]
    assert sha(batch / "runs.json") == summary["runs_sha256"]
    rows = read(batch / "runs.json")
    assert len(rows) == 78 and not any(r["failures"] for r in rows)
    assert sum(r["outcome"] == "capture" for r in rows) == 64
    for cohort, caps, expected in (
        ("s_", "1p10_0p85", (0, 4)),
        ("s_", "1p10_1p90", (6, 6)),
        ("s_", "1p40_0p85", (6, 6)),
        ("s_", "1p40_1p90", (6, 6)),
        ("b_", "", (10, 10)),
    ):
        for policy, captures in zip(("fixed_adaptive", "corrected_adaptive"), expected):
            group = [r for r in rows if r["policy"] == policy and r["refinement"] == 1
                     and r["case"].startswith(cohort) and caps in r["case"]]
            count = sum(r["outcome"] == "capture" for r in group)
            assert count == captures
            print(cohort, caps, policy, count, len(group),
                  sum(r["requests"] - 1 for r in group),
                  round(sum(r["receiver_occupied_s"] for r in group), 4))
    audit = read(batch / "contact-audit.json")
    assert audit["captures"] == 64 and audit["max_abs_reference_error_s"] < 0.0001
    assert all(r["matches"] for r in summary["convergence"])
    assert read(REV / "model-checks.json")["scalar"]["violations"] == 0
    ipc = read(REV / "ipc-repeat/ipc_results.json")
    assert len(ipc) == 12 and not any(r["score_failures"] for r in ipc)
    print("Verified run counts, refinement, contact audit and IPC scores.")


if __name__ == "__main__":
    main()
