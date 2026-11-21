"""Reviewer-only check of a shared 62.5 ms processing floor in the existing simulator."""
from pathlib import Path
import json
import sys
import shutil
from concurrent.futures import ThreadPoolExecutor
import run_comparison as rc
from run_study import score_telemetry

workspace = rc.ROOT
review = Path(__file__).resolve().parent
engine = workspace.parents[2] / "build/worktrees/v1.1.0"
rc.verify_pinned_source(engine)
(review / "fsw").mkdir()
for f in (workspace / "fsw").glob("*.py"):
    shutil.copy2(f, review / "fsw" / f.name)
shutil.copy2(workspace / "oel-workspace.yaml", review / "oel-workspace.yaml")
f = review / "fsw/comparator_stack.py"
src = f.read_text()
a = src.index('        if self.policy not in ("corrected_adaptive", "elapsed_adaptive"):')
b = src.index('    def options(', a)
src = src[:a] + '''        if self.policy != "floor_adaptive":
            return fixed
        elapsed_upper = (current_local - self.requests_sent[grant.number] + 1) / NS / (1 - DRIFT)
        age_upper = max(0.0, elapsed_upper - 0.0625)
        return max(fixed, grant.lease_ns / NS / (1 + DRIFT) - age_upper)

''' + src[b:]
f.write_text(src)
rc.ROOT = review
out = review / "outputs"
out.mkdir()
def run(case):
    folder = out / case.name
    folder.mkdir()
    cfg = rc.make_config(folder, case, "floor_adaptive")
    assert case.turnaround_s >= .0625
    rc.run_oel(engine, folder, cfg)
    r = rc.read_result(folder, cfg, case, "floor_adaptive")
    r["occupancy_s"] = score_telemetry(folder)
    print(case.name, r["outcome"], r["requests"], r["failures"], flush=True)
    return r
with ThreadPoolExecutor(max_workers=2) as pool:
    rows = list(pool.map(run, rc.CASES))
assert all(not r['failures'] for r in rows)
(out / "review-scores.json").write_text(json.dumps(rows,indent=2))
rc.verify_pinned_source(engine)
