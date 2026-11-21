"""Replay the fixed-time Lambert experiment using a pinned OEL source installation."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def engine_identity():
    import sim
    from importlib.metadata import version
    from sim.installation.workspace import load_workspace
    from sim.project_version import source_project_version
    root = Path(sim.__file__).resolve().parent.parent
    pin = json.loads((ROOT / "engine-revision.json").read_text())
    workspace = load_workspace(ROOT)
    actual_version = source_project_version(source_root=root) or version("orbital-engineering-lab")
    if workspace["engine"]["locked_version"] != pin["package_version"] or actual_version != pin["package_version"]:
        raise ValueError("Workspace and engine revision lock disagree.")
    found = source_digest(root)
    if found not in pin["accepted_runtime_sha256"].values():
        raise ValueError("The installed OEL Python source differs from the tested development and export builds.")
    return found


def source_digest(root):
    hashes = hashlib.sha256()
    for path in sorted((root / "sim").rglob("*.py")):
        if "tests" in path.relative_to(root / "sim").parts:
            continue
        hashes.update(path.relative_to(root).as_posix().encode() + b"\0")
        hashes.update(hashlib.sha256(path.read_bytes()).digest())
    return hashes.hexdigest()


def propagate(problem, branch, steps):
    from sim.dynamics.orbit.two_body import propagate_two_body_rk4
    r1, r2 = np.array(problem["r1_km"]), np.array(problem["r2_km"])
    axis = r1 / np.linalg.norm(r1)
    normal = np.cross(r1, r2)
    normal /= np.linalg.norm(normal)
    transverse = np.cross(normal, axis)
    state = np.r_[r1, branch.v1_km_s]
    dt = problem["time_of_flight_s"] / steps
    zero = np.zeros(3)
    angles = [0.0]
    for _ in range(steps):
        state = propagate_two_body_rk4(state, dt, problem["mu_km3_s2"], zero)
        angles.append(math.atan2(state[:3] @ transverse, state[:3] @ axis))
    return {"steps":steps, "position_error_km":float(np.linalg.norm(state[:3]-r2)),
            "velocity_error_km_s":float(np.linalg.norm(state[3:]-branch.v2_km_s)),
            "signed_turns":float(np.unwrap(angles)[-1]/math.tau)}


def compute(problem):
    from sim.dynamics.orbit import solve_multi_revolution_lambert
    identity = engine_identity()
    result = solve_multi_revolution_lambert(problem["r1_km"], problem["r2_km"], problem["time_of_flight_s"],
        revolutions=problem["revolutions"], mu_km3_s2=problem["mu_km3_s2"], short_way=problem["short_way"],
        tolerance_s=problem["solver_tolerance_s"])
    rows = [{"branch":b.branch,"v1_km_s":b.v1_km_s.tolist(),"v2_km_s":b.v2_km_s.tolist(),
             "time_residual_s":b.residual_s,"checks":[propagate(problem,b,n) for n in problem["propagation_steps"]]}
            for b in result.branches]
    return {"engine_sha256":identity,"minimum_time_s":result.minimum_time_s,"branches":rows}


def compare(problem, result):
    from itertools import permutations
    from reference import reference
    reference_result = reference(problem)
    branches = result["branches"]
    expected = reference_result["branches"]
    if len(branches) != len(expected):
        return {"passed":False, "reason":"solution count differs", "reference":reference_result}
    def discrepancy(order):
        return max(np.linalg.norm(np.subtract(b[k],e[k])) for b,e in zip(branches,order)
                   for k in ("v1_km_s","v2_km_s"))
    velocity_error = float(min(discrepancy(order) for order in permutations(expected)))
    minimum_error = abs(result["minimum_time_s"]-reference_result["minimum_time_s"])
    limits = problem["limits"]
    r1, r2 = np.array(problem["r1_km"]), np.array(problem["r2_km"])
    angle = math.atan2(np.linalg.norm(np.cross(r1,r2)), r1@r2)/math.tau
    winding = problem["revolutions"]+angle if problem["short_way"] else -(problem["revolutions"]+1-angle)
    checks = [row for b in branches for row in b["checks"]]
    passed = (velocity_error < 1e-8 and minimum_error < 0.001 and
              all(row["position_error_km"] < limits["position_km"] and
                  row["velocity_error_km_s"] < limits["velocity_km_s"] and
                  abs(row["signed_turns"]-winding) < limits["winding_turns"] for row in checks))
    return {"passed":bool(passed),"reference":reference_result,"max_reference_velocity_error_km_s":velocity_error,
            "minimum_time_difference_s":minimum_error}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    problem = json.loads((ROOT / "problem.json").read_text())
    result = compute(problem)
    result["comparison"] = compare(problem,result)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(result,stream,allow_nan=False,indent=2)
        stream.write("\n")
    return 0 if result["comparison"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
