"""Run the paper's reference encounter with its pinned local OEL source."""

import argparse
import json
import os
import subprocess
import sys
from hashlib import sha256
from pathlib import Path
from tempfile import mkdtemp

import yaml

ROOT = Path(__file__).resolve().parent
SCENARIO = ROOT / "configs" / "cooperative_docking_paper_replay.yaml"
ORIGINAL_PLAN_SHA256 = "9ae81cd9219230cf98e74f9851c717a278353f084ccbb906917222382086aa25"
PLAN_ENV = "OEL_DOCKING_ARCHIVED_PLAN"
STUDY = "research/rpo-contact-leases-20261116"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def verify_pinned_source(engine):
    location = os.environ.get(PLAN_ENV)
    if not location:
        raise ValueError(f"Set {PLAN_ENV} to the original archived holdout plan for a trajectory replay.")
    original_plan = Path(location)
    if digest(original_plan) != ORIGINAL_PLAN_SHA256:
        raise ValueError("Original holdout plan digest mismatch.")
    plan = json.loads(original_plan.read_text())
    expected = plan["source_identity"]
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=engine, text=True).strip()
    if head != expected["head"]:
        raise RuntimeError(f"Expected OEL commit {expected['head']}, got {head}.")
    different = []
    for name, wanted in expected["engine_files"].items():
        file = engine / name
        got = digest(file) if file.is_file() else None
        if got != wanted:
            different.append(name)
    if different:
        raise RuntimeError(f"The pinned OEL source changed: {', '.join(different[:8])}")
    for name in ("certificate.py", "evaluate.py", "lease_stacks.py"):
        if digest(ROOT / "fsw" / name) != expected["study_files"][name]:
            raise RuntimeError(f"The pinned flight software changed: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", type=Path, required=True)
    args = parser.parse_args()
    engine = args.engine.expanduser().resolve()
    verify_pinned_source(engine)
    replay_root = ROOT / "outputs" / "replays"
    replay_root.mkdir(parents=True, exist_ok=True)
    run = Path(mkdtemp(prefix="reference-", dir=replay_root))
    evidence = run / "evidence"
    config = yaml.safe_load(SCENARIO.read_text())
    config["outputs"]["output_dir"] = str(evidence)
    selected = run / "config.yaml"
    selected.write_text(yaml.safe_dump(config, sort_keys=False))
    env = os.environ.copy()
    env.pop("OEL_ALLOW_EXTERNAL_CONFIG_PATHS", None)
    env.pop("OEL_CLI_ALLOW_EXTERNAL_CONFIG_PATHS", None)
    env.update(
        PYTHONPATH=os.pathsep.join((str(engine), str(ROOT / "fsw"))),
        OEL_WORKSPACE_ROOT=str(ROOT),
        OEL_OUTPUT_ROOT=str(ROOT / "outputs"),
        OEL_ENGINE_VERSION="1.0.0",
        OEL_INSTALLATION_DISPOSITION="developer",
        MPLCONFIGDIR=str(ROOT / ".oel" / "local-cache" / "matplotlib"),
    )
    for key in ("OEL_ENGINE_EDITION", "OEL_RELEASE_MANIFEST_SHA256", "OEL_INSTALLATION_TRANSACTION_ID"):
        env.pop(key, None)
    cmd = [sys.executable, "-B", str(engine / "run_simulation.py"), "--config", str(selected)]
    for name, flags in (("safe-inspection", ["--safe-validate"]),
                        ("validation", ["--validate-only"]), ("simulation", [])):
        log = run / f"{name}.log"
        with log.open("w") as stream:
            result = subprocess.run(cmd + flags, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{name} failed. See {log}")
    fresh = json.loads((evidence / "docking.json").read_text())
    case = ROOT / STUDY / "runs" / "evidence-w4" / "early_valid" / "output" / "docking.json"
    reference = json.loads(case.read_text())
    if fresh["events"] != reference["events"] or fresh["events"][0]["event"] != "first_contact":
        raise RuntimeError(f"Reference event comparison failed. Evidence is in {run}.")
    review = evidence / "review" / "run.sqlite"
    receipt = {"output": str(evidence), "review_sha256": digest(review),
               "reference_events_match": True, "first_contact_s": fresh["events"][0]["time_s"],
               "terminal_event": fresh["events"][-1]["event"]}
    (run / "replay-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
