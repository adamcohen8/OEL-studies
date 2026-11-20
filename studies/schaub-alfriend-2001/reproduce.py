"""Propagate the two scenarios with the pinned OEL engine and measure phase drift."""
import argparse
import json
import math
from pathlib import Path
from statistics import mean
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
SAMPLES = 64
ORBITS = 45
PERIOD = 6021
REFERENCE_DEG = {
    "published": (0.00533648, -0.0000174839),
    "variation": (-0.00458144, -0.0000165002),
}
QUERY = """SELECT object_id, sample_index, time_s, ecc, argp_deg, true_anomaly_deg,
                  conversion_status FROM object_orbital_elements ORDER BY object_id, sample_index"""


def measure(rows, method):
    phase = {}
    for member in ("chief", "deputy"):
        track = [r for r in rows if r["object_id"] == f"{method}_{member}"]
        if len(track) != ORBITS * SAMPLES + 1:
            raise ValueError(f"Incomplete history for {method}_{member}")
        angles = []
        for index, row in enumerate(track):
            if (row["sample_index"] != index or row["time_s"] != index * PERIOD / SAMPLES
                    or row["conversion_status"] != "ok"):
                raise ValueError(f"Invalid orbital elements or sampling grid at row {index}")
            f = math.radians(row["true_anomaly_deg"])
            e = row["ecc"]
            E = math.atan2(math.sqrt(1 - e * e) * math.sin(f), e + math.cos(f))
            angles.append(row["argp_deg"] + math.degrees(E - e * math.sin(E)))
        phase[member] = angles
    relative = [(d - c + 180) % 360 - 180 for c, d in zip(phase["chief"], phase["deputy"])]
    averages = [mean(relative[i:i + SAMPLES]) for i in range(0, ORBITS * SAMPLES, SAMPLES)]
    center = (ORBITS - 1) / 2
    return ORBITS * sum((i - center) * y for i, y in enumerate(averages)) / sum(
        (i - center) ** 2 for i in range(ORBITS))


def replay(cli, case, folder):
    def invoke(*args):
        return subprocess.check_output([*cli, *args], cwd=ROOT, text=True)

    config = str(ROOT / "configs" / f"{case}.yaml")
    out = str(folder / case)
    for flag in ("--safe-validate", "--validate-only"):
        invoke("sim", "--config", config, flag)
    # The lifecycle owner also writes the accepted and terminal records for later inspection.
    receipt = folder / f"{case}.jsonl"
    with receipt.open("w") as log, subprocess.Popen(
        [*cli, "runs", "start", "--config", config, "--output-dir", out, "--jsonl"],
        stdout=subprocess.PIPE, cwd=ROOT, text=True,
    ) as owner:
        line = owner.stdout.readline()
        log.write(line)
        log.flush()
        if not line:
            raise RuntimeError("OEL did not return an accepted run handle")
        handle = json.loads(line)
        identity = [handle["run_id"], "--expect-manifest-ref", handle["manifest_ref"],
                    "--expect-config-sha256", handle["normalized_config_sha256"]]
        observation = json.loads(invoke("runs", "await", *identity, "--timeout", "600"))
        log.write(owner.stdout.read())
        if owner.wait() or observation["outcome"] not in {"terminal", "still_running"}:
            raise RuntimeError(f"Inspect the run receipt at {receipt}")
    state = json.loads(invoke("runs", "inspect", *identity))["state"]
    if state["state"] != "completed" or not state["review"]["available"]:
        raise RuntimeError(f"Incomplete OEL run at {out}")
    result = json.loads(invoke("review", out, "--json", "--query", QUERY,
                               "--max-rows", "11524", "--max-vm-steps", "3000000"))
    if result["truncated"] or result["row_count"] != 11524:
        raise ValueError("Expected 2,881 records for each of four spacecraft")
    direct, mapped = (measure(result["rows"], method) for method in ("direct", "mapped"))
    target_direct, target_mapped = REFERENCE_DEG[case]
    factor = abs(direct / mapped)
    checks = [math.isclose(direct, target_direct, rel_tol=.05),
              math.isclose(mapped, target_mapped, rel_tol=.20),
              math.isclose(factor, abs(target_direct / target_mapped), rel_tol=.25)]
    return {"direct_deg": direct, "mapped_deg": mapped, "factor": factor, "passed": all(checks)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oel", default="oel", help="OEL launcher on PATH or absolute path")
    cli = [parser.parse_args().oel, "--workspace", str(ROOT)]
    version = json.loads(subprocess.check_output([*cli, "version"], text=True))
    if version["version"] != "1.1.0":
        raise SystemExit("This workspace requires OEL 1.1.0")
    outputs = ROOT / "outputs"
    outputs.mkdir(exist_ok=True)
    folder = Path(tempfile.mkdtemp(prefix="replay-", dir=outputs))
    print(f"Saving OEL evidence in {folder}", flush=True)
    results = {}
    for case in REFERENCE_DEG:
        results[case] = replay(cli, case, folder)
        print(case, json.dumps(results[case]), flush=True)
    (folder / "comparison.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0 if all(r["passed"] for r in results.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
