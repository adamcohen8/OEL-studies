"""Rebuild the compact numerical record from the replay output directories."""
from pathlib import Path
import json
import numpy as np


def main():
    root = Path(__file__).resolve().parent
    outputs = root / "outputs"
    results = {}
    for name in ("production", "refine", "targeted", "long", "halfstep"):
        folder = outputs / (name + "-c16")
        summary = json.loads((folder / "summary.json").read_text())
        batch = {key:summary[key] for key in ("physics","settings","calibration","tracks","periodic_orbits")}
        if name in ("long", "halfstep"):
            batch["windows"] = []
            with np.load(folder / summary["tracks"][0]["data"]) as track:
                for i in range(0, len(track["log_growth"]), 96):
                    points = track["peri"][i + 1:i + 97]
                    cells = np.histogram2d(points[:,0], points[:,1], bins=(12,12),
                                           range=((-np.pi/2,np.pi/2),(-1,3)))[0]
                    batch["windows"].append(dict(first_orbit=i+1, last_orbit=i+len(points),
                        mean_growth_per_orbit=float(np.mean(track["log_growth"][i:i+96])),
                        min_spin_over_n=float(min(points[:,1])), max_spin_over_n=float(max(points[:,1])),
                        occupied_of_144_cells=int(np.count_nonzero(cells))))
        results[name] = batch
    m = json.loads((outputs / "multireturn.json").read_text())
    results["multiple_returns"] = m["trials"]
    supplemental = json.loads((outputs / "paper_targets.json").read_text())
    results["source_guided_followup"] = {k:v for k,v in supplemental.items() if k not in ("source","provenance","lock","process_cpu_s")}
    (root / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
