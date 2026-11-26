"""Draw the periapse return map from saved native samples. Run after the study batches."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
import numpy as np


def main():
    root = Path(__file__).resolve().parent
    outputs = root / "outputs"
    fig, ax = plt.subplots(figsize=(8.2, 5.7), layout="constrained")
    def load_run(name):
        path = outputs / name
        return path, json.loads((path / "summary.json").read_text())

    survey, batch = load_run("production-c16")
    for i, track in enumerate(batch["tracks"]):
        with np.load(survey / track["data"]) as data:
            points = data["peri"][1:]
            ax.scatter(*points[:, :2].T, s=6, color="#a0a8ab", alpha=.48, lw=0,
                       label="14-seed survey" if i == 0 else None)

    for name, seed, label, color in [
        ("long-c16", [0,1.5], "384-orbit irregular trajectory", "#b13b24"),
        ("targeted-c16", [0,.5], "Bounded synchronous trajectory", "#114b83"),
    ]:
        directory, batch = load_run(name)
        track = next(t for t in batch["tracks"] if t["seed"] == seed)
        with np.load(directory / track["data"]) as data:
            ax.scatter(*data["peri"][1:, :2].T, s=9, color=color, alpha=.84, lw=0, label=label)

    _, batch = load_run("halfstep-c16")
    roots = [t["seed"] for t in batch["periodic_orbits"] if t["classification"] == "elliptic_linearization"]
    angles, spins = np.array(roots).T
    # Place a symmetry-line point on the displayed side of the modulo-pi seam.
    angles = (angles + np.pi/2) % np.pi - np.pi/2
    ax.scatter(angles, spins, s=54, marker="D", facecolors="none", edgecolors="#152c3b",
               linewidths=1.25, label="Elliptic one-orbit fixed point", clip_on=False)
    supplement = json.loads((outputs / "paper_targets.json").read_text())
    target = next(t for t in supplement["targets"] if t["winding_turns_per_orbit"] == "9/4")
    for row in target["trials"]:
        if row.get("classification") == "elliptic_linearization" and row["status"] == "matched_two_cycle":
            cycle = np.array(row["fine_cycle"]["cycle_theta_w"])
            theta = (cycle[:, 0] + np.pi/2) % np.pi - np.pi/2
            ax.scatter(theta, cycle[:, 1], marker="*", color="#126249", s=98,
                       label="Paper-guided 9/4 two-cycle", clip_on=False)
            break
    ax.set_xlim(-np.pi/2, np.pi/2)
    ax.xaxis.set_major_locator(FixedLocator(np.linspace(-np.pi/2,np.pi/2,5)))
    ax.xaxis.set_major_formatter(FixedFormatter(["−π/2","−π/4","0","π/4","π/2"]))
    ax.set(xlabel="Long-axis angle at periapse (radians, modulo π)", ylabel="Spin / mean orbital motion")
    ax.set_title("Hyperion’s planar spin at successive periapses", loc="left", weight="bold")
    ax.spines[["top","right"]].set_visible(False)
    ax.grid(axis="y", alpha=.17)
    ax.legend(loc="upper center", bbox_to_anchor=(.5,-.17), ncol=2, frameon=False, fontsize=9)
    fig.savefig(root/"return-map.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
