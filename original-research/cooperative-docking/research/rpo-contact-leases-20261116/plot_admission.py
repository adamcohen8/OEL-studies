"""Plot the analytical admission region and matched workload coverage."""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from certificate import certify
from evaluate import HANDOFF_S, approach

ROOT = Path(__file__).resolve().parent


def main():
    summary = json.loads((ROOT / "evidence/analytical_verification.json").read_text())
    gaps = np.linspace(20, 90, 141)
    windows = np.linspace(0.40, 1.85, 146)
    values = np.full((len(windows), len(gaps)), np.nan)
    for i, window in enumerate(windows):
        for j, gap in enumerate(gaps):
            cert = certify(approach(gap / 1000, 0.040), window - HANDOFF_S, times=False)
            if cert:
                values[i, j] = 1000 * cert.acceleration

    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42})
    fig, (ax, bars) = plt.subplots(1, 2, figsize=(9.1, 3.4), layout="constrained",
                                   gridspec_kw={"width_ratios": (1.24, 1)})
    ax.set_facecolor("#eeeeee")
    image = ax.pcolormesh(gaps, windows, values, shading="nearest", cmap="viridis", vmin=0, vmax=50)
    ax.set(xlabel="Estimated port gap (mm)", ylabel="Guaranteed readiness remaining (s)")
    ax.set_title("(a) Admission at 40 mm/s", loc="left", pad=8)
    ax.text(72, 0.60, "Decline", color="#4c5257", ha="center")
    legend = fig.colorbar(image, ax=ax, pad=0.025)
    legend.set_label("Selected acceleration (mm/s²)")

    keys = ("fixed_round_trip", "measured_round_trip", "fixed_reverse",
            "turnaround_corrected", "combined", "oracle")
    labels = ("Fixed round trip", "Measured round trip", "Fixed reverse bound",
              "Corrected round trip", "Corrected + fixed reverse", "Exact-time oracle")
    heights = [summary["admissions"][key + "_admits"] for key in keys]
    colors = ["#8b9ca6", "#8b9ca6", "#4d7187", "#167477", "#167477", "#262b2e"]
    positions = np.arange(len(keys))
    bars.barh(positions, heights, color=colors, height=0.64)
    bars.set_yticks(positions, labels)
    bars.invert_yaxis()
    bars.set(xlim=(0, 625), xlabel="Admitted cases (out of 1,080)")
    bars.set_title("(b) Timing policy comparison", loc="left", pad=8)
    bars.set_axisbelow(True)
    bars.grid(axis="x", linewidth=0.5, alpha=0.25)
    for row, count in enumerate(heights):
        bars.text(count + 8, row, str(count), va="center")

    figures = ROOT / "figures"
    figures.mkdir(exist_ok=True)
    for suffix in ("pdf", "png"):
        fig.savefig(figures / f"admission.{suffix}", dpi=230)


if __name__ == "__main__":
    main()
