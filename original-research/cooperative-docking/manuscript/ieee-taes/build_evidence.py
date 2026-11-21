"""Rebuild TAES figures and supplement rows from the verified OEL batches."""

import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
WORKSPACE = HERE.parents[1]
HOLDOUT = WORKSPACE / "research/rpo-contact-leases-20261116/runs/mac-holdout-v1"
COMPARATOR = WORKSPACE / "outputs/adaptive-comparator/matched-v1"
BLUE = "#185477"
ORANGE = "#a14825"
GRAY = "#535f69"
POLICIES = ("fixed_one_shot", "fixed_adaptive", "fixed_renew_coast", "corrected_adaptive")


def read_json(path, expected_digest=None):
    data = path.read_bytes()
    if expected_digest is not None and hashlib.sha256(data).hexdigest() != expected_digest:
        raise ValueError(f"Digest mismatch for {path}")
    return json.loads(data)


def plot_admission():
    """Display the scalar certificate and the measured information ablation."""
    import csv
    import sys
    import numpy as np

    sys.path.insert(0, str(WORKSPACE / "fsw"))
    from certificate import certify
    from evaluate import HANDOFF_S, approach
    sys.path.insert(0, str(WORKSPACE / "research/implementation-review"))
    from verify_information import grid_minimum_admits

    folder = WORKSPACE / "research/rpo-contact-leases-20261116/evidence"
    with (folder / "decision_grid.csv").open() as stream:
        grid = list(csv.DictReader(stream))
    policies = (
        ("Fixed round trip", ("fixed_round_trip",)),
        ("Measured round trip", ("measured_round_trip",)),
        ("Fixed reverse", ("fixed_reverse",)),
        ("Fixed + chaser elapsed", ("fixed_reverse", "measured_round_trip")),
        ("Turnaround alone", ("turnaround_corrected",)),
        ("Fixed + turnaround", ("combined",)),
        ("Exact-time oracle", ("oracle",)),
    )
    totals = [sum(any(row[k + "_admits"] == "True" for k in keys) for row in grid)
              for _, keys in policies]
    latest = read_json(WORKSPACE / "outputs/implementation-review/matched-v2/verification.json")
    assert len(grid) == 1080 and totals[3] == latest["grid"]["elapsed_plus_fixed"]
    labels = [label for label, _ in policies]
    labels.insert(4, "Fixed + elapsed + 31.25 ms min.")
    totals.insert(4, sum(map(grid_minimum_admits, grid)))
    gaps, times = np.linspace(20, 90, 141), np.linspace(.40, 1.85, 146)
    values = np.full((len(times), len(gaps)), np.nan)
    for j, gap in enumerate(gaps):
        state = approach(gap / 1000, .040)
        for i, time in enumerate(times):
            result = certify(state, time - HANDOFF_S, times=False)
            if result:
                values[i, j] = 1000 * result.acceleration
    fig, (ax, bars) = plt.subplots(1, 2, figsize=(8.8, 3.25), layout="constrained",
                                   gridspec_kw={"width_ratios": (1.13, 1)})
    ax.set_facecolor("#eeeeee")
    mesh = ax.pcolormesh(gaps, times, values, shading="nearest", cmap="viridis",
                         vmin=0, vmax=50, rasterized=True)
    ax.set(xlabel="Estimated port gap (mm)", ylabel="Guaranteed readiness remaining (s)")
    ax.set_title("(a) Admission at 40 mm/s", loc="left")
    ax.text(72, .60, "Decline", color=GRAY, ha="center")
    fig.colorbar(mesh, ax=ax, pad=.025, label="Selected acceleration (mm/s²)")
    pos = np.arange(len(totals))
    bars.barh(pos, totals, height=.65, color=["#8b9ca6", "#8b9ca6", BLUE, BLUE, BLUE,
                                             "#167477", "#167477", "#262b2e"])
    bars.set_yticks(pos, labels)
    bars.invert_yaxis()
    bars.set(xlim=(0, 625), xlabel="Admitted cases (out of 1,080)")
    bars.set_title("(b) Timing policy comparison", loc="left")
    for i, value in enumerate(totals):
        bars.text(value + 7, i, str(value), va="center")
    fig.savefig(HERE / "admission.pdf", dpi=250)
    plt.close(fig)


def plot_timeline():
    fig, ax = plt.subplots(figsize=(7.1, 1.84))
    fig.subplots_adjust(left=0.02, right=0.985, bottom=0.03, top=0.98)
    ax.set(xlim=(0, 1), ylim=(0, 4.5))
    ax.set_axis_off()
    sent, received, start, reply, decision = 0.045, 0.198, 0.341, 0.490, 0.612
    earliest, latest, handoff, guaranteed, expiry = 0.689, 0.787, 0.857, 0.904, 0.974

    def span(left, right, y, label, color, *, offset=0.19, align="bottom"):
        ax.hlines(y, left, right, color=color, lw=0.85)
        ax.vlines((left, right), y - 0.10, y + 0.10, color=color, lw=0.85)
        ax.text((left + right) / 2, y + offset, label, ha="center", va=align, color=color)

    span(sent, decision, 3.94, r"$E_C$ (chaser clock)", BLUE, offset=0.13)
    span(received, start, 3.20, r"$P_R$ (receiver clock)", ORANGE, offset=0.13)
    span(start, expiry, 3.20, r"$H_R$ (receiver clock)", ORANGE, offset=0.13)
    ax.hlines(2.24, sent, 0.995, color=GRAY, lw=0.8)
    ax.annotate("", (0.995, 2.24), (0.984, 2.24),
                arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 0.8})
    ax.text(0.995, 2.43, "Physical time", ha="right", va="bottom", color=GRAY)
    for x, height, name in (
        (sent, 3.94, "Request\nsent"),
        (received, 3.20, "Request\nreceived"),
        (start, 3.20, "Grant\nstarts"),
        (reply, 2.43, "Reply\nreceived"),
        (decision, 3.94, "Chaser\ndecides"),
        (expiry, 3.20, "Actual\nexpiry"),
    ):
        ax.vlines(x, 2.22, height, color="#abb6bd", lw=0.6, linestyles="dotted")
        ax.plot(x, 2.24, "o", ms=2.6, color=GRAY)
        ax.text(x, 2.07, name, ha="center", va="top", linespacing=1.05)
    ax.vlines(guaranteed, 0.86, 2.24, color=BLUE, lw=0.65, linestyles="dashed")
    ax.vlines(decision, 0.91, 1.43, color=BLUE, lw=0.65, linestyles="dotted")
    span(decision, guaranteed, 1.02, r"$R$ (guaranteed remaining)", BLUE, offset=0.13)
    span(earliest, latest, 0.43, "", BLUE)
    span(latest, handoff, 0.43, r"$q$", ORANGE, offset=-0.10, align="top")
    ax.text(earliest - 0.013, 0.43, "Certified contact window", ha="right", va="center",
            color=BLUE)
    fig.savefig(HERE / "readiness-timeline.pdf")
    plt.close(fig)


def plot_holdout():
    verified = read_json(HOLDOUT / "verification.json")
    plan = read_json(HOLDOUT / "frozen_plan.json", verified["frozen_plan_sha256"])
    rows = read_json(HOLDOUT / "results.json", verified["oel_results_sha256"])
    cases = {record["case"]["name"]: record for record in plan["cases"]}
    results = {row["case"]: row for row in rows}
    assert len(cases) == len(rows) == len(results) == 36
    assert cases.keys() == results.keys()
    assert not any(row["score_failures"] for row in rows)

    fig, (timing, geometry) = plt.subplots(1, 2, figsize=(7.1, 2.66))
    fig.subplots_adjust(left=0.079, right=0.982, wspace=0.40, bottom=0.32, top=0.90)

    for pair in range(1, 9):
        inside, outside = (f"timing_{pair:02d}_{side}" for side in ("inside", "outside"))
        margin_ms = [1e3 * cases[key]["horizon_displacement_from_boundary_s"]
                     for key in (outside, inside)]
        assert margin_ms[0] < 0 < margin_ms[1]
        assert [results[key]["outcome"] for key in (outside, inside)] == ["retreat", "capture"]
        timing.plot(margin_ms, (pair, pair), color="#b9c2c8", lw=0.85)
        timing.plot(margin_ms[0], pair, marker="x", ls="none", color=ORANGE, ms=5, mew=1.1,
                    label="Retreat" if pair == 1 else None)
        timing.plot(margin_ms[1], pair, marker="o", ls="none", color=BLUE, ms=4.5,
                    label="Soft capture" if pair == 1 else None)
    timing.axvline(0, color=GRAY, ls="--", lw=0.75)
    timing.set(title="(a) Readiness boundary", xlabel="Horizon minus admission threshold (ms)",
               ylabel="Timing profile", yticks=range(1, 9), ylim=(0.4, 8.6), xlim=(-18, 18),
               xticks=(-15, -10, -5, 0, 5, 10, 15))
    timing.invert_yaxis()
    timing.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.27),
                  handletextpad=0.25, columnspacing=1.15)

    models = (("ring", "Finite radius", "o", BLUE), ("center", "Center only", "s", ORANGE))
    early_misses = 0
    for model, label, marker, color in models:
        for bracket, dx in (("lower", -0.14), ("higher", 0.14)):
            offsets = []
            for profile in range(1, 6):
                key = f"geometry_{profile:02d}_{bracket}_{model}"
                record, row = cases[key], results[key]
                frozen = record["prediction"]
                lower = frozen["decision_s"] + frozen["certificate"]["earliest"]
                upper = frozen["decision_s"] + frozen["certificate"]["latest"]
                assert row["outcome"] == "capture" and row["contact_s"] <= upper
                difference_ms = (row["contact_s"] - lower) * 1e3
                assert (difference_ms < 0) == record["planned_early_miss"]
                early_misses += difference_ms < 0
                offsets.append(difference_ms)
            geometry.plot([i + dx for i in range(1, 6)], offsets, ls="none",
                          color=color, marker=marker, ms=5.1, mew=1,
                          mfc="white" if bracket == "lower" else color,
                          label=f"{label}, {bracket}")
    assert early_misses == 5
    geometry.axhline(0, color=GRAY, ls="--", lw=0.75)
    geometry.set(title="(b) Geometry boundary", xlabel="Geometry profile",
                 ylabel="Contact minus earliest bound (ms)", xticks=range(1, 6),
                 xlim=(0.6, 5.4), ylim=(-33, 137), yticks=(-25, 0, 25, 50, 75, 100, 125))
    handles, labels = geometry.get_legend_handles_labels()
    geometry.legend(handles, labels, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.27),
                    handletextpad=0.25, columnspacing=0.85, fontsize=8)

    fig.savefig(HERE / "boundary-holdout.pdf")
    plt.close(fig)


def contact_cell(row, highlight=False):
    time = row["first_contact_s"]
    if row["outcome"] == "capture":
        assert time is not None and row["commit_grant"] in (1, 2)
        label = f"{time:.3f}"
        if highlight:
            label = rf"\textbf{{{label}}}"
    else:
        assert row["outcome"] == "retreat" and time is None and row["commit_grant"] is None
        label = "R"
    return f"{label} ({row['requests']})"


def write_comparator_rows():
    verified = read_json(COMPARATOR / "verification.json")
    plan = read_json(COMPARATOR / "plan.json", verified["plan_sha256"])
    rows = read_json(COMPARATOR / "runs.json", verified["results_sha256"])
    names = [case["name"] for case in plan["cases"]]
    results = {(row["case"], row["policy"]): row for row in rows}
    assert len(rows) == len(results) == verified["cases_verified"] == 84
    assert len(names) == len(set(names)) == 21
    assert set(results) == {(name, policy) for name in names for policy in POLICIES}
    assert not any(row["failures"] for row in rows)

    current = WORKSPACE / "outputs/implementation-review/matched-v2"
    identity = read_json(current / "identity.json")
    elapsed = {r["case"]: r for r in read_json(current / "runs.json", identity["runs_sha256"])
               if r["case"] in names and r["policy"] == "elapsed_adaptive"}
    assert len(elapsed) == 21 and not any(r["failures"] for r in elapsed.values())
    special = {
        "first_reply_lost": "First reply lost",
        "renewal_reply_lost": "Renewal reply lost",
        "slow_processing": "Extended turnaround",
    }
    counts = Counter()
    lines = ["% Generated from verified OEL results by build_evidence.py."]
    for index, case in enumerate(plan["cases"]):
        name = case["name"]
        cells = [results[name, policy] for policy in POLICIES]
        fixed, corrected = cells[1:3], cells[3]
        fixed_times = [r["first_contact_s"] for r in fixed if r["outcome"] == "capture"]
        contact = corrected["first_contact_s"]
        gain = contact is not None and not fixed_times
        if gain:
            cause = "cap" if case["second_cap_s"] < 1.9 else "loss"
            assert cause == "cap" or case["drop"][1]
            counts[cause] += 1
        delta = None if contact is None or not fixed_times else (min(fixed_times) - contact) * 1e3
        delta_label = "0" if delta == 0 else f"{delta:+.0f}" if delta is not None else "--"
        if name in special:
            profile, delay, cap = special[name], "--", f"{case['second_cap_s']:.2f}"
        else:
            group, delay, _ = name.split("_")
            profile, delay, cap = group.title(), delay.title(), f"{case['second_cap_s']:.2f}"
        if index in (6, 12, 18):
            lines.append(r"\addlinespace[2pt]")
        fields = [profile, delay, cap]
        cells.insert(3, elapsed[name])
        fields.extend(contact_cell(row, highlight=gain and row is corrected) for row in cells)
        fields.append(delta_label)
        lines.append(" & ".join(fields) + r" \\")

    assert counts == {"cap": 4, "loss": 1}
    assert [sum(results[name, p]["outcome"] == "capture" for name in names) for p in POLICIES] == [
        10, 15, 15, 20,
    ]
    lines.append(r"\bottomrule")
    (HERE / "comparator-rows.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.7,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "legend.frameon": False,
        "pdf.fonttype": 42,
    })
    plot_admission()
    plot_timeline()
    plot_holdout()
    write_comparator_rows()
