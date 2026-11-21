"""Rescore physical evidence and paired outcomes using the standard library."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
import json
from math import fsum, isclose, sin
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
POLICIES = ("fixed_adaptive", "elapsed_adaptive", "corrected_adaptive")


def load(path):
    return json.loads(path.read_bytes())


def historical(row, obs):
    contact = next((e for e in obs["events"] if e[0] == "first_contact"), None)
    capture = next((e for e in obs["events"] if e[0] == "soft_capture"), None)
    assert not row["failures"]
    assert row["review_sha256"] == obs["review_sha256"] and row["config_sha256"] == obs["config_sha256"]
    if row["outcome"] == "capture":
        assert contact and capture and contact[1] == row["first_contact_s"] and contact[2]
        return None, capture[1] == contact[1]
    assert not contact and not capture and row["outcome"] == "retreat"
    samples = obs["retreat_samples"]
    center = min(g for _, g, _, _ in samples)
    edge = min(g - .2 * sin(angle) for _, g, angle, _ in samples)
    assert center == row["min_gap_m"] and edge >= .002 and samples[-1][-1] > 0
    return (center, edge), False


def score(row, obs):
    grants = {int(k): (start / 1e9, start / 1e9 + duration / 1e9 / .999260)
              for k, (start, duration) in obs["grants"].items()}
    occupied = end = 0.
    for start, stop in sorted(grants.values()):
        occupied += max(0., stop - max(end, start))
        end = max(end, stop)
    assert isclose(occupied, row["receiver_occupied_s"], abs_tol=1e-12)
    assert obs["requests"] == row["requests"]
    if row["outcome"] == "retreat":
        assert obs["mode"].startswith("abort")
        return
    assert obs["mode"] == "commit"
    c = obs["certificate"]
    contact = next(e for e in obs["events"] if e[0] == "first_contact")
    captured = next(e for e in obs["events"] if e[0] == "soft_capture")
    t = obs["decision_ns"] / 1e9
    assert t + c["earliest_s"] - 1e-8 <= contact[1] <= t + c["latest_s"] + 1e-8
    assert c["contact_speed_min_m_s"] - 1e-8 <= -contact[3] <= c["contact_speed_max_m_s"] + 1e-8
    assert c["chosen_grant"] == row["commit_grant"]
    assert grants[c["chosen_grant"]][1] >= t + c["latest_s"] + .04 - 1e-8
    assert contact[1] <= captured[1] <= contact[1] + .04 <= obs["end_s"]
    assert not any(e[0] in ("rebound", "latch_release") and e[1] <= contact[1]+.04 for e in obs["events"])


def onset(row, obs):
    capture = row["policy"] == "corrected_adaptive"
    assert row["outcome"] == ("capture" if capture else "retreat")
    acceleration = obs["certificate"]["selected_acceleration_m_s2"] if capture else -.08
    assert acceleration != 0
    t = next(t for t, a in obs["axial_command"] if a == acceleration)
    # OEL stores the applied acceleration at the right endpoint of its step.
    start, stop = next((a, b) for a, b in zip(obs["thrust"], obs["thrust"][1:])
                       if b[1] is not None and isclose(-1000*b[1], acceleration, abs_tol=1e-12))
    assert start[0] == t / 1e9 and start[1] == 0 and stop[0] - start[0] == .03125
    delay = t - obs["decision_ns"]
    nominal = int(row["case"].split("_")[1])
    assert delay == ((nominal + 31249999) // 31250000) * 31250000
    assert 31250000 <= delay <= 62500000
    return delay


def verify(batch):
    observations = json.loads(gzip.decompress((batch / "observations.json.gz").read_bytes()))
    for key, expected in load(batch / "identity.json").items():
        assert sha256((batch / (key.removesuffix("_sha256") + ".json")).read_bytes()).hexdigest() == expected
    source = {"original": ROOT / "outputs/adaptive-comparator/matched-v1",
              "revision": ROOT / "outputs/peer-revision/sensitivity-v1", "current": batch}
    original = {}
    groups = defaultdict(lambda: defaultdict(list))
    margins = {}
    delays = {}
    repeats = 0
    immediate = {}
    for name, root in source.items():
        rows = load(root / "runs.json")
        lookup = {(r["case"], r["policy"], r.get("refinement", 1)): r for r in rows}
        ids = [(o["case"], o["policy"], o["refinement"]) for o in observations[name]]
        assert len(ids) == len(set(ids)) == len(rows) and set(ids) == set(lookup)
        gaps = []
        immediate[name] = 0
        for key, obs in zip(ids, observations[name]):
            row = lookup[key]
            margin, same_time = historical(row, obs)
            immediate[name] += same_time
            if margin:
                gaps.append(margin)
            if name != "current":
                if key[2] == 1:
                    original[key[:2]] = row
                continue
            score(row, obs)
            if margin:
                assert margin[1] == row["min_sampled_edge_gap_m"]
            if key[0].startswith("onset_"):
                delays["/".join(key[:2])] = onset(row, obs)
                continue
            if key[:2] in original:
                reference = original[key[:2]]
                for field in ("outcome", "requests", "commit_grant", "first_contact_s", "abort_s"):
                    assert reference[field] == row[field]
                repeats += 1
            cohort = "sensitivity" if key[0].startswith("s_") else "broader" if key[0].startswith("b_") else "matched"
            groups[cohort][row["policy"]].append(row)
        margins[name] = dict(retreats=len(gaps), center_m=min(x[0] for x in gaps),
                            edge_m=min(x[1] for x in gaps))
    grouped = {}
    wins = {}
    for cohort, policies in groups.items():
        grouped[cohort] = {}
        for p, values in policies.items():
            grouped[cohort][p] = dict(n=len(values), captures=sum(v["outcome"]=="capture" for v in values),
                                     first=sum(v["commit_grant"]==1 for v in values),
                                     renewals=sum(v["requests"]-1 for v in values),
                                     occupancy_s=fsum(v["receiver_occupied_s"] for v in values))
        wins[cohort] = {}
        for p, q in zip(POLICIES, POLICIES[1:]):
            a,b = ({r["case"] for r in policies[k] if r["outcome"]=="capture"} for k in (p,q))
            wins[cohort][q+"_over_"+p] = {"added": sorted(b-a), "lost": sorted(a-b)}
    grid = ROOT / "research/rpo-contact-leases-20261116/evidence/decision_grid.csv"
    with grid.open() as stream:
        rows = list(csv.DictReader(stream))
    counts = Counter()
    for r in rows:
        flags = [r[k+"_admits"]=="True" for k in ("fixed_reverse","measured_round_trip","turnaround_corrected","combined")]
        counts["fixed"] += flags[0]
        counts["elapsed_plus_fixed"] += flags[0] or flags[1]
        counts["turnaround_plus_fixed"] += flags[3]
        counts["turnaround_only_gain"] += flags[2] and not flags[0]
        counts["turnaround_only_loss"] += flags[0] and not flags[2]
    assert repeats == 114 and len(delays) == 10 and len(rows) == 1080
    return dict(groups=grouped, paired=wins, realized_onset_ns=delays, margins=margins,
                simultaneous_contact_capture=immediate, matched_archival_repeats=repeats, grid=counts)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = verify(args.batch)
    summary = args.batch / "verification.json"
    if args.write:
        summary.write_text(json.dumps(result, indent=2) + "\n")
    else:
        assert result == load(summary)
    print(json.dumps(result, indent=2))
