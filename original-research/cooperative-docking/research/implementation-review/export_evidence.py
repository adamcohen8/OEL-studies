"""Export only the archived measurements needed for independent score checks."""
import argparse
import gzip
import json
from pathlib import Path

from run_comparison import ROOT, QUERY, diagnostic_frames, digest
from sim.review import ReviewWorkspace


def records(batch, new=False):
    for row in json.loads((batch / "runs.json").read_text()):
        folder = batch / row.get("folder", row["case"] + "_" + row["policy"])
        result = folder / "result"
        assert digest(folder / "config.yaml") == row["config_sha256"]
        assert digest(result / "review/run.sqlite") == row["review_sha256"]
        data = {"case": row["case"], "policy": row["policy"], "refinement": row.get("refinement", 1),
                "config_sha256": row["config_sha256"], "review_sha256": row["review_sha256"]}
        review = ReviewWorkspace.open(result)
        events = json.loads((result / "docking.json").read_text())["events"]
        data["events"] = [[e["event"], e["time_s"], e.get("capture_armed"),
                           e.get("normal_speed_m_s")] for e in events]
        sample = review.query("SELECT time_s, gap_m, angle_rad, normal_speed_m_s "
                              "FROM docking_samples ORDER BY time_s", max_rows=10000)
        assert sample.rows and not sample.truncated
        data["end_s"] = sample.rows[-1]["time_s"]
        if not any(e["event"] == "first_contact" for e in events):
            data["retreat_samples"] = [[s["time_s"], s["gap_m"], s["angle_rad"],
                                        s["normal_speed_m_s"]] for s in sample.rows]
        if new:
            diagnostic = review.query(QUERY, max_rows=5000, max_vm_steps=750000)
            assert not diagnostic.truncated
            chaser = diagnostic_frames(diagnostic.rows, "chaser")
            receiver = diagnostic_frames(diagnostic.rows, "target")
            data["grants"] = {}
            for t, f in receiver:
                if f["issued_grant"]:
                    data["grants"].setdefault(f["issued_grant"], [t, f["issued_duration_ns"]])
            data["requests"] = max(f["requests"] for t, f in chaser)
            decision_ns, decision = next((t, f) for t, f in chaser if
                                         f["mode"] == "commit" or f["mode"].startswith("abort"))
            data["decision_ns"] = decision_ns
            data["mode"] = decision["mode"]
            if decision["mode"] == "commit":
                data["certificate"] = {k: decision[k] for k in
                                       ("chosen_grant", "earliest_s", "latest_s", "contact_speed_min_m_s",
                                        "contact_speed_max_m_s", "selected_acceleration_m_s2")}
            if row["case"].startswith("onset_"):
                thrust = review.query("SELECT time_s, accel_z_eci_km_s2 FROM thrust "
                                      "WHERE object_id='chaser' ORDER BY time_s", max_rows=1000)
                assert not thrust.truncated
                data["axial_command"] = [[t, f["commanded_axial_acceleration_m_s2"]] for t, f in chaser]
                data["thrust"] = [[s["time_s"], s["accel_z_eci_km_s2"]] for s in thrust.rows]
        yield data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch", type=Path)
    batch = parser.parse_args().batch
    sources = {"original": ROOT / "outputs/adaptive-comparator/matched-v1",
               "revision": ROOT / "outputs/peer-revision/sensitivity-v1", "current": batch}
    data = {name: list(records(source, name == "current")) for name, source in sources.items()}
    (batch / "observations.json.gz").write_bytes(gzip.compress(json.dumps(data).encode(), mtime=0))


if __name__ == "__main__":
    main()
