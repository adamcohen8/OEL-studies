# Cooperative Spacecraft Docking with Time-Limited Readiness Guarantees

**Adam Cohen · Orbital Engineering Lab · November 2026**

OEL Studies original research. [Read the paper](manuscript/ieee-taes/paper.pdf) · [Read the supplement](manuscript/ieee-taes/supplement.pdf)

How can an approaching spacecraft decide whether a docking port will stay ready long enough for contact when the two vehicles have separate clocks? The receiver reports how long it spent processing a readiness request. The chaser uses that measurement, bounds on clock drift, and the remaining distance and closing speed to decide whether it can commit.

Two exchanges take exactly the same time on the chaser’s clock. In one, the receiver’s measurement lets the chaser certify an acceleration of 42 mm/s². In the other, an allowed disturbance leaves the vehicles 10.34 mm apart at the contact deadline even at maximum acceleration. Across 57 simulated encounters, an alternative policy using a known 62.5 ms minimum processing time matched all capture and retreat outcomes. The paper establishes when the measurement can change a guaranteed decision. The frequency of that advantage depends on the encounters and timing information available.

## Verify the results

From this folder, run the following with Python 3.11 and its standard library. The commands work offline.

```sh
python3 -B research/implementation-review/verify_information.py
python3 -B research/implementation-review/verify.py outputs/implementation-review/matched-v2
python3 -B research/peer-revision/verify_tables.py
```

The first check recomputes the two-exchange example and compares the policies profile by profile. It reports zero differences in capture outcome across groups of 21, 24 and 12 profiles. The second independently rescores recorded contact events, retreat clearance, receiver occupancy and actuator onset. The third checks sensitivity results, numerical refinement and the single-host timing experiments.

## Materials

- [Reference scenario](configs/cooperative_docking_paper_replay.yaml), [study flight software](fsw/) and [recorded reference encounter](research/rpo-contact-leases-20261116/runs/evidence-w4/early_valid/output/docking.json)
- [Decision grid](research/rpo-contact-leases-20261116/evidence/decision_grid.csv), [holdout protocol](research/rpo-contact-leases-20261116/HOLDOUT_PROTOCOL.md) and [frozen cases and predictions](research/rpo-contact-leases-20261116/runs/mac-holdout-v1/frozen_plan.json)
- [Closed-loop comparison](outputs/implementation-review/matched-v2/), [known-minimum comparator](outputs/referee-pass-20261118/) and [sensitivity results](outputs/peer-revision/sensitivity-v1/)
- [Paper and supplement sources, figures and bibliography](manuscript/ieee-taes/)

## Rebuild the figures and PDFs

To regenerate the figures and the matched-profile table, install the [analysis dependencies](requirements.txt) in a Python 3.11 environment, then run `python3 -B manuscript/ieee-taes/build_evidence.py` from this folder. The figure PDFs and table are also included. With Tectonic installed, rebuild the documents from those files:

```sh
cd manuscript/ieee-taes
tectonic --untrusted paper.tex
tectonic --untrusted supplement.tex
```

## OEL build and trajectory replay

The [workspace manifest](oel-workspace.yaml) records the runtime version. The trajectories used an experimental OEL checkout reporting version 1.0.0, based on commit `8421e35773876e98c5c09a1e13ffc33d9cdfe4dc` with local source overrides. Supplement S2 identifies the build, scenario and event times. The public holdout plan contains the case inputs and predictions. Its [verification record](research/rpo-contact-leases-20261116/runs/mac-holdout-v1/verification.json) identifies both the public copy and the original frozen plan by digest.

Regenerating the original trajectories requires the author's archived checkout, original frozen plan and compatible runtime. With those materials, set `OEL_DOCKING_ARCHIVED_PLAN` to the original plan path and run `python -B replay.py --engine /path/to/pinned/oel` from this folder. The runner checks the source hashes and compares the complete docking-event sequence. In the reference encounter, first contact and modeled soft capture both occurred at 1.026535727 s. The public checks above operate on the recorded evidence. Physical docking hardware and independent physical clocks remain untested.
