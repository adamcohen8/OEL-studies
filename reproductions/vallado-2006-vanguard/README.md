# Vanguard 1 propagation over 72 hours

Vallado, Crawford, Hujsak and Kelso distributed reference trajectories with *Revisiting Spacetrack Report #3* (AIAA 2006-6753). This study reproduces the historical Vanguard 1 verification case, catalog 00005. The two-line element set describes an epoch on June 27, 2000. The comparison covers 13 native TEME position and velocity vectors spaced six hours apart. [Paper and verification materials](https://celestrak.org/publications/AIAA/2006-6753/)

## Results

The acceptance limits, set before the run, were 1 m for the Euclidean position difference and 1 mm/s for the Euclidean velocity difference at every epoch. Both comparisons use public OEL 1.1.0 and the same input.

| Trajectory | Maximum position difference | Maximum velocity difference | Verdict |
| --- | ---: | ---: | --- |
| Analyst’s submitted Python backend | 7.891 m at 6 h | 6.263 mm/s at 66 h | Exceeded both limits |
| Reviewer’s subsequent comparison of the default Rust backend | 0.00000682 m | 0.0007705 mm/s | Within both limits at all 13 epochs |

GPT-6.2 Sol designed and wrote the replay, corrected its initial assumptions about the API after a failed run, and chose a finite-difference velocity check, reverse-order replay and a Python/Rust comparison. Sol recognized the historical test input. It saw its own computed trajectory and both sets of diagnostics, with the published vectors withheld until it had submitted its finding. Sol flagged the disagreement and retained the Python trajectory it had originally selected. The reviewer then compared both paths against the published data. Agreement here concerns numerical reproduction of one reference trajectory. It provides no estimate of Vanguard’s actual orbit uncertainty. [Analyst’s method](METHOD.md) · [Finding as submitted before scoring](FINDINGS.md)

## Inputs and scoring

The supplied TLE is the 69-character orbit record from the authors’ input file. The sampling interval comes from the trailing test controls on its second line. The archived ephemeris gives elapsed seconds, TEME kilometers, and kilometers per second. The scorer converts its epochs to minutes before comparing all six coordinates. The ephemeris header declares 146 points, while the actual data block contains the 13 expected epochs. The scorer checks the rows against the prescribed grid. [Pinned input collection](https://github.com/CelesTrak/fundamentals-of-astrodynamics/blob/faa0d45f504e51f9d0b946594c39f18b311ab992/datalib/SGP4-VER.TLE) · [Pinned reference ephemeris](https://github.com/CelesTrak/fundamentals-of-astrodynamics/blob/faa0d45f504e51f9d0b946594c39f18b311ab992/software/cpp/TestSGP4/TestSGP4/00005.e)

## Install and verify

The tested runtime was CPython 3.11 on Apple Silicon. From the downloaded study folder, the following commands install the pinned public wheel with its bundled dependencies. The network is used only to obtain the signed release bundle. [Public OEL 1.1.0 release](https://github.com/adamcohen8/orbital-engineering-lab/releases/tag/v1.1.0)

```sh
curl -fLO https://github.com/adamcohen8/orbital-engineering-lab/releases/download/v1.1.0/oel-public-1.1.0-arm64-py311.bundle.zip
echo '355a09916cfccb4ef3eade15e799040e785c7ea57b1e8612bce70784f3a1d20e  oel-public-1.1.0-arm64-py311.bundle.zip' | shasum -a 256 -c -
unzip -q oel-public-1.1.0-arm64-py311.bundle.zip -d release-bundle
python3.11 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links release-bundle/wheelhouse orbital-engineering-lab==1.1.0
.venv/bin/python -I replay.py
.venv/bin/python -I grade.py
```

The replay generates `states.csv`. The separate scorer reads those rows, verifies the archived reference file by SHA-256, then computes Euclidean position and velocity differences on the frozen epoch grid. Its exit code is zero when both tolerances hold at every epoch and one for a numerical mismatch.

For the committed Python trajectory the grader prints `"passed": false` and exits 1. To reproduce the reviewer’s second comparison, run `.venv/bin/python -I check_default.py`. It generates a temporary trajectory using OEL’s default backend and passes it to the same scorer. The expected exit code is zero. Both scoring results are included in this folder.
