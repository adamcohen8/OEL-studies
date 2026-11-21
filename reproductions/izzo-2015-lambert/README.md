# One-revolution Lambert transfers

Can two different Keplerian ellipses connect the same positions in the same time, in the same direction, after one extra revolution? This workspace reproduces the branch result in section 4.2.2 of Dario Izzo’s *Revisiting Lambert’s Problem* (2015) on an independently specified heliocentric geometry. The paper’s algorithms 1 and 2 give a second formulation for checking the endpoint states. [Paper](https://doi.org/10.1007/s10569-014-9587-y) · [Open manuscript](https://arxiv.org/abs/1403.2705)

## Problem and result

The force is point-mass solar gravity. The departure radius is 1 au and the arrival radius is 1.2 au. Their inertial directions are perpendicular. The specified short-way transfer sweeps 90° after one full revolution in 550 days. The Cartesian inputs, tolerances and step counts are in [problem.json](problem.json).

OEL found a minimum flight time of **479.0704693 days** and two distinct transfers at 550 days. Both trajectories swept **1.25 turns**. Velocities below use the inertial axes in the problem file. The endpoint errors come from the finer independent integration.

| Branch by eccentric-anomaly increment | Departure velocity (km/s) | Arrival velocity (km/s) | Arrival position error |
| --- | --- | --- | ---: |
| Lower phase | (1.780097, 25.261531, 18.946148) | (−26.314095, 2.786177, 2.089633) | 0.063 m |
| Higher phase | (16.316241, 19.419832, 14.564874) | (−20.228992, −9.816354, −7.362265) | 0.218 m |

OEL brackets the minimum of the universal-variable time-of-flight curve and searches each monotone side. The replay integrates each departure state with OEL’s two-body RK4 propagator at 8,192 and 16,384 steps. It compares final position and velocity, and unwraps the in-plane angle along the path to check the revolution count. A separate calculation evaluates Izzo’s equation (18) and reconstructs the velocities using algorithm 1. The reference uses a different root parameter and pairs solutions by their velocities.

The [verified replay output](result.json) includes full-precision vectors and errors. Both step counts met the preset 1 km position, 10⁻⁶ km/s velocity and 10⁻⁸ turn tolerances. The independent formulations agreed within 2 × 10⁻¹² km/s on all endpoint velocity vectors. This case tests the two-branch result for one specified geometry. The paper’s large random-sample accuracy, iteration-count and timing benchmarks are outside its scope.

## Replay

The workspace requires OEL **1.2.0+preview.1**. Its exact revision and runtime source digests are in [engine-revision.json](engine-revision.json). The baseline OEL 1.1.0 release supplies dependency wheels for the tested Apple Silicon / CPython 3.11 installation. The OEL engine itself is built from the pinned preview source. [Public OEL preview](https://github.com/adamcohen8/orbital-engineering-lab/tree/5a3d1eaf91775a304efccb85274571700578a0cf)

Run these commands in an empty installation directory with Git, curl, unzip and CPython 3.11 available.

```sh
git init -q engine
git -C engine fetch -q --depth=1 https://github.com/adamcohen8/orbital-engineering-lab.git 5a3d1eaf91775a304efccb85274571700578a0cf
git -C engine checkout -q --detach FETCH_HEAD
curl --fail --location --output deps.zip https://github.com/adamcohen8/orbital-engineering-lab/releases/download/v1.1.0/oel-public-1.1.0-arm64-py311.bundle.zip
printf '%s  deps.zip\n' 355a09916cfccb4ef3eade15e799040e785c7ea57b1e8612bce70784f3a1d20e | shasum -a 256 -c -
unzip -q deps.zip 'wheelhouse/*.whl'
python3.11 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links "$PWD/wheelhouse" -c engine/constraints/py311.txt ./engine
.venv/bin/python -m pip check
```

Download this study folder, retaining its files together. From the installation directory, substitute the path to the downloaded folder and a fresh output filename.

```sh
.venv/bin/python /path/to/izzo-2015-lambert/reproduce.py --output replay.json
```

A successful replay exits with status zero and writes both branches, checks at both step counts, independently calculated reference states and a passing comparison. The replay runs offline. [Inputs](problem.json) · [Replay](reproduce.py) · [Reference calculation](reference.py)
