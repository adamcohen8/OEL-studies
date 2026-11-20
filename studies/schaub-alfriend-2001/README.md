# Schaub & Alfriend (2001) — relative phase drift over 45 orbits

Does applying the paper’s orbital-element offsets in mean elements reduce the long-term drift between two passive spacecraft? The target is the 48° example in section 4.1, Figures 3B and 4B, of *J₂ Invariant Relative Orbits for Spacecraft Formations*, Celestial Mechanics and Dynamical Astronomy 79, 77–95. The curves show approximately +0.0053° of drift for offsets applied directly to osculating elements and −0.000019° after mean-element initialization. Those values are visual estimates. [Author’s copy](https://hanspeterschaub.info/PapersPrivate/MeanElemKluv.pdf)

## Replay

The workspace locks OEL **1.1.0**. The verified signed installation uses macOS on Apple Silicon and CPython 3.11. Install the core profile from the pinned release. If several Python launchers are installed, place the real CPython 3.11 provider directory on PATH for the installer. [Release and platform instructions](https://github.com/adamcohen8/orbital-engineering-lab/releases/tag/v1.1.0)

```sh
release=https://github.com/adamcohen8/orbital-engineering-lab/releases/download/v1.1.0
curl -fLsS "$release/install.sh" -o oel-install-1.1.0.sh
echo 'c32689b3a1d99c14848f628b8306bd769d97b0c63f0f410410e25b2a25e3f305  oel-install-1.1.0.sh' | shasum -a 256 -c -
# Inspect the downloaded script before running it.
OEL_INSTALL_PROFILE=core OEL_UPDATE_CHANNEL_URL="$release/public-stable.json" sh oel-install-1.1.0.sh
```

Download this folder, or use Git to select it. Then run the replay with Python 3 and the installed OEL launcher.

```sh
git clone --depth 1 --filter=blob:none --sparse https://github.com/adamcohen8/OEL-studies.git
cd OEL-studies
git sparse-checkout set studies/schaub-alfriend-2001
cd studies/schaub-alfriend-2001
python3 reproduce.py
```

If the launcher is outside PATH, use `python3 reproduce.py --oel /absolute/path/to/oel`. The script validates each YAML, runs OEL, checks the accepted run identity and complete time grid, and writes the measured drifts and pass flags in a fresh `outputs/replay-*/comparison.json`. Both runs use only OEL and Python’s standard library.

## Initial states and dynamics

In the published case the chief starts at 7,153 km semimajor axis, eccentricity 0.05, inclination 48°, ascending node 0°, argument of perigee 30° and mean anomaly 0°. The deputy offsets are −0.351765 m, +0.0001, +0.001035°, +0.005°, +0.01° and −0.01° in the same order. The variation changes the chief’s initial argument of perigee to 75° and keeps those offsets. It is a sensitivity check devised for this reproduction.

Each YAML supplies osculating elements to the OEL Numerical Propagator. The direct pair applies the offsets in osculating elements. The mapped pair applies them in mean elements and converts both spacecraft using the first-order Brouwer mapping. The optional `initial_elements.py` regenerates the initial elements with Basilisk 2.12.1. The paper supplies the offsets and the replay starts from the prepared elements. [Basilisk orbital-element conversion](https://hanspeterschaub.info/basilisk/Documentation/utilities/orbitalMotion.html)

OEL propagates both pairs in an inertial frame with point-mass gravity and static Earth zonals J₂ through J₅. Each YAML lists the unnormalized coefficients. The reproduction uses a 6,378.1366 km reference radius and J₂ of 0.001082616. The paper quotes rounded values of 6,378.2 km and 0.001082. The assumed epoch is JD 2451545.0 UTC. The nominal period is rounded to 6,021 s, the RK4 substep is 6,021/1,024 s and the output cadence is 6,021/64 s.

## Measurement and comparison

The review query in `reproduce.py` retrieves osculating elements at 2,881 times per spacecraft. For each pair, compute deputy minus chief in ω + M, wrapping the difference at ±180°. Average the first 2,880 differences in 45 equal one-orbit bins. Fit a line to the bin means and multiply its slope by 45 to measure the secular change in degrees. The final sample verifies coverage. The reduction factor is the ratio of absolute drifts.

An independently fitted reference reconstructed mean elements from simulated Cartesian states. Replay acceptance tolerances relative to that reference are 5% for the direct drift, 20% for the smaller mapped residual and 25% for the reduction factor. Both signed changes must agree in sign. The paper supplies a graphical comparison for the 30° case only.

| Case | Reference direct (°) | Reference mapped (°) | Reference reduction |
| --- | ---: | ---: | ---: |
| Published, 30° | +0.00533648 | −0.0000174839 | 305× |
| Variation, 75° | −0.00458144 | −0.0000165002 | 278× |

A fresh signed OEL 1.1.0 installation produced the following results from this folder on Apple Silicon. Both cases passed all three comparisons.

| Case | Replay direct (°) | Replay mapped (°) | Replay reduction |
| --- | ---: | ---: | ---: |
| Published, 30° | +0.00533650 | −0.0000174655 | 305.5× |
| Variation, 75° | −0.00458136 | −0.0000164307 | 278.8× |

Changing the initial argument of perigee reverses the direct drift in this experiment. Initializing in mean elements still leaves a small negative residual in both cases.

In the original blinded experiment, GPT-6.2 Sol inspected the OEL simulation results. The host supplied initial-state conversion and ran the selected configurations. Sol chose simulation settings and evidence queries, measured drift and recovered the smaller residual in both cases without the paper title or reference answers. The script here replays the numerical comparison without model access. The analytical derivation, near-polar example and performance in other orbit regimes remain outside this reproduction.
