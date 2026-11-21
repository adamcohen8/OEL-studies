# Reproducing the Juice one-exposure results

Read the [paper](manuscript/juice-recovery.pdf) and [computation supplement](manuscript/juice-supplement.pdf). Section 4 of the supplement follows the limiting box from its six incoming coordinate ranges through the four final camera-face bounds. The recorded observations, complete partition, LPs, original trajectories, alternate trajectory and NASA inputs are included. No external service is needed during numerical replay once the environment is installed.

## Offline checks without OEL

From this study folder, using Python with NumPy and mpmath installed, run

```sh
shasum -a 256 -c SHA256SUMS
python revision/check_geometry.py
python revision/check_cell_duals.py
```

The first script checks the least-squares pointing and the nominal and modified-model roll inequalities from serialised directions. It evaluates 4,096 roll intervals with mpmath at 36 decimal digits. It also verifies which original trajectory fails terrain, which two original trajectories fit together, and the three-orbit replacement. The second evaluates all four stored detector-support duals using exact fractions, including the dual residual in each bounded variable. Expected limiting face total is −4.826692 km. These quick checks use the archived outputs of the nonlinear dynamics.

## Environment for regenerating those outputs

The original computation used macOS arm64, CPython 3.14.10, NumPy 2.4.6, SciPy 1.17.1, Numba 0.66.0, llvmlite 0.48.0 and mpmath 1.3.0. The terrain replay adds SpiceyPy 8.4.0. Matplotlib 3.10.9 is needed only to redraw Figure 1. The original local OEL backend was oel-rust-orbit 0.12.0. The public replay uses the tagged engine and backend below.

The corresponding OEL v1.0.0 checkout is available at commit `108aa6412e6875bbfa5d9741ebee7200b362badc`. The public build uses oel-rust-orbit 0.18.4. Follow the public release's source installation instructions for the native toolchain. The following commands specify the orbit module and numerical packages exercised by the reviewer scripts. The commands assume an activated Python 3.14 environment. Run them from this study folder in your clone. Keep the study and the OEL source checkout in separate directories.

```sh
git clone --branch v1.0.0 --depth 1 https://github.com/adamcohen8/orbital-engineering-lab.git ../juice-oel
git -C ../juice-oel rev-parse HEAD  # 108aa6412e6875bbfa5d9741ebee7200b362badc
export JUICE_OEL="$(cd ../juice-oel && pwd)"
python -m pip install 'maturin>=1.10,<2' numpy==2.4.6 scipy==1.17.1 numba==0.66.0 llvmlite==0.48.0 mpmath==1.3.0 spiceypy==8.4.0 PyYAML
python -m maturin build --release --locked --manifest-path "$JUICE_OEL/rust/oel-orbit/Cargo.toml" --interpreter "$(command -v python)" --out ../juice-wheels
python -m pip install --no-index --find-links ../juice-wheels oel-rust-orbit==0.18.4
export PYTHONPATH="$JUICE_OEL"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
```

## Regenerate the four nominal witnesses and check the observing arc

With the environment above active, run each command from this study folder.

```sh
python review-round-2/recheck_witnesses.py
python review-round-2/previous_night.py
python review-round-2/refit_pointing.py
```

The witness command starts from the four stored six-component states and the published least-squares boresight. Every invocation computes fresh trajectory enclosures, training residual bounds, lunar clearance floors and forecast directions and radii. It tests the fixed-boresight miss and all 4,096 roll intervals for A, B and C. On the original stack the miss is 4.22639 arcminutes and the minimum chord excess is 0.0001288955. The four regenerated certificates and both detector checks are recorded in `review-round-2/rechecked-witnesses.json`.

The preceding-night command writes great-circle residuals in arcseconds to stdout as CSV. Rows 62–74 are the 13 August 14 observations. Row 75 is the first training exposure and serves as an epoch check. The last C23 observation (row 74) gives 1.59138″ for the fitted orbit, 2.06116″ for the missed orbit and 3.35181″, 2.34393″ and 2.83485″ for nominal A, B and C. The archived comparison is `review-round-2/previous-night.csv`. The script includes the additional station and Earth-orientation constants, with their sources identified in the supplement.

The refit command writes three local fits as JSON. It uses equal-weight tangent coordinates, iterated emission states and centered quadratic interpolation of Earth's displacement. On the original stack, changing the finite-difference step scale gives boresight shifts of 0.0067–0.0071″ at a radial tangent-plane RMS of 0.269058″. See `review-round-2/refit-pointing.json` for the original numerical values. Small changes across OEL builds and numerical differentiation steps are expected.

## Rebuild the nominal enclosure

```sh
python revision/replay_cover.py ../juice-rerun/study --oel "$JUICE_OEL"
python revision/worked_cell.py ../juice-rerun/worked-linear.npz > ../juice-rerun/worked-cell.json
```

The replay requires a new destination. It rebuilds the 175-variable incoming optical program and both generations of reference-relative constraints from the 25 observations. It then reconstructs the 606 cuts of the proposed tree, recomputes its contraction bounds, widens the boxes, performs fresh nonlinear classifications of all 638 terminal regions, and checks the final placement and 14 witness memberships. The supplied chart, tree topology, camera and witness states are inputs. There is no need to repeat the heuristic search for them to test their validity. To isolate the box and propagation stages, add `--frozen-optics` to the first command. Failures and unresolved boxes are scientific outputs. Inspect individual rows in the two derived shard files if the result differs.

On the original stack a fresh run from the observations gave 394 contained regions, 176 collision exclusions, 68 infeasible regions and zero unresolved. Its least favorable face was −4.826691853351395 km. A separate rebuild from the observations with public OEL v1.0.0 and its source-built native backend returned the same four counts and −4.826459 km on the limiting face. The four high-degree lunar replays agreed in forecast direction to 0.000014 arcsecond. See `revision/results/public-cover.json` and `public-lunar.json`. Both rebuilds share the research enclosure algorithms. An independent derivation would also have to justify construction of the floating-point inequalities and the analytic and arithmetic allowances in the nonlinear propagation.

## Recompute terrain and gravity

```sh
python revision/replay_physical.py ../juice-lunar-replay
```

This propagates A, B, C and the replacement from their initial states with point-mass gravity, truncated GRGM1200A and the near-encounter 1200-degree field. It repeats the higher-degree run with a second integrator. It computes the terrain clearances, all training residuals and the forecast directions, then checks the geometry using the newly computed directions. Expected signs on the high-resolution grid are negative for A and positive for C and the replacement. B is thousands of kilometres above the lunar terrain. B has the smallest margin to the training cap, 0.004907″ at entry 79. The 5″ perturbations test only the forecast directions in the three-orbit geometry calculation. The high-degree propagations have empirical numerical comparisons rather than validated continuous-error radii.

For one trajectory run `python revision/lunar_replay.py A 1200 10`. The arguments are witness, harmonic degree and encounter knot spacing in seconds. Zero selects the original point-mass model. A value such as `A+.003` moves the initial state that fraction of the way from A to the archived least-squares state. An optional fourth argument selects the OEL integrator for the harmonic case. The full coefficient table, both product labels, full 16-pixel/degree terrain map, 128-pixel/degree regional crop, lunar frame kernels and sliced DE430 libration data are in `revision/inputs/`. The range and digest of the source image band are in `terrain128-source.json`. The input preparation scripts allow regeneration from the NASA source products using the URLs and checksums in `sources.json`.

## Evidence map

| File | Contents |
|---|---|
| `revision/results/observations.csv` | All 25 reception epochs, site codes, site positions and observed unit rays |
| `followup/corrected/missions/` | Archived observations, least-squares state, reference directions, ephemeris comparison inputs |
| `followup/corrected/de440-geocentric.npz` | Corrected compact translational ephemeris |
| `field-certificate/nominal-LP.npz` | Base optical inequalities and physical coordinate bounds |
| `square-decision/contracted.npz` | Tightened incoming optical and dynamics constraints |
| `square-decision/chart.json`, `partition-routes.json`, `audited-boxes.json` | Chart, full split schedule and expanded leaf domains |
| `square-decision/final-shard-0.json`, `final-shard-1.json` | Frozen classification per terminal region |
| `revision/results/rebuilt-cover.json` | Result regenerated from the observations with the original OEL environment |
| `revision/results/worked-cell.json`, `worked-linear.npz` | Limiting route, incoming box, budgets, sparse LP and four sets of dual multipliers |
| `revision/results/geometry.json` | Least-squares miss, boresight displacement, original and replacement roll checks |
| `revision/results/*1200*.json` | Initial state, 25 residuals, sky direction, terrain minimum and integration settings for each orbit |
| `revision/decision_figure.py` | Published Figure 1, generated from the archived directions |

The paper, computation supplement and `revision/results/` contain the final results. The numerical search products record intermediate steps.
