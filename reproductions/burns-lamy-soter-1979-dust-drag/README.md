# Inspiral of a dust grain under Poynting–Robertson drag

Burns, Lamy and Soter (1979) give a lifetime of approximately **400 R²/β years** for a grain on a slowly evolving circular orbit. Here R is the initial heliocentric radius in astronomical units and β is the ratio of direct radiation pressure to solar gravity. This workspace reproduces the normalization and both exponents in Eq. (50), then checks how accurately the approximation predicts a grain's first inward crossing of the solar photosphere. [Original paper, p. 31](https://www-users.cse.umn.edu/~gehrz001/ASTRO_2001_Spring_Semester_2020/10_Burns_Lamy_Soter.pdf#page=31)

## Evidence

Across 20 circular-orbit cases spanning 0.03–10 AU and β of 0.01–0.9, a fit to the force-derived, cutoff-corrected lifetimes gave **400.494 years**, a radius exponent of **2.000**, and a β exponent of **−1.000**. The coefficient differs from the paper's rounded value by **+0.123%**. This fit uses circular orbital power sampled from the force model. Independent unaveraged trajectories test whether the circular approximation describes actual motion. [Measurements](results.json) · [Numerical method](METHOD.md) · [Paper comparison](comparison.json)

| Initially circular launch at 1 AU | Circular estimate to photosphere | Cartesian trajectory to photosphere |
| --- | ---: | ---: |
| β = 0.85 | 471.159 years | 471.160 years |
| β = 0.9999 | 400.525 years | 402.516 years |
| β = 0.999999 | 400.486 years | 560.895 years |

Equation (50) integrates to zero radius. For the photosphere comparison, the calculation stops at 695,700 km and uses the same endpoint in the circular estimate. As β approaches one, radiation pressure nearly cancels solar gravity and the circular averaging assumption loses accuracy. These are idealized dynamical times. Grain destruction, solar wind, collisions and planetary forces are outside the study. [Analyst's findings](FINDINGS.md) · [Trajectory checks](summary.json)

## Reproduce

Start in this folder on Apple Silicon macOS with CPython 3.11. The public release bundle pins the complete offline installation, and the lock lists the direct dependency hashes. The force bridge and all study inputs are included. [OEL v1.1.0 release](https://github.com/adamcohen8/orbital-engineering-lab/releases/tag/v1.1.0)

```sh
curl -fLO https://github.com/adamcohen8/orbital-engineering-lab/releases/download/v1.1.0/oel-public-1.1.0-arm64-py311.bundle.zip
echo '355a09916cfccb4ef3eade15e799040e785c7ea57b1e8612bce70784f3a1d20e  oel-public-1.1.0-arm64-py311.bundle.zip' | shasum -a 256 -c -
unzip -q oel-public-1.1.0-arm64-py311.bundle.zip 'wheelhouse/*' -d bundle
python3.11 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links bundle/wheelhouse orbital-engineering-lab==1.1.0
.venv/bin/python -m unittest -v test_dust_bridge
.venv/bin/python verify_bridge.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python study.py
.venv/bin/python summary.py
.venv/bin/python compare_paper.py
```

## Contributions

GPT-6.2 Sol authored the study-specific radiation force, its tests, both trajectory formulations, experiment design, measurement summaries and scientific interpretation. Sol reported prior familiarity with the approximate lifetime formula before seeing the numerical results. OEL supplies the central gravity acceleration. The study bridge applies direct radiation pressure and velocity-dependent drag, and SciPy integrates the trajectories. The supervisor set the force interface, verified it before the study, withheld the paper and grading criteria, relayed execution, and checked the final comparison and replay. [Bridge](dust_bridge.py) · [Independent force checks](verify_bridge.py) · [Experiment](study.py)
