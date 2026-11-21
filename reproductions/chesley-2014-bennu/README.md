# Bennu’s orbital drift from optical and radar astrometry

Chesley et al. (2014) measured a small inward drift in the semimajor axis of asteroid (101955) Bennu. This workspace independently fits the 1999–2013 observation arc and checks the signed drift inferred from the fitted transverse acceleration. It reproduces the published **direction and magnitude** using a different optical selection and force model. The authors’ manual solution 87, the bulk density estimate and future encounter predictions are outside this replay. [Original paper, sections 2–3 and Tables 1, 3 and 5](https://echo.jpl.nasa.gov/asteroids/chesley_etal_2014_Bennu.pdf)

| Quantity | Published | Workspace replay |
| --- | ---: | ---: |
| Semimajor-axis drift (metres per Julian year) | −284 ± 1.5 | −282.90 ± 1.48 |
| Transverse acceleration at 1 au (10⁻¹⁴ au/day²) | −4.618 ± 0.024 | −4.603 ± 0.024 |
| Optical positions included | 478 | 561 |
| Radar delays / Doppler measurements | 22 / 7 | 22 / 7 |

The estimates differ by 1.1 m/year. The fits share observations, so that difference should not be interpreted as a statistical test between independent measurements. The uncertainties are formal one-standard-deviation estimates under each fit’s assumptions.

## Replay

Use an Apple Silicon Mac with CPython 3.11, Git, curl and Apple’s Command Line Tools. Install the signed OEL 1.1.0 core release and put its launcher on PATH. The workspace manifest selects that release. [Public OEL release](https://github.com/adamcohen8/orbital-engineering-lab/releases/tag/v1.1.0)

From this folder, run the commands below. Setup downloads the pinned GRSS source, build dependency and about 218 MiB of planetary and orientation kernels. It rebuilds GRSS’s Python extension from source because the tested macOS wheel alone failed to resolve its integration routine. Use a dedicated virtual environment because setup installs study-specific station data and kernels there.

```sh
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install --only-binary=:all: -r requirements.lock
python setup_bridge.py
python reproduce.py
```

Use `python reproduce.py --oel /absolute/path/to/oel` for a launcher outside PATH. Subsequent replays use local inputs and require no model credentials. Each replay writes both fitted states, covariance matrices, per-observation residuals and a `comparison.json` under a new `outputs/replay-*` directory. A failed comparison produces a nonzero exit code.

## Which observations went into the fit?

We matched 569 optical positions from 43 stations in the paper’s date range by station, UTC time and sky coordinates. Chesley et al. manually removed 91 positions, sometimes rejecting a whole night from one station. We could not identify the authors’ complete row-level selection and final weights in the accessible sources. This replay uses ESA’s published selection flags, angular weights and catalog biases for those observations. Eight positions are excluded. All 22 center-of-mass round-trip delays and seven Doppler observations are included with the reported radar uncertainties. [Minor Planet Center observations](https://data.minorplanetcenter.net/api/get-obs) · [ESA NEO Coordination Centre observation data](https://neo.ssa.esa.int/PSDB-portlet/download?file=101955.rwo) · [JPL radar astrometry](https://ssd-api.jpl.nasa.gov/sb_radar.api?des=101955)

The frozen numerical inputs include the MPC observation identifiers, UTC epochs, angular positions in degrees, and tangent-plane biases and standard deviations in arcseconds. The radar records retain their original units, transmitter, receiver, carrier frequency and center-of-mass flag. Station longitudes and parallax constants come from the MPC. The source endpoints and input digests are recorded in `inputs/provenance.json`. [MPC observatory codes](https://data.minorplanetcenter.net/api/obscodes)

## What GRSS, OEL and the analyst each did

**GRSS 4.6.0** fits six heliocentric orbital elements to the optical and radar data. A second fit also estimates a transverse acceleration proportional to heliocentric distance raised to −2.25. Both fits use the same observations and weights and start from a rounded reference geometry. The acceleration starts at zero and has no fitted-parameter prior. GRSS uses DE440 and 16 asteroid perturbers here. The published reference used DE424 and 25. For the blinded analyst comparison, the direction had to agree and the inferred drift had to lie within 4.5 m/year of the published value. We set this criterion after the bridge pilot and before contacting the analyst. The gravity-only radar delay residuals have an RMS of 40.3 µs, falling to 0.605 µs with the transverse term. [GRSS source at the tested revision](https://github.com/rahil-makadia/grss/tree/b0bc1124795493ad88589af2988ff0502cacc188)

**OEL 1.1.0** converts the fitted osculating elements into Cartesian positions and velocities with the solar gravitational parameter. The replay computes the work done by the transverse acceleration, averages the resulting semimajor-axis rate over orbital time, and propagates the fitted covariance into the drift uncertainty. It compares this with a separate eccentric-anomaly average. Agreement tests the conversion from acceleration to drift. GRSS performs the underlying orbit determination.

**GPT-6.2 Sol** requested the fit evidence and received the acceleration, uncertainty, orbital elements and goodness of fit for each model. The paper, object name, reference drift and replay’s drift calculation were withheld. Sol derived the orbit average and reported −282.9 ± 1.48 m/year. It correctly identified the inward drift and qualified the uncertainty as conditional on the fit. The host selected and reduced the observations and ran both fits. This tests interpretation of fitted evidence, including its sign and units.
