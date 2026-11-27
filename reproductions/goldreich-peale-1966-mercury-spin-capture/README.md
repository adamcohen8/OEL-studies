# Mercury’s capture into 3:2 spin–orbit resonance

The averaged single-passage estimate for 3:2 spin–orbit capture was 17.83%. The original switch-band criterion classified one trajectory as captured and left 17 unresolved. Extending the certified band classified those 17 as escapes and retained the same finite-time capture. At both numerical resolutions, the two periapsis launches yielded 0 of 6 captures each, and the apoapsis launch yielded 1 of 6. These fractions describe the tested orientations and neither validate nor refute the averaged probability, which assumes a different initial-phase distribution.

## Comparison with the paper

Goldreich and Peale studied a single approach to the 3:2 spin–orbit resonance. At eccentricity 0.2 their constant-angular-lag case used a triaxiality of 0.0001 and a dimensionless tidal coefficient of 0.00001. They quote a capture probability of 0.16 and report one capture among six numerical trials. Figures 10 and 11 show one capture and one passage. [Original paper, pp. 431–432](https://adsabs.harvard.edu/pdf/1966AJ.....71..425G)

Our averaged calculation gives **0.1783**, or **1.83 percentage points above** the value quoted in the text. The numerical trials tested six equally spaced initial orientations from 0 to π/2. The paper describes the initial spin only as slightly above resonance. The following three launch choices were fixed before inspecting outcomes. The first was the primary experiment.

| Starting orbital position | Initial spin / mean motion | Captures | Escapes |
| --- | ---: | ---: | ---: |
| Perihelion | 1.525 | 0 | 6 |
| Perihelion | 1.522 | 0 | 6 |
| Aphelion | 1.525 | 1 | 5 |

In the original run, 17 trials reached the lower limit of the range used to check tidal torque reversals. Those 17 are recorded as unresolved in `study.json`. The follow-up extended that lower limit from 1.468 to 1.455 times mean motion. It retained the same launches, time steps and classification criteria and produced the table above. Every follow-up outcome and stopping orbit agreed at 64 and 128 steps per orbit. The captured aphelion trial began at relative orientation π/5 and met the 256-orbit libration criterion at orbit 894.

The historical count falls within the range obtained here. The primary launch gave a different count, and the initial conditions in the paper are insufficient to identify its individual trajectories. The averaged estimate assumes a distribution of incident energies. The six deterministic orientations at any one launch do not measure that distribution. This comparison addresses the prescribed historical tidal law and a single passage through resonance.

## Replay

This workspace uses Python 3.12 and the public OEL preview at commit `76ad5422bb476dd6cd3e162fc6346b3f7af6c069`. The numerical model uses the OEL orbit and gravity-gradient functions with the constant-angular-lag tide in `spin_bridge.py`. The run uses a fixed Kepler ellipse. The choice of semimajor axis sets the conversion between dimensionless time and seconds.

With Git and a Rust toolchain installed, start in the downloaded study folder. The following commands build the public engine and run the numerical checks.

```sh
study_dir="$PWD"
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install 'maturin>=1.10,<2'
git clone --single-branch --branch preview/hyperion-passive-attitude https://github.com/adamcohen8/orbital-engineering-lab.git oel-source
git -C oel-source checkout --detach 76ad5422bb476dd6cd3e162fc6346b3f7af6c069
for component in oel-orbit oel-game; do
    python -m maturin build --release --locked --manifest-path "oel-source/rust/$component/Cargo.toml" --interpreter python --out "$study_dir/.wheels" || exit 1
done
python -m pip install -c oel-source/constraints/py312.txt "$study_dir"/.wheels/*.whl -e ./oel-source -r requirements.lock
python verify_engine.py
python check_bridge.py
python study.py --output replay.json
python followup.py --output followup-replay.json
```

`check_bridge.py` independently tests the torque sign, scale and zero at equal angular rates, the conservative torque, and a full orbit against Equation 34 in true anomaly. `study.json` holds the reference numerical output. The original orbit and spin are planar. The capture classifications apply over the integration and observation windows specified in that file.

The follow-up changes only the lower bound used to audit torque reversals. Its reference output is `followup.json`. The bound on crossing slope comes from the analytic orbit. Agreement with the interpolated orbit was tested numerically.
