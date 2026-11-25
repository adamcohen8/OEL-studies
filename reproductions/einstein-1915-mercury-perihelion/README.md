# Mercury’s relativistic perihelion advance

Einstein (1915), *Explanation of the Perihelion Motion of Mercury from the General Theory of Relativity*, pp. 831–839. Equations (13)–(14) give the advance in an isolated solar field. [Paper and English translation](https://etienneklein.fr/wp-content/uploads/2016/01/Relativit%C3%A9-g%C3%A9n%C3%A9rale.pdf)

The measured excess of +42.98045 arcseconds per Julian century rounds to Einstein’s published prediction of +43. The calculation uses J2000 orbital elements. [JPL input elements, Table 1](https://ssd.jpl.nasa.gov/planets/approx_pos.html)

## Conclusion
The leading nonrotating relativistic correction advances Mercury’s perihelion by **+42.980 arcseconds per Julian century** relative to Newtonian gravity in this isolated Sun–Mercury model. The finest native fits give +42.98045377359536 arcsec per Julian century with the Schwarzschild correction and +7.534150136240728 × 10⁻⁷ arcsec per Julian century for the Newtonian control. Their signed difference, Schwarzschild minus Newtonian, is +42.98045302018034 arcsec per Julian century. Positive means prograde. The Newtonian control is consistent with the zero secular advance expected for an inverse-square central force.

## Model and measurement
Mercury was treated as a test particle orbiting a stationary point Sun in a Sun-relative inertial frame. Both branches started at perihelion from the same state, constructed using JPL J2000 Mercury mean elements (semimajor axis 0.38709927 AU and eccentricity 0.20563593). The JPL elements furnished initial conditions only; no ephemeris was fitted. The only force-model change was enabling the leading nonspinning solar Schwarzschild correction in the second branch. The integrations covered 20 Julian years (7,305 days or 631,152,000 s). Interior perihelia were located at negative-to-nonnegative radial-velocity crossings using native partial-step root finding. The rate was obtained by fitting unwrapped perihelion angle against elapsed time and converting the slope from radians per second to arcseconds per Julian century (3,155,760,000 s). The unrounded equivalent is +0.4298045302018034 arcsec per Julian year, a unit conversion rather than a separate fit.

## Convergence and reconstruction checks
Each branch yielded 83 perihelia, split 41 and 42 across the two 10-Julian-year halves, as expected from the Kepler period. The second-half minus first-half rates were −6.900 × 10⁻⁷ for Newtonian gravity and +1.214 × 10⁻⁶ for the corrected branch, in arcsec per Julian century. The full-span angular-fit RMS residuals were 2.33 × 10⁻⁸ and 3.38 × 10⁻⁸ arcsec, respectively. Halving the maximum outer propagation-call span down to 1,855.456 s gave final two signed changes in the branch difference of +7.928 × 10⁻⁷ and −2.329 × 10⁻⁶ arcsec per Julian century. At that spacing, tightening the relative tolerance from 10⁻¹¹ to 10⁻¹² and then 10⁻¹³ left the reported difference unchanged. As a separate event-reconstruction check, cubic Hermite roots on sampled outputs from the same finest trajectories gave +42.98037078231322 arcsec per Julian century at the finest output spacing. The last output-spacing halving shifted that estimate by +1.4543493156082832 × 10⁻⁴, and the finest sampled estimate was −8.223786711880621 × 10⁻⁵ below the native-root result, in the same units. These diagnostics pass the 0.001 arcsec per Julian century criterion. The largest shift among the final checks, about 1.45 × 10⁻⁴ arcsec per Julian century, is an observed numerical scale, not an error bound.

## Analytic comparison and limits
The leading weak-field formula, 6πGM☉/[a(1−e²)c²] per orbit, gives 5.01866043967899 × 10⁻⁷ rad per orbit. Using the Newtonian Kepler period gives +42.980475409730644 arcsec per Julian century. The native numerical difference minus that analytic estimate is **−2.238955030264833 × 10⁻⁵ arcsec per Julian century**. This reproduces the familiar roughly 43-arcsecond-per-century Mercury benchmark. The formula is a leading-order comparison using Newtonian elements and period, not an exact solution of the integrated truncated force law. Planetary perturbations, solar oblateness, solar spin and frame dragging, finite Mercury mass and solar recoil, higher-order relativistic terms, and reference-frame precession are excluded. The century unit rescales a slope measured over 20 years; no century-long evolution or long-term constancy was tested. Independent replay on the public engine reproduced the complete numerical report. A separate fit to the 83 saved perihelion times and angles in each branch recovered the same difference. The physical comparison concerns the solar monopole contribution to precession.

## Research roles

GPT-6.2 Sol designed the comparison, authored the driver and chose the numerical checks. It reported prior familiarity with Mercury’s benchmark. After seeing an epoch-coverage validation error it chose a 20-year baseline. When sampled perihelia failed its original convergence criterion, it refined the grid and retained that criterion. The supervisor implemented the OEL force selection, supplied the exact setting and solar GM, ran the driver, and independently checked the results.

## Replay

The tested public build uses CPython 3.11, Rust/Cargo and an editable OEL installation. From the downloaded study folder, the following POSIX commands create the environment and build both public OEL extensions. The revision lock records both the development commit used for the experiment and its public export.

```sh
study_dir="$PWD"
python3.11 -m venv .venv
. .venv/bin/activate
python -m pip install "maturin>=1.10,<2"
git clone https://github.com/adamcohen8/orbital-engineering-lab.git oel-source
git -C oel-source checkout --detach e3b0e7620d4ae616d01da10b7ee23ff282d17800
cd oel-source
python -m maturin build --release --locked --manifest-path rust/oel-orbit/Cargo.toml --interpreter python --out "$study_dir/.wheels"
python -m maturin build --release --locked --manifest-path rust/oel-game/Cargo.toml --interpreter python --out "$study_dir/.wheels"
python -m pip install -c constraints/py311.txt "$study_dir"/.wheels/*.whl -e . -r "$study_dir/requirements.lock"
cd "$study_dir"
python replay.py
```

The replay checks the OEL Git revision, tracked files, package version and native extension capability before starting. It writes `outputs/mercury.json` with the input and normalized configurations for each run, perihelion times and angles, both fitted slopes, and convergence diagnostics. The committed `results.json` retains the measured rates and checks for comparison. A Julian century here is 36,525 days.
