# One-exposure recovery bounds after the 2024 Juice lunar flyby

**Adam Cohen · OEL original research · November 2026**

[Read the paper (PDF)](manuscript/juice-recovery.pdf) · [Computational supplement (PDF)](manuscript/juice-supplement.pdf) · [Replay the detector checks](#replay-the-detector-checks) · [LaTeX source](manuscript/)

Could a single Catalina exposure cover every direction consistent with a short optical arc after Juice passed the Moon? We fitted 25 observations from August 15–17, 2024 and evaluated an ideal, freely rotatable 4.4° square at Catalina station 703 on August 20 at 06:00 UTC.

At a 2.00″ per-observation error cap, the least-squares pointing misses one compatible trajectory by at least **4.226 arcminutes** at every roll. A **0.20557°** shift of the boresight covers the surviving family in the nominal numerical enclosure. At a cap of **2.148″**, three nominal trajectories exclude every pointing and roll. Lunar terrain removes one member of that trio. A nearby trajectory restores the exclusion in a numerical replay with rotating lunar gravity.

The choice of observing arc matters. A single C23 exposure from August 14 rules out the four trajectories used to demonstrate the two failures. The sufficient pointing at 2.00″ still covers the smaller family. Whether other trajectories demonstrate either failure on the longer arc remains open. The paper treats an ideal field and reports a numerical enclosure for the specified force model.

## Replay the detector checks

From this study folder, run the following commands in a Python environment.

```sh
python -m pip install numpy==2.4.6 mpmath==1.3.0
python reproduce.py
```

The script reads the [numerical inputs](inputs/detector.json). It checks the fixed-boresight miss and evaluates 4,096 roll intervals for each three-orbit obstruction. The nominal miss is 4.22639 arcminutes. The minimum positive chord excess is 0.0001288955 for the nominal trio and 0.0000508640 for the terrain-replacement trio under the paper’s illustrative 5″ forecast perturbations.

This compact replay starts from the reported directions and radii. The [supplement](manuscript/juice-supplement.pdf) describes the separate trajectory integration, lunar terrain checks and 638-region enclosure. The larger dynamics archive referenced in the manuscript is separate from this edition. The public numerical engine used for verification was [OEL v1.0.0](https://github.com/adamcohen8/orbital-engineering-lab/tree/108aa6412e6875bbfa5d9741ebee7200b362badc).

The [25 training observations](inputs/observations.csv) give reception times, station codes, observed directions and site vectors. The [August 14 check](inputs/previous-night.csv) gives angular residuals for archive entries 62–74, with entry 75 repeated as a timing check. [Project Pluto observation archive](https://www.projectpluto.com/pluto/mpecs/23053a.htm)
