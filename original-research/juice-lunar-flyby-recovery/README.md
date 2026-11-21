# One-exposure recovery bounds after the 2024 Juice lunar flyby

**Adam Cohen · OEL original research · November 2026**

[Read the paper (PDF)](manuscript/juice-recovery.pdf) · [Computational supplement (PDF)](manuscript/juice-supplement.pdf) · [Replay instructions](REPRODUCING.md) · [LaTeX source](manuscript/)

Could a single Catalina exposure cover every direction consistent with a short optical arc after Juice passed the Moon? We fitted 25 observations from August 15–17, 2024 and evaluated an ideal, freely rotatable 4.4° square at Catalina station 703 on August 20 at 06:00 UTC.

At a 2.00″ per-observation error cap, the least-squares pointing misses one compatible trajectory by at least **4.226 arcminutes** at every roll. A **0.20557°** shift of the boresight covers the surviving family in the nominal numerical enclosure. At a cap of **2.148″**, three nominal trajectories exclude every pointing and roll. Lunar terrain removes one member of that trio. A nearby trajectory restores the exclusion in a numerical replay with rotating lunar gravity.

The choice of observing arc matters. A single C23 exposure from August 14 rules out the four trajectories used to demonstrate the two failures. The sufficient pointing at 2.00″ still covers the smaller family. Whether other trajectories demonstrate either failure on the longer arc remains open. The paper treats an ideal field and reports a numerical enclosure for the specified force model.

## Replay the detector checks

With NumPy and mpmath installed, run `python reproduce.py` in this folder to recover the 4.22639 arcminute fixed-boresight miss and both 4,096-interval roll checks from the [published detector inputs](inputs/detector.json). The full dynamics replay below starts from observations and trajectory states.

## Explore the workspace

- [Run the geometry and exact-arithmetic checks](REPRODUCING.md#offline-checks-without-oel) using Python, NumPy and mpmath.
- [Regenerate the witness trajectories and preceding-night comparison](REPRODUCING.md#regenerate-the-four-nominal-witnesses-and-check-the-observing-arc) using the pinned public OEL v1.0.0 source and native backend.
- [Rebuild the 638-region enclosure](REPRODUCING.md#rebuild-the-nominal-enclosure) or [repeat the lunar terrain and gravity calculations](REPRODUCING.md#recompute-terrain-and-gravity).
- Inspect the [25 training observations](revision/results/observations.csv), [preceding-night residuals](review-round-2/previous-night.csv), [public-engine enclosure results](revision/results/public-cover.json) and [lunar replay results](revision/results/public-lunar.json).

The study includes the fixed-width astrometry records used by the analysis and the public gravity and terrain inputs. The ephemeris records its source kernel by filename. The published PDFs match the revised manuscript and supplement. [Project Pluto astrometry](https://www.projectpluto.com/pluto/mpecs/23053a.htm) · [Lunar data sources and checksums](revision/inputs/sources.json) · [Workspace checksums](SHA256SUMS)
