# Exterior attraction of a cube and an equal-mass spherical mascon

## Setup

The study uses a homogeneous, axis-aligned cube of side 1 km and total GM 1 km³/s², centered at the origin. Distances below are in cube side lengths. In the medial plane `z = 0`, we evaluate the cube exterior and include exact boundary values as exterior limits. The cube field comes from the validated `cube_bridge.py`; the spherical exterior field comes from the installed OEL 1.1.0 `accel_two_body` routine. The revised replay and all five bridge self-tests passed, as did the study controls for inverse-square values, rotational symmetry, and inward projection.

Equal mass does not determine the sphere's radius. We chose the inscribed radius of 0.5 so the exterior law applies throughout the cube exterior. The same law applies to an equal-volume sphere only at and beyond its radius, 0.62035049; the near-face values inside that radius are not predictions for the equal-volume sphere's interior. The signed error is `100 × (|M| − |C|)/|C|`, where `M` and `C` are the mascon and cube acceleration vectors. Positive means the mascon overpredicts attraction. The map and primary threshold use the absolute value of this signed error. I did not recognize a specific paper from the description and used no external paper values.

## Near the face

At the face center `(0.5, 0, 0)`, the exterior-limit attraction magnitudes are 2.596897 km/s² for the cube and 4.000000 km/s² for the mascon, an overprediction of **54.03%**. Along the outward normal, the overprediction falls to 41.54% at radius 0.55, 32.33% at 0.60, and 6.045% at 1.00. This is not representative of the whole face: at `(0.5, 0.4, 0)` the mascon overpredicts magnitude by only 1.515%, but the vectors differ in direction by 8.95° and the relative vector error is 15.79%. At the medial edge `(0.5, 0.5, 0)` the mascon **underpredicts** attraction by 8.86%; just outside the edge on the diagonal, at `(0.501, 0.501, 0)`, it underpredicts by 8.35%. See `results/near.csv`.

## Where the error stays below 0.1%

The largest measured last downcrossing was **2.9092 side lengths from the center**, along a face normal, at both scan resolutions. At 15°, 30°, and 45° from that normal the fine-scan last crossings were 2.6543, 1.6835, and 1.9379, respectively. Nine of the 181 fine-scan rays also crossed back up through 0.1% at smaller radii, so a first downcrossing is not a reliable cutoff. No sampled ray was unresolved at radius 4. The sampled angular maxima were 0.101374% at radius 2.8992, 0.098649% at 2.9192, 0.088530% at 3, and 0.028222% at 4. Using vector error instead of magnitude error gave the same largest sampled crossing; switching the magnitude denominator to the mascon moved it to 2.9085.

About **2.91 side lengths** is therefore the numerical estimate for the outer 0.1% threshold in this plane. The crossing itself is at equality. The ray and shell checks do not rule out every possible recrossing between samples below radius 4, nor do they establish the same cutoff outside the plane. Separately, the multipole bound derived in `METHOD.md` bounds the cube-normalized vector error, and hence the magnitude error, by 0.078407% for all radii at least 4, in every direction.

## Reconstructed panel and reproducibility

`results/contours.png` covers `[-2, 2]` on each axis and masks the cube interior. Each absolute level combines the positive and negative contours of the signed field. Tracing the signs separately avoids losing paired crossings when the signed error passes through zero within a grid cell. The host's visual check confirmed that the external legend is unclipped, the four levels are distinguishable, and the square mask and axes are intact.

| Absolute error | Segments in the 800-interval panel (+ / −) | Farthest plotted vertex radius | Largest vertex residual (percentage points) |
|---:|:---:|---:|---:|
| 40% | 4 / 0 | 0.55742 | 0.00383 |
| 10% | 4 / 0 | 0.86626 | 0.000716 |
| 1% | 4 / 4 | 1.61748 | 0.00257 |
| 0.1% | 8 / 4 | 2.19227 | 0.00237 |

The signed segment counts were unchanged across 200, 400, and 800 grid intervals. The largest 0.1% vertex residual fell from 0.03304 to 0.00945 to 0.00237 percentage points. The positive 0.1% contour meets the panel boundary at eight endpoints. Its largest plotted radius is therefore not the outer threshold; the ray scans extend beyond the panel for that measurement.

Use the commands in `README.md` to regenerate the figure and tables. The three contour CSVs preserve signed levels and vertex coordinates; the two ray CSVs and `results/shells.csv` preserve the numerical threshold checks.
