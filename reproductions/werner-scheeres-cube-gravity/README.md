# Gravity near a uniform cube

Werner and Scheeres compare a homogeneous cube with a spherical mass of equal mass in Figure 15 of *Exterior Gravitation of a Polyhedron Derived and Compared with Harmonic and Mascon Gravitation Representations of Asteroid 4769 Castalia*. This standalone OEL workspace tests the cube example in the central plane. The geometry is specified by the paper and requires no asteroid shape file. [Paper, Figure 15 on p. 339](https://doi.org/10.1007/BF00053511) · [Accessible copy](https://adsabs.harvard.edu/pdf/1997cemda..65..313w)

**The near-face discrepancy and distance scale were recovered. The blinded analyst’s primary map only partially reproduces the published contours.** In the analyst’s calculation, the spherical approximation gives **54.03% stronger attraction at the face center**. The authors describe roughly 50%. The last sampled 0.1% crossing is **2.9092 cube side lengths** from the center along the face normal. This satisfies the paper’s statement that more than 2.5 units are needed. It exceeds the reviewer’s preselected figure-scale comparison interval of 2.25–2.85. The paper gives a lower bound and the displayed axes stop at ±2.

## Which error is plotted?

Sol used the absolute difference of the two acceleration magnitudes, divided by the cube attraction. It also computed the magnitude of the **vector difference** along rays as a sensitivity check. The quantities agree on a face normal. Away from that axis, the vector difference includes the change in direction.

After the blind conclusion was frozen, the reviewer plotted Sol’s alternate vector metric and sampled the lower positive-y crossing of the right-hand 1% lobe. Its ordinates agree with the printed figure to about 0.004 side lengths at the four sampled locations. This is evidence for the vector interpretation. The four readings were selected during review and do not constitute a pointwise verification of the entire panel.

| x / side | Paper y / side, approximate | Difference of magnitudes | Magnitude of vector difference |
| ---: | ---: | ---: | ---: |
| 0.8 | 0.684 | 0.508 | 0.682 |
| 1.0 | 0.695 | 0.530 | 0.693 |
| 1.2 | 0.640 | 0.504 | 0.637 |
| 1.4 | 0.501 | 0.410 | 0.501 |

For the digitization, the reviewer rendered PDF page 27 (printed page 339) at 2,400 pixels high. The bounding axes of Figure 15 fell near columns 534 and 1,087 and rows 212 and 760.5. The crossing estimates use linear interpolation between the −2 and +2 axis limits. Allow roughly 0.01 side lengths for reading the printed lines. [Blind magnitude-error map](results/contours.png) · [Reviewer’s vector-error map](results/vector_comparison.png) · [Numerical comparison](results/figure_check.json)

## Authorship and checks

GPT-6.2 Sol authored the cube bridge, its tests, the study script, grid and ray refinements, and the blinded findings. It modeled the sphere as an inscribed equal-mass sphere, so its exterior acceleration is available everywhere on and outside the cube. OEL 1.1.0 supplied that spherical acceleration. Sol’s bridge computed the homogeneous cube acceleration from analytic face integrals. [Analyst’s method](METHOD.md) · [Analyst’s findings](FINDINGS.md)

The reviewer tested the bridge before the study could proceed. All five author tests passed. An independent volume-integral calculation agreed for all 38 vectors tested, including interior points, exact faces, edges, corners, near-boundary points and rescaled inputs. The portable verifier checks convergence at two quadrature orders. Its maximum absolute vector discrepancy was below 7 × 10⁻¹⁵ km/s² on the pinned installation. These checks cover the bridge’s documented domain and sampled conditions. The grader, source comparison and packaging were reviewer contributions.

## Install the pinned runtime

From the study folder, on the tested Apple Silicon macOS and CPython 3.11 configuration, run the following commands. The bundle digest pins the complete wheelhouse. The lock records the direct wheel hashes. The analytical cube bridge ships in this folder. [OEL 1.1.0 public release](https://github.com/adamcohen8/orbital-engineering-lab/releases/tag/v1.1.0)

```sh
curl -fLO https://github.com/adamcohen8/orbital-engineering-lab/releases/download/v1.1.0/oel-public-1.1.0-arm64-py311.bundle.zip
echo '355a09916cfccb4ef3eade15e799040e785c7ea57b1e8612bce70784f3a1d20e  oel-public-1.1.0-arm64-py311.bundle.zip' | shasum -a 256 -c -
unzip -q oel-public-1.1.0-arm64-py311.bundle.zip 'wheelhouse/*' -d bundle
python3.11 -m venv .venv
.venv/bin/python -m pip install --no-index --find-links bundle/wheelhouse orbital-engineering-lab==1.1.0
.venv/bin/python -m unittest -v test_cube_bridge
.venv/bin/python verify_bridge.py
```

From the same study folder, regenerate the analyst’s results and the reviewer’s separate comparison:

```sh
.venv/bin/python cube_study.py --out results
.venv/bin/python compare_figure.py
```

The programs write the figures and numeric tables to `results/`. The study and verification run offline once the public bundle is installed.
