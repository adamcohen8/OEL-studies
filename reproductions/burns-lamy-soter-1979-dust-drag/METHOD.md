# Poynting–Robertson inspiral study

Run `python study.py` with NumPy, SciPy, and the supplied `dust_bridge` importable. The script writes `results.json` in the working directory. It uses inexpensive quadrature and short trajectory windows for the broad grid, reserving full integrations for seven targeted cases.

## Prediction

For constant radiation-to-gravity ratio β, set m = μ(1−β) and k = βμ/c. On a circular orbit the drag gives specific power −km/r³. Combining this with E = −m/(2a) gives da/dt = −2k/a. Integrating to a cutoff q times the initial radius gives t = c r₀²(1−q²)/(4βμ). Use q = R☉/r₀ for the photosphere; dropping the cutoff overestimates the surface time by q² of the zero-cutoff prediction.

## Tests and interpretation

- **Circular drift:** At 0.03, 0.1, 1, 5, and 10 AU and β = 0.01, 0.1, 0.5, and 0.9, reconstruct circular drift from bridge power at several phases in a tilted plane. Integrate its inverse using 4 and 8 Gauss–Legendre nodes and compare to the prediction. Fit the prefactor and exponents after dividing out the cutoff factor. This checks the bridge's circular power and scaling, not the accuracy of circular averaging on a real trajectory.
- **Local unaveraged drift:** Evolve circular launches for at most 16 initial orbital periods, ending earlier if the secular law predicts an 8% decrease in a². Compare the measured change in osculating a² to the prediction at two tolerances. The β = 0 run measures numerical drift on the same diagnostic; its decay ratio is undefined and recorded as null.
- **Full trajectories:** Compare first inward radius crossings from Cartesian DOP853 using the bridge and from the independently coded angle-domain RK45 model, each at two tolerances. The latter uses the exact scaled torque identity j = j₀−λθ for this force law. Examine changes under refinement and agreement between formulations before interpreting deviations from the secular law. Multiple crossings are possible; stored counts are detections, not a guarantee that every crossing was resolved. The gravity-only-speed release is initially eccentric under reduced gravity and is not compared to the circular lifetime formula.

The model assumes a point-source radiation field, constant β, and first-order drag in v/c; the photosphere is a stopping surface, not a dust-survival model. A small sampled v/c does not ensure slow inspiral relative to orbital motion, especially near β = 1 where effective gravity is weak. Use the unaveraged comparisons to assess that regime; numerical agreement tests the equations, not the physical assumptions.
