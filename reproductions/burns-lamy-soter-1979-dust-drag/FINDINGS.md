# Poynting–Robertson inspiral: quantitative findings

The circular secular law reproduces the simulated photosphere-crossing times to within 12 parts per million in the four targeted circular cases with β = 0.6 or 0.85. At 1 AU it underpredicts the crossing time by 0.497% at β = 0.9999 and 40.05% at β = 0.999999, with percentages relative to the secular prediction. Those discrepancies are much larger than either the tolerance-refinement shifts or the differences between the two trajectory formulations.

## Analytic prediction and grid checks

For constant β and circular motion under reduced gravity, the analytic secular law is da/dt = −2βμ/(ca). Its integrated prediction is t = [cr₀²/(4βμ)](1−q²), with q = R☉/r₀ for the photosphere. Thus the r₀²/β scaling and the cutoff correction are analytic relationships, not measured lifetime trends. Dropping the cutoff adds q² of the zero-cutoff time: 0.8651% at 0.05 AU, 0.05407% at 0.2 AU, and 0.002163% at 1 AU.

Across the 20-case grid (0.03–10 AU, β = 0.01–0.9), quadrature of the reconstructed circular drift returned a cutoff-corrected prefactor of 400.49393 yr and fitted exponents of 2.0 in r₀/AU and −1.0 in β, matching the analytic values at reported precision. The maximum absolute log residual was 1.86×10⁻¹¹. Both the maximum 4-to-8-node relative change and the maximum fine-quadrature deviation from the analytic time were 2.46×10⁻¹¹. The largest relative mean-drift error at a node was 5.88×10⁻¹¹, and the largest relative phase span was 5.06×10⁻¹⁰. This checks circular power and its scaling; the quadrature times are inferred lifetimes, not integrations of unaveraged trajectories to the photosphere.

The short unaveraged runs directly integrated at most 16 initial periods, stopping earlier when the predicted fractional decrease in osculating a² reached 0.08. The largest fine-run relative deviation in the decrease was −4.57×10⁻⁴ at 0.1 AU, β = 0.9: 0.079963426 measured versus 0.08 predicted over 0.3560 yr, with a refinement shift of 1.60×10⁻⁸ in the ratio. At the weakest signal (10 AU, β = 0.01), the directly integrated span was 508.523 yr; the decrease was 0.00012697398 versus 0.00012697396 predicted, but the ratio shifted by 3.56×10⁻⁵ under refinement, the largest shift on the grid. The 16-period β = 0 control's signed fractional a² drift fell from 4.54×10⁻⁹ to 2.43×10⁻¹¹ under refinement; its decay ratio is undefined. These runs test local drift, not accumulated lifetime error.

## Directly integrated first inward photosphere crossings

| Launch | Cartesian time (yr) | Circular secular time (yr) | (Cartesian − secular)/secular |
|---|---:|---:|---:|
| 0.05 AU, β = 0.6, circular | 1.6543041 | 1.6542890 | +9.14×10⁻⁶ |
| 0.2 AU, β = 0.6, circular | 26.685219 | 26.685159 | +2.24×10⁻⁶ |
| 0.2 AU, β = 0.85, circular | 18.836808 | 18.836583 | +1.20×10⁻⁵ |
| 1 AU, β = 0.85, circular | 471.16027 | 471.15914 | +2.41×10⁻⁶ |
| 1 AU, β = 0.9999, circular | 402.51561 | 400.52532 | +0.00496921 |
| 1 AU, β = 0.999999, circular | 560.89467 | 400.48567 | +0.400536 |
| 0.05 AU, β = 0.35, gravity-only launch speed | 8.4617223 | Not applicable | — |

At the surface, the Cartesian and angle-domain times agreed within 2.70×10⁻¹⁰ relative across the six circular launches; the largest within-solver refinement shift was 4.71×10⁻⁸. At interior cuts, the first-crossing deviations from the secular law were larger and of either sign: up to 8.51×10⁻⁴ in magnitude for the four moderate-β cases, 5.82×10⁻³ at β = 0.9999, and 0.801 at β = 0.999999. At the last β, the two formulations agreed within 3.09×10⁻¹² at the reported interior cuts, with within-solver shifts below 6.24×10⁻¹⁰. The near-unity breakdown is consistent with loss of separation between orbital and inspiral timescales: the circular-reference drag parameter at the surface rises from 0.146 to 1.46 between the two near-unity cases, even as the maximum sampled v/c falls from 1.41×10⁻⁵ to 5.00×10⁻⁷.

The gravity-only-speed launch is eccentric under reduced gravity and was not compared with the circular formula. Its angle-domain surface time was 8.4617224 yr, a relative solver difference of 1.02×10⁻⁸; the Cartesian and angle-domain refinement shifts were 1.44×10⁻⁶ and 2.23×10⁻⁷, respectively. Its first interior crossings at q = 0.7, 0.3, and 0.1 occurred at 6.1024124, 8.1888019, and 8.4578789 yr. Repeated inward crossings were reported in 12 entries, with as many as 147 detections in one run at q = 0.7. These are detection counts, not proof that all crossings were resolved. Solver agreement and refinement support the reported first-crossing times; they are empirical numerical checks, not rigorous error bounds or statistical confidence intervals.

## Scope and familiarity

The results apply to the tested initial conditions in a point-source, constant-β, first-order v/c model. The photosphere is a stopping radius, not a dust-survival prediction. The cross-formulation comparison checks the numerical solution of the same force law; it does not independently validate the physical model or establish a general β threshold for failure of circular averaging. I was already familiar with the standard Poynting–Robertson secular scaling. The numerical findings come from the supplied study records, not from a separate run or independent external check.
