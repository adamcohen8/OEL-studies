# Cube versus spherical mascon: medial-plane study

## Run

Use Python 3.11 with NumPy, Matplotlib, and OEL 1.1.0 installed. Place the accepted, unmodified `cube_bridge.py` and `test_cube_bridge.py` beside `cube_study.py`. The comparator calls `sim.dynamics.orbit.accelerations.accel_two_body` directly; there is no replacement implementation.

```sh
python3.11 -m unittest -v test_cube_bridge
python3.11 cube_study.py --out results
```

## Geometry and conventions

The homogeneous cube has side 1 km, center at the origin, and total GM 1 km^3/s^2. Positions in the figures are expressed in side lengths. Accelerations in the tables are in km/s^2 for that normalization. Fractional errors are independent of the chosen positive GM and physical side length. The medial plane is z=0. The plotted exterior is max(|x|,|y|)>0.5; the boundary is included by continuity and the cube interior is masked. Reflection and coordinate exchange reduce the polar scan to 0 <= angle <= 45 degrees, measured from a face normal.

Equal mass does not specify the sphere radius. We use an inscribed homogeneous spherical mascon of radius 0.5 so every strict cube-exterior query is also in its exterior. The exterior field is independent of this choice of radius wherever the exterior laws of two choices both apply. The equal-volume radius is (3/(4*pi))^(1/3); `near.csv` labels points where that alternative sphere's exterior law would also apply. No interior spherical field is called, and no claim is made that equal-volume and inscribed spheres agree at near-face points inside the equal-volume sphere.

Let C and M denote cube and mascon acceleration vectors. The signed primary error is 100*(|M|-|C|)/|C|; positive means the mascon overpredicts attraction. Contours use the absolute value of this error at 40, 10, 1, and 0.1 percent. The alternatives are 100*|M-C|/|C| (which includes direction error) and 100*abs(|M|-|C|)/|M| (which changes the denominator). `near.csv` also includes the angle between C and M in degrees.

## Diagnostics and interpretation

`near.csv` samples normal offsets from a face center, other positions along a face up to its edge, and the diagonal just outside the edge. `contours.png` shows the finest reconstruction over [-2,2] in each coordinate. The contour CSVs at 200, 400, and 800 grid intervals contain vertices grouped by absolute level and segment ID, with a signed level for each segment. An absolute-error contour is the union of the positive and negative contours of the signed field. Interpolating the signed field avoids losing paired crossings when the magnitude error passes through zero within a grid cell. The printed diagnostics report segment counts, radial extents, panel and cube boundary endpoints, and direct signed-field residuals at the interpolated vertices. The plotted panel is not used to infer the farthest crossing.

`rays_90.csv` and `rays_180.csv` contain per-ray diagnostics, with the suffix denoting angular intervals in the fundamental sector. Each ray starts on the cube boundary and ends at radius 4. For each metric and level, the scan records the number of sampled below/not-below transitions and refines the last transition from not-below to below by bisection. If the endpoint is not below, the ray is flagged as unresolved. A blank crossing on an unflagged ray means it was below throughout the sampled ray. `shells.csv` gives angular envelopes on common radii. Strictly below 0.1% excludes the crossing itself. The intended summary is a radial distance from the center, not the distance from a face or the half-width of a square. Mesh refinement and outward-suffix diagnostics assess numerical consistency; they do not prove that there are no unresolved narrow recrossings between samples. The range from radius 4 outward is treated separately by the bound below.

## Analytic far-field control

Normalize side and GM to 1. Write the exterior multipole expansion of the potential kernel with coefficients A_l(n) = E[|s|^l P_l(n dot s/|s|)], where s is uniform in the cube and n is a unit field direction. Odd degrees vanish by inversion and l=2 vanishes by cubic symmetry. Using E[s_i^4]=1/80 and E[s_i^2 s_j^2]=1/144 gives A_4 = -(7/192)*(sum n_i^4 - 3/5). Thus |A_4| <= 7/480. Its spherical gradient has magnitude at most 7/96: the gradient of sum n_i^4 has norm 4*sqrt(sum n_i^6 - (sum n_i^4)^2), and the expression under the square root is the variance of n_i^2 with weights n_i^2, at most 1/4. The degree-four acceleration correction divided by the monopole magnitude is bounded by (5*(7/480)+7/96)/r^4.

For even l >= 6, use |P_l| <= 1 and |P_l'| <= l*(l+1)/2. Radial plus angular differentiation then bounds the relative acceleration correction by ((l+1)*(1+l/2))*E[|s|^l]/r^l. Here E[|s|^6]=3/448+18/(80*12)+6/1728 and |s|^2 <= 3/4. At fixed r >= 4 the sequence of these even-degree upper bounds has successive ratio at most (45/28)*(3/4)/r^2; the code sums the resulting geometric majorant. If B bounds |C-M|/|M|, then both cube-normalized metrics in fractional units are bounded by B/(1-B), and the mascon-normalized metric by B. The bounds decrease with r and apply in every direction, not only in the study plane. This control does not certify the interpolated crossing radius inside r=4.
