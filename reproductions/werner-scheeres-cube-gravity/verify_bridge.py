"""Independent numerical acceptance check for the homogeneous cube bridge."""
import json
import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from cube_bridge import cube_accel


def reference(point, side, mu, order):
    """Integrate Newton's volume kernel via its separable Gaussian transform."""
    # 1/r^3 = (4/sqrt(pi)) * integral_0^infty t^2 exp(-t^2 r^2) dt.
    # Integrating each physical coordinate over the cube gives two Gaussian
    # interval factors and one exponential difference per vector component.
    center = np.asarray(point) / side
    lo, hi = -0.5 - center, 0.5 - center
    nodes, weights = leggauss(order)
    acc = np.zeros(3)
    for start in range(-24, 42, 2):
        for u, weight in zip(nodes + start + 1, weights):
            t = math.exp(u)
            intervals = [math.sqrt(math.pi) / (2 * t) *
                         (math.erf(t * h) - math.erf(t * l)) for l, h in zip(lo, hi)]
            for i in range(3):
                a, b = (t * lo[i]) ** 2, (t * hi[i]) ** 2
                diff = math.exp(-min(a, b)) * -math.expm1(-abs(a - b))
                if a > b:
                    diff = -diff
                acc[i] += weight * t * diff * intervals[(i + 1) % 3] * intervals[(i + 2) % 3]
    return acc * (2 / math.sqrt(math.pi)) * mu / side**2


def main():
    normalized = [(0,0,0), (.13,-.21,.07), (.499999,-.19,.09),
                  (.499999,.499999,.499999), (-.37,.38,-.48),
                  (.5,0,0), (.5,.5,0), (.5,.5,.5), (.5,.17,-.26),
                  (.500001,.499999,.05), (.500000001,.13,.09),
                  (.52,.11,-.32), (.51,.51,.51), (.7,0,0), (-1.3,.39,-.17),
                  (.01,-.63,1.19), (0,2.4,-.11), (-2.8,1.4,1.9), (4,0,0)]
    worst = 0.0
    for point in normalized:
        for side, mu in [(1., 1.), (2.3, .065)]:
            p = np.array(point) * side
            coarse = reference(p, side, mu, 14)
            fine = reference(p, side, mu, 28)
            limit = 5e-8 * mu / side**2 + 2e-8 * np.linalg.norm(fine)
            error = np.linalg.norm(cube_accel(p, side, mu) - fine)
            if not (np.linalg.norm(coarse - fine) < limit / 10 and error <= limit):
                raise ValueError(f"Verification failed at {p.tolist()}, side={side}, mu={mu}")
            worst = max(worst, float(error))
    print(json.dumps({"cases": 2*len(normalized), "max_error_km_s2": worst, "passed": True}))


if __name__ == "__main__":
    main()
