"""Independent checks in a radial/transverse/normal basis, using km and s."""
import math
import numpy as np
from dust_bridge import acceleration

def main():
    mu = 132712440018.0
    au = 149597870.7
    c = 299792.458
    count = 0
    for direction in ([1.,0.,0.], [2.,-3.,5.]):
        radial = np.array(direction) / np.linalg.norm(direction)
        tangent = np.cross(radial, [0.,0.,1.])
        tangent /= np.linalg.norm(tangent)
        normal = np.cross(radial, tangent)
        for distance in (0.3*au, au, 2*au):
            g = mu / distance**2
            for beta in (0., 0.03, 0.4):
                for vr, vt, vn in ((0.,0.,0.), (0.,30.,0.), (0.,-30.,0.), (12.,0.,0.), (-12.,0.,0.), (-8.,21.,-11.)):
                    r = distance * radial
                    v = vr * radial + vt * tangent + vn * normal
                    actual = acceleration(r, v, mu, beta, c)
                    components = np.array([np.dot(actual, axis) for axis in (radial,tangent,normal)]) / g
                    expected = [-1+beta-2*beta*vr/c, -beta*vt/c, -beta*vn/c]
                    np.testing.assert_allclose(components, expected, rtol=1e-12, atol=1e-14)
                    count += 1
    print(f"Passed {count} independent force-vector checks.")
if __name__=='__main__':
    main()
