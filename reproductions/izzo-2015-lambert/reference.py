"""Independent evaluation of Izzo (2015), Eq. 18 and Algorithm 1, for the assigned case."""
import math
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def reference(problem):
    r1, r2 = (np.array(problem[key]) for key in ("r1_km", "r2_km"))
    mu = problem["mu_km3_s2"]
    turns = problem["revolutions"]
    n1, n2, chord = np.linalg.norm(r1), np.linalg.norm(r2), np.linalg.norm(r2-r1)
    s = (chord + n1 + n2) / 2
    lam = math.sqrt(1 - chord / s) * (1 if problem["short_way"] else -1)
    time_unit = math.sqrt(s**3 / (2*mu))

    def y(x):
        return math.sqrt(1 - lam**2 * (1 - x*x))

    def flight(x):
        q = 1 - x*x
        psi = math.acos(np.clip(x*y(x) + lam*q, -1, 1))
        return ((psi + turns*math.pi)/math.sqrt(q) - x + lam*y(x)) / q

    lo, hi = -1+1e-9, 1-1e-9
    minimum = minimize_scalar(flight, bounds=(lo, hi), method="bounded", options={"xatol":1e-14})
    assert minimum.success
    target = problem["time_of_flight_s"] / time_unit
    roots = [brentq(lambda x: flight(x)-target, *bracket, xtol=1e-14) for bracket in
             [(lo, minimum.x), (minimum.x, hi)]]
    rho = (n1 - n2) / chord
    sigma = math.sqrt(1-rho*rho)
    gamma = math.sqrt(mu*s/2)
    i1, i2 = r1/n1, r2/n2
    normal = np.cross(i1,i2)
    normal /= np.linalg.norm(normal)
    if not problem["short_way"]:
        normal *= -1
    states = []
    for x in roots:
        vminus, vplus = lam*y(x)-x, lam*y(x)+x
        tangential = gamma*sigma*(y(x)+lam*x)
        v1 = (gamma*(vminus-rho*vplus)*i1 + tangential*np.cross(normal,i1))/n1
        v2 = (-gamma*(vminus+rho*vplus)*i2 + tangential*np.cross(normal,i2))/n2
        states.append({"x": x, "v1_km_s": v1.tolist(), "v2_km_s": v2.tolist()})
    return {"minimum_time_s": float(minimum.fun*time_unit), "branches": states}
