"""Exploratory constrained search for distant pairs of admissible post-encounter directions."""
import json
import numpy as np
from scipy.optimize import least_squares, minimize
from itertools import combinations
import forecast as f
from forecast import o, HERE

arc = f.Arc("juice")
data = json.loads((HERE / "juice-forecast.json").read_text())
origin = np.array(data["interior_state"])
n = len(arc.obs)
jac = least_squares(lambda p: arc.residual(origin + arc.delta(p)).ravel(), np.zeros(6), jac="3-point", max_nfev=1).jac
_, scales, axes = np.linalg.svd(jac, full_matrices=False)
def state(v):
    return origin + arc.delta((v / scales) @ axes)
def pred(v):
    return arc.sky(state(v), n + 1)[-1]
def cons(v):
    return 1.98**2 - np.sum(arc.residual(state(v))**2, axis=1)
def inspect(x):
    ray = arc.sky(x,n+1)[-1]
    enc = arc.encounter(x)
    if not enc["resolved"] or enc["moon_radius_km"] < 1937.4: return None
    earth = arc.geometry(x)["min_earth_km"]
    moon = o.hpop_de440_positions_km(arc.JD+arc.times[-1]/86400,arc.ephem)["moon"] - arc.sites[-1]
    limb = o.distance(ray, o.unit(moon))/3600 - np.rad2deg(np.arcsin(1737.4/np.linalg.norm(moon)))
    alt = np.rad2deg(np.arcsin(o.unit(arc.sites[-1]) @ ray))
    fitmax = max(o.distance(arc.sky(x,n),arc.obs))
    if min(alt,limb) <= 0 or earth < 6500 or fitmax > 1.995: return None
    return dict(state=x.tolist(), los=ray.tolist(), lunar_limb_deg=limb, elevation_deg=alt,
                fit_max_as=fitmax, min_earth_sample_km=earth, encounter=enc)
center = pred(np.zeros(6))
east, north = o.basis(center)
candidates = [inspect(np.array(row["state"])) for row in data["witnesses"]]
for theta in np.arange(0,360,45):
    axis = np.cos(np.deg2rad(theta))*east + np.sin(np.deg2rad(theta))*north
    sol = minimize(lambda v: -np.rad2deg(pred(v)@axis), np.zeros(6), method="SLSQP",
                   bounds=[(-20.,20.)]*6, constraints=[dict(type="ineq", fun=cons)],
                   options=dict(ftol=5e-8,eps=1e-4,maxiter=60))
    candidate = inspect(state(sol.x))
    if candidate is None:
        low,high=0.,1.
        for _ in range(11):
            factor=(low+high)/2
            trial=inspect(state(factor*sol.x))
            if trial is None:
                high=factor
            else:
                low=factor
                candidate=trial
    if candidate is not None:
        candidates.append(candidate)
    print(theta, sol.success, len(candidates), flush=True)
pair = max(combinations(candidates, 2), key=lambda p:o.distance(*[np.array(q["los"]) for q in p]))
arc.engine=o.OrbitPropagator(integrator="dopri5",plugins=arc.engine.plugins,adaptive_rtol=5e-13,adaptive_atol=3e-12)
pair=[inspect(np.array(q["state"])) for q in pair]
assert all(pair)
separation=o.distance(*[np.array(q["los"]) for q in pair])/3600
(HERE/"juice-clearance-support.json").write_text(json.dumps(dict(diameter_lower_deg=separation, independent_replay=pair),indent=2)+"\n")
print(separation,flush=True)
