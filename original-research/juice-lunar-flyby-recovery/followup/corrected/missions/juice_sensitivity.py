"""Explore cap sensitivity on a line between the interior fit and each accepted witness."""
import json
import numpy as np
from scipy.optimize import brentq
from forecast import Arc, o, HERE
a=Arc("juice")
w=json.loads((HERE/"juice-clearance-support.json").read_text())["independent_replay"]
interior=np.array(json.loads((HERE/"juice-forecast.json").read_text())["interior_state"])
def maxres(x):
    return max(o.distance(a.sky(x,len(a.rows)),a.obs))
results=[]
for cap in (0.65,1.,1.25,1.5,1.75,2.):
    dirs=[]
    geometries=[]
    for w0 in w:
        x0=np.array(w0["state"])
        def state(v):
            return interior+v*(x0-interior)
        alpha=1. if maxres(x0)<cap-0.01 else brentq(lambda v:maxres(state(v))-cap+.01,0,1.)
        x=state(alpha)
        ray=a.sky(x,len(a.rows)+1)[-1]
        dirs.append(ray)
        enc=a.encounter(x)
        geom=a.geometry(x)
        moon=o.hpop_de440_positions_km(a.JD+a.times[-1]/86400,a.ephem)["moon"]-a.sites[-1]
        elev=np.rad2deg(np.arcsin(ray@o.unit(a.sites[-1])))
        limb=o.distance(ray,o.unit(moon))/3600-np.rad2deg(np.arcsin(1737.4/np.linalg.norm(moon)))
        assert enc["resolved"] and enc["moon_radius_km"]>1937.4 and geom["min_earth_km"]>6500 and min(elev,limb)>0
        geometries.append(dict(state=x.tolist(),los=ray.tolist(),fit_max_as=maxres(x),encounter=enc,
                               min_earth_sample_km=geom["min_earth_km"],limb_deg=limb,el_deg=elev))
    results.append(dict(cap_as=cap,pair_separation_deg=o.distance(*dirs)/3600,witnesses=geometries))
print([(r["cap_as"],r["pair_separation_deg"]) for r in results])
(HERE/"juice-cap-sensitivity.json").write_text(json.dumps(results,indent=2)+"\n")
