"""Independent witness audit with OEL DOPRI5 and successive emission-time solutions."""
import json
from pathlib import Path
import numpy as np
import forecast as f
from forecast import HERE, o

def validate(target,pair):
    a=f.Arc(target)
    a.engine=o.OrbitPropagator(integrator="dopri5", plugins=a.engine.plugins, adaptive_rtol=5e-13,adaptive_atol=3e-12)
    def bodies(t):
        return o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)
    observations = []
    for witness in pair:
        start=np.array(witness["state"])
        x=start.copy()
        rays=[]
        last=0.
        for i,t in enumerate(a.times):
            if t>last:
                # Advance to within 15 s of reception. The full light-time is < 6 s for these witnesses.
                pre=t-15.
                if pre>last:
                    x=a.propagate(x,pre-last,last)
                    last=pre
                at_t=a.propagate(x,t-last,last)
            else:
                at_t=x
            lag=np.linalg.norm(at_t[:3]-a.sites[i])/o.LIGHT_SPEED
            for _ in range(4):
                emit=t-lag
                if emit<0:
                    acc=a.engine.acceleration_at(t_s=0.,x_eci=start,command_accel_eci_km_s2=o.zero,env=a.env,ctx=o.ctx)
                    position=start[:3]+emit*start[3:]+0.5*emit**2*acc
                else:
                    assert emit>=last
                    position=a.propagate(x,emit-last,last)[:3]
                ray=position-a.sites[i]+(bodies(emit)["earth"]-bodies(t)["earth"])
                lag=np.linalg.norm(ray)/o.LIGHT_SPEED
            rays.append(o.unit(ray))
            x=at_t
            last=t
        rays=np.array(rays)
        lunar=np.array([bodies(t)["moon"] for t in a.times])-a.sites
        limb=o.distance(rays,o.unit(lunar))/3600-np.rad2deg(np.arcsin(1737.4/np.linalg.norm(lunar,axis=1)))
        alt=np.rad2deg(np.arcsin(np.sum(rays*o.unit(a.sites),axis=1)))
        observations.append(dict(max_fit_as=max(o.distance(rays[:-1],a.obs)),score_los=rays[-1].tolist(),
                                 min_limb_clearance_deg=min(limb),min_elevation_deg=min(alt)))
    return dict(witnesses=observations,lower_diameter_deg=o.distance(*[np.array(v["score_los"]) for v in observations])/3600)

out={}
for name,filename,key in [("juice","juice-clearance-support.json","independent_replay"),("tess","tess-forecast.json","witnesses")]:
    result=validate(name,json.loads((HERE/filename).read_text())[key])
    assert all(w["max_fit_as"]<2 and min(w["min_limb_clearance_deg"],w["min_elevation_deg"])>0 for w in result["witnesses"]),result
    out[name]=result
    print(name,result,flush=True)
(HERE/"emission-checks.json").write_text(json.dumps(out,indent=2)+"\n")
