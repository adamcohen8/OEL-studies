"""Independent corner and axial spot checks of the conservative local tubes."""
from common import *
from itertools import product
a=Arc("juice")
a.engine=o.OrbitPropagator(integrator="dopri5",plugins=a.engine.plugins,adaptive_rtol=5e-13,adaptive_atol=3e-12)
bounds=json.loads((HERE/"tubes.json").read_text())
out=[]
for j,w in enumerate(load("juice-clearance-support.json")["independent_replay"]):
    x=np.array(w["state"])
    for r in (1.,10.):
        b=next(q for q in bounds if q["witness"]==j and q["initial_radius_km"]==r and q["acceleration_km_s2"]==0 and q["step_scale"]==.5)
        axes=np.r_[np.eye(3),-np.eye(3)]
        # Boundary of the Cartesian product of two Euclidean balls.
        vals=[o.distance(a.sky(x+np.r_[r*p,r*v/100000],len(a.rows)+1)[-1],np.array(b["center"]))/3600 for p,v in product(axes,axes)]
        assert max(vals)<b["field_radius_deg"]
        out.append(dict(witness=j,r=r,tested=len(vals),largest_offset_deg=max(vals),bound_deg=b["field_radius_deg"]))
save("cell-checks.json",out)
