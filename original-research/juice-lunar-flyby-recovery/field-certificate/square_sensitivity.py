"""Seek witnesses ruling out every roll angle of the 703 square under alternative caps."""
from common import *
from scipy.optimize import minimize,least_squares
from itertools import combinations
a=Arc("juice");n=len(a.rows)
x0=np.array(load("juice-forecast.json")["interior_state"])
data=json.loads((HERE/"juice-comparison.json").read_text())
J=least_squares(lambda p:a.residual(x0+a.delta(p)).ravel(),np.zeros(6),jac="3-point",max_nfev=1).jac
_,s,vt=np.linalg.svd(J,full_matrices=False)
def state(v):return x0+a.delta((v/s)@vt)
def params(x):
    qp=np.r_[a.axes@(x[:3]-x0[:3]),a.axes@(x[3:]-x0[3:])]/a.scale
    return (vt@qp)*s
def valid(x,cap):
    err=max(o.distance(a.sky(x,n),a.obs))
    if err>cap-.002:return None
    enc=a.encounter(x)
    if not enc["resolved"] or enc["moon_radius_km"]<=1750:return None
    ray=a.sky(x,n+1)[-1];moon=o.hpop_de440_positions_km(a.JD+a.times[-1]/86400,a.ephem)["moon"]-a.sites[-1]
    limb=o.distance(ray,o.unit(moon))/3600-np.rad2deg(np.arcsin(1737.4/np.linalg.norm(moon)))
    if limb<=.02 or ray@a.sites[-1]<=0:return None
    return dict(state=x.tolist(),los=ray.tolist(),max_as=err,encounter=enc)
out={}
for cap in (2,3):
    orig=data["grid"][str(cap)]["witnesses"]
    starts=[np.array(w["state"]) for w in orig]
    dirs=[np.array(w["los"]) for w in orig];axis=o.unit(dirs[1]-dirs[0])
    cand=[valid(x,cap) for x in starts];cand=[r for r in cand if r]
    # Constrain the trajectory to stay outside the lunar surface at periapsis. During optimization
    # evaluate distance at a moving closest-approach time and enforce zero radial velocity there.
    for sign,seed in zip((-1,1),starts):
        from datetime import datetime
        enc=a.encounter(seed);tc=datetime.fromisoformat(enc["utc"]).timestamp()-a.rows[0][0]
        def peri(v):
            t=tc+100*v[-1];x=a.propagate(state(v[:6]),t,0)
            moon=lambda tt:o.hpop_de440_positions_km(a.JD+tt/86400,a.ephem)["moon"]
            d=x[:3]-moon(t);vrel=x[3:]-(moon(t+1)-moon(t-1))/2
            return np.linalg.norm(d),o.unit(d)@vrel
        con=[dict(type="ineq",fun=lambda v:(cap-.018)**2-np.sum(a.residual(state(v[:6]))**2,axis=1)),
             dict(type="ineq",fun=lambda v:(peri(v)[0]-1751)/100),
             dict(type="eq",fun=lambda v:peri(v)[1])]
        opt=minimize(lambda v:-sign*np.rad2deg(a.sky(state(v[:6]),n+1)[-1]@axis),np.r_[params(seed),0.],method="SLSQP",
                     bounds=[(-32,32)]*6+[(-60,60)],constraints=con,options=dict(maxiter=80,eps=3e-5,ftol=2e-8))
        for f in (1.,.998,.99,.98,.95,.9):
            r=valid(seed+f*(state(opt.x[:6])-seed),cap)
            if r:
                cand.append(r);break
    pair=max(combinations(cand,2),key=lambda p:o.distance(np.array(p[0]["los"]),np.array(p[1]["los"])))
    out[cap]=dict(witnesses=pair,diameter_deg=o.distance(np.array(pair[0]["los"]),np.array(pair[1]["los"]))/3600)
    save("square-sensitivity.json",out);print(cap,out[cap]["diameter_deg"],flush=True)
