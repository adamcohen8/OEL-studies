from common import *
from comparison import second_tess
from scipy.optimize import minimize
from itertools import combinations
# Only definitions and setup before the direction search, with no operations on archived files.
ns={}
exec((HERE/"nuisance.py").read_text().split("for sign in (-1,1):")[0],ns)
a=ns["a"];ii=ns["ii"];stations=ns["stations"];n=len(a.rows)
r0=json.loads((HERE/"square-sensitivity.json").read_text())
out={}; union=[]
for cap,rec in r0.items():
    rows=[]
    for w in rec["witnesses"]:
        x=np.array(w["state"])
        rays=ns["directions"](x,np.zeros(4),True)
        g=a.geometry(x);enc=a.encounter(x)
        assert a.physical(g) and enc["moon_radius_km"]>1750
        residual=max(o.distance(rays[:-1],a.obs))
        moon=np.array([o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)["moon"] for t in a.times])-a.sites
        limb=o.distance(rays,o.unit(moon))/3600-np.rad2deg(np.arcsin(1737.4/np.linalg.norm(moon,axis=1)))
        elevation=np.rad2deg(np.arcsin(np.sum(rays*o.unit(a.sites),axis=1)))
        assert residual<float(cap) and min(limb)>0 and min(elevation)>0
        rr=dict(state=x.tolist(),los=rays[-1].tolist(),max_as=residual,min_limb_deg=min(limb),min_elevation_deg=min(elevation),encounter=enc,coast=g)
        rows.append(rr)
        starts=rays
        def fit(x):
            d0=ns["directions"](x,np.zeros(4),True)
            rr=ns["project"](d0)
            rates=(ns["project"](ns["directions"](x,np.ones(4)))-ns["project"](ns["directions"](x,-np.ones(4))))/2
            b=[];dt=[]
            for j in range(4):
                def cons(q):
                    return np.r_[q[3]-np.linalg.norm(rr[ii==j]+rates[ii==j]*q[2]-q[:2],axis=1),.999**2-q[:2]@q[:2]]
                opt=minimize(lambda q:q[3],np.array([0,0,0,4.]),method="SLSQP",bounds=[(-1,1)]*3+[(0,100)],
                             constraints=dict(type="ineq",fun=cons),options=dict(ftol=1e-10,maxiter=80))
                b.append(opt.x[:2]);dt.append(opt.x[2])
            b=np.array(b);dt=np.array(dt)
            u=ns["directions"](x,dt,True)
            ar=(b[ii,0,None]*a.east+b[ii,1,None]*a.north)/o.ASEC
            mag=np.linalg.norm(ar,axis=1)[:,None];obs=np.cos(mag)*a.obs+np.sinc(mag/np.pi)*ar
            err=max(o.distance(u[:-1],obs))
            return dict(state=x.tolist(),station_bias_as=b.tolist(),station_clock_s=dt.tolist(),max_as=err,los=u[-1].tolist())
        if cap=="3":
            middle=ns["origin"];candidate=None
            for alpha in (1,.99,.98,.96,.94,.9):
                xnew=middle+alpha*(x-middle);fitrow=fit(xnew)
                if fitrow["max_as"]<1.999:
                    fitrow["encounter"]=a.encounter(xnew)
                    if fitrow["encounter"]["moon_radius_km"]>1750:
                        union.append(fitrow);break
    angle=o.distance(*[np.array(w["los"]) for w in rows])/3600
    out[cap]=dict(witnesses=rows,diameter_deg=angle)
prior=json.loads((HERE/"nuisance.json").read_text())["witnesses"]
pair=max(combinations(prior+union,2),key=lambda p:o.distance(np.array(p[0]["los"]),np.array(p[1]["los"])))
out["shared"]=dict(stations=stations,witnesses=pair,diameter_deg=o.distance(*[np.array(w["los"]) for w in pair])/3600)
save("expanded-witnesses.json",out)
print({k:v["diameter_deg"] for k,v in out.items()})
