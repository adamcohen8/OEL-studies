"""Shared clock and angular bias witnesses. Validate fixed shifts with fresh OEL emission times."""
from scipy.optimize import least_squares,minimize
from itertools import combinations
from common import *
a=Arc("juice");n=len(a.rows)
stations=sorted(set(r[1] for r in a.rows));ii=np.array([stations.index(r[1]) for r in a.rows])
origin=np.array(load("juice-forecast.json")["interior_state"])
J=least_squares(lambda p:a.residual(origin+a.delta(p)).ravel(),np.zeros(6),jac="3-point",max_nfev=1).jac
_,s,V=np.linalg.svd(J,full_matrices=False)
def state(v):return origin+a.delta((v[:6]/s)@V)
def unshifted(x):return a.residual(x)
def directions(x,offset,exact=False):
    old=a.engine
    if exact:a.engine=o.OrbitPropagator(integrator="dopri5",plugins=old.plugins,adaptive_rtol=5e-13,adaptive_atol=3e-12)
    res=[]
    last=0.;curr=x.copy()
    for j,t0 in enumerate(a.times):
        sitecode=a.rows[j][1] if j<n else "703"
        t=t0+(offset[ii[j]] if j<n else 0.)
        site=a.site(sitecode,t)
        def position(clock):
            if clock<0:
                ac=a.engine.acceleration_at(t_s=0,x_eci=x,command_accel_eci_km_s2=o.zero,env=a.env,ctx=o.ctx)
                return x[:3]+clock*x[3:]+.5*clock**2*ac
            if clock==last:return curr[:3]
            return a.propagate(curr,clock-last,last)[:3]
        if t>last+15:
            curr=a.propagate(curr,t-15-last,last);last=t-15
        at=position(t)
        earth=lambda t:o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)["earth"]
        lag=np.linalg.norm(at-site)/o.LIGHT_SPEED
        for _ in range(3):
            ray=position(t-lag)-site+earth(t-lag)-earth(t)
            lag=np.linalg.norm(ray)/o.LIGHT_SPEED
        res.append(o.unit(ray))
        if t>last:curr=a.propagate(curr,t-last,last);last=t
    a.engine=old
    return np.array(res)
def project(d):
    return o.ASEC*np.column_stack((np.sum(d[:n]*a.east,axis=1),np.sum(d[:n]*a.north,axis=1)))/np.sum(d[:n]*a.obs,axis=1)[:,None]
rates=(project(directions(origin,np.ones(len(stations))))-project(directions(origin,-np.ones(len(stations)))))/2
def nuisance(v):
    q=v[6:].reshape(len(stations),3)
    return q[:,:2],q[:,2]
def cons(v):
    b,t=nuisance(v)
    rr=unshifted(state(v))+rates*t[ii,None]-b[ii]
    return np.r_[1.97**2-np.sum(rr**2,axis=1),.995**2-np.sum(b**2,axis=1)]
base=a.sky(origin,n+1)[-1];e,nn=o.basis(base)
axis=o.unit(np.array(load("juice-clearance-support.json")["independent_replay"][1]["los"])-base)
axis=o.unit(axis-(axis@base)*base)
start=np.zeros(6+3*len(stations));out=[]
for sign in (-1,1):
    opt=minimize(lambda v:-sign*np.rad2deg(a.sky(state(v),n+1)[-1]@axis),start,method="SLSQP",
        bounds=[(-30,30)]*6+[(-1,1)]*(3*len(stations)),constraints=[dict(type="ineq",fun=cons)],
        options=dict(maxiter=75,ftol=1e-8,eps=1e-4))
    def verify(v):
        x=state(v);enc=a.encounter(x)
        if not enc["resolved"] or enc["moon_radius_km"]<1760:return None
        b,dt=nuisance(v)
        rays=directions(x,dt,True)
        # Exact biased reference unit direction via the exponential map in each tangent frame.
        angular=(b[ii,0,None]*a.east+b[ii,1,None]*a.north)/o.ASEC
        mag=np.linalg.norm(angular,axis=1)[:,None]
        observed=np.cos(mag)*a.obs+np.sinc(mag/np.pi)*angular
        err=o.distance(observed,rays[:-1])
        moon=np.array([o.hpop_de440_positions_km(a.JD+(t+(dt[ii[j]] if j<n else 0))/86400,a.ephem)["moon"] for j,t in enumerate(a.times)])
        sites=np.array([a.site(a.rows[j][1],t+dt[ii[j]]) if j<n else a.sites[-1] for j,t in enumerate(a.times)])
        limb=o.distance(rays,o.unit(moon-sites))/3600-np.rad2deg(np.arcsin(1737.4/np.linalg.norm(moon-sites,axis=1)))
        elev=np.rad2deg(np.arcsin(np.sum(rays*o.unit(sites),axis=1)))
        geom=a.geometry(x)
        if max(err)>2 or min(limb)<=0 or min(elev)<=0 or not a.physical(geom):return None
        return dict(state=x.tolist(),station_bias_as=b.tolist(),station_clock_s=dt.tolist(),max_as=max(err),
                    min_elevation_deg=min(elev),min_lunar_limb_deg=min(limb),los=rays[-1].tolist(),encounter=enc,geometry=geom)
    low,hi=0,1;winning=verify(opt.x)
    if winning is None:
        for _ in range(12):
            mid=(hi+low)/2
            w=verify(mid*opt.x)
            if w is None:hi=mid
            else:low=mid;winning=w
    if winning:out.append(winning)
    print(sign,opt.success,len(out),flush=True)
if len(out)==2:
    angle=o.distance(*[np.array(w["los"]) for w in out])/3600
    save("nuisance.json",dict(stations=stations,witnesses=out,diameter_lower_deg=angle))
    print(angle,flush=True)
