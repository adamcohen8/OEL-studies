"""Shared-time OEL RK4 histories for evaluating many affine verification cells."""
from affine_cells import *
from functools import lru_cache
engine=o.OrbitPropagator(integrator='rk4',numeric_backend='rust',plugins=a.engine.plugins)
@lru_cache(maxsize=1)
def grid():
    tt=[a0]
    for stop,step in ((tr-13*3600,120.),(tr-5*3600,3.),(tr-1,30.)):
        tt.extend(np.linspace(tt[-1],stop,int(np.ceil((stop-tt[-1])/step))+1)[1:])
    t=np.array(tt)
    body=[o.hpop_de440_positions_km(a.JD+z/86400,a.ephem) for z in t]
    moon=np.array([b['moon'] for b in body]);sun=np.array([b['sun'] for b in body])
    dr=body[-1]['earth']-o.hpop_de440_positions_km(a.JD+tr/86400,a.ephem)['earth']
    return t,moon,sun,dr
def fasttrack(last):
    t,moon,sun,dr=grid();h=np.diff(t)
    x=engine.try_propagate_sampled_steps(last,h,range(1,len(t)),t[0],o.zero,a.env,o.ctx,step_times=t[:-1])
    if x is None:raise RuntimeError('OEL sampled history unavailable')
    r=x[:,:3];bodies=[r,r-moon,r-sun,moon,sun]
    d=np.column_stack([np.linalg.norm(u,axis=1) for u in bodies])
    # Proxy acceleration for the a posteriori quadrature test. J2 is included in OEL's
    # propagation and bounded in the residual. The stored knots need not be exact.
    acc=np.zeros_like(r);J=np.zeros((len(x),3,3))
    for k,mu in enumerate([E,M,S,M,S]):
        rel=bodies[k];rd=d[:,k]
        acc-=mu*rel/rd[:,None]**3
        if k<3:
            J+=mu*(3*rel[:,:,None]*rel[:,None,:]/rd[:,None,None]**5-np.eye(3)/rd[:,None,None]**3)
    return t-a0,x,acc,d,J,r[-1]-site+dr,(tr-t[-1])*o.LIGHT_SPEED
if __name__=='__main__':
    from outer_test import path
    out=[]
    for j in (0,11):
        last=path(np.array(nom['accepted'][j]['state']))[6*len(a.obs)-6:6*len(a.obs)]
        start=perf_counter();b=fasttrack(last);elapsed=perf_counter()-start
        other=track(last);q=test(b,np.zeros(6));p=test(other,np.zeros(6))
        print(len(b[0]),elapsed,q,p,flush=True)
        out.append(dict(reference=j,seconds=elapsed,samples=len(b[0]),rk4=q,rkf78=p))
    save('reference-check.json',out)
