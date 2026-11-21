"""Enclose each historical measurement directly as a function of the terminal cell."""
from fast_reference import *
from outer_sky import reverse
from outer_test import ref,n,N
from affine_cells import tube
@lru_cache(maxsize=1)
def optical_grid():
    tr=a.times[:n];rho=ref[6*n:];te=tr-rho/o.LIGHT_SPEED
    t=np.unique(np.r_[np.arange(0,tr[-1],150.),tr,te,-10.])
    b=[o.hpop_de440_positions_km(a.JD+q/86400,a.ephem) for q in t]
    return t,te,np.array([q['moon'] for q in b]),np.array([q['sun'] for q in b]),np.array([q['earth'] for q in b])
def constraints(last,width):
    t,te,moon,sun,earth=optical_grid()
    x0=state(reverse(last))
    k0=np.searchsorted(t,0.);times=t[k0:];h=np.diff(times)
    hist=engine.try_propagate_sampled_steps(x0,h,range(1,len(times)),0.,o.zero,a.env,o.ctx,step_times=times[:-1])
    # Reference only for the two short pre-epoch nodes. The a posteriori defect below
    # accounts for the approximation at those nodes and for shooting error at reception.
    acc0=a.engine.acceleration_at(t_s=0.,x_eci=x0,command_accel_eci_km_s2=o.zero,env=a.env,ctx=o.ctx)
    neg=t[:k0,None];pre=np.column_stack([x0[:3]+neg*x0[3:]+neg**2*acc0/2,x0[3:]+neg*acc0])
    x=np.vstack([pre,hist]);x[-1]=last
    r=x[:,:3];rel=[r,r-moon,r-sun,moon,sun];d=np.column_stack([np.linalg.norm(z,axis=1) for z in rel])
    acc=np.zeros_like(r);J=np.zeros((len(t),3,3))
    for k,mu in enumerate([E,M,S,M,S]):
        z=rel[k];rd=d[:,k]
        acc-=mu*z/rd[:,None]**3
        if k<3:J+=mu*(3*z[:,:,None]*z[:,None,:]/rd[:,None,None]**5-np.eye(3)/rd[:,None,None]**3)
    X=x[::-1].copy();X[:,3:]*=-1;B=inverse@np.diag(width);B[3:]*=-1
    where=np.searchsorted(t,te);ids=(len(t)-1-where)[::-1].copy()
    Ps=np.zeros((n,3,6));errs=np.full(n,-1.)
    tube((t[-1]-t[::-1]).copy(),X,acc[::-1].copy(),d[::-1].copy(),J[::-1].copy(),B,ids,Ps,errs)
    assert min(errs)>=0
    Ps=Ps[::-1];errs=errs[::-1]
    box=np.load(D/'contracted.npz')['bounds'];alpha=2./o.ASEC;rows=[];rhs=[]
    ang=np.linspace(0,2*np.pi,64,endpoint=False)
    for j,i in enumerate(where):
        pos=x[i,:3]-a.sites[j]+earth[i]-earth[np.searchsorted(t,a.times[j])]
        rho0=(a.times[j]-t[i])*o.LIGHT_SPEED
        eps=errs[j]+.021
        lag=(abs(np.linalg.norm(pos)-rho0)+maxball(Ps[j])+eps)/(o.LIGHT_SPEED-50)
        maxrange=rho0+lag*o.LIGHT_SPEED
        vel=x[i,3:]+a.earth_v[j]
        dv=np.linalg.norm(np.max(abs(box[j*6+3:j*6+6]-x[i,3:,None]),axis=1))+.004
        assert dv<2
        tangent=a.east[j]*np.cos(ang)[:,None]+a.north[j]*np.sin(ang)[:,None]
        beta=(tangent@vel)/(o.LIGHT_SPEED+a.obs[j]@vel)
        kvec=tangent-beta[:,None]*a.obs[j]
        slack=(np.sin(alpha)+abs(beta)*alpha**2/2)*maxrange+np.linalg.norm(kvec,axis=1)*(eps+dv*lag+.00004*lag**2)+1e-5
        rows.extend(kvec@Ps[j]);rhs.extend(slack-kvec@pos-beta*rho0)
    return np.array(rows),np.array(rhs),dict(max_reference_and_affine_radius_km=max(errs),observations=n)
if __name__=='__main__':
    from adaptive import *
    z=json.loads((D/'joint-cover.json').read_text())['leaves'];q=next(k for k in z if k['result']['status']=='unresolved')
    b=np.array(q['box']);last=origin+inverse@b.mean(axis=1)
    A,h,info=constraints(last,np.ptp(b,axis=1)/2);print(A.shape,info,min(h),flush=True)
