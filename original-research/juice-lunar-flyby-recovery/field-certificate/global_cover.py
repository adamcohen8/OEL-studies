"""Incoming optical LP, exact dual checks, and a force-independent terminal speed cone."""
from refinement import add_light, add_forces
from affine import add as affine_gravity
from math import asin, degrees
from scipy.optimize import linprog
from common import *
from tubes import E,M,S
def norm(v):
    return np.linalg.norm(v)
def fit_cover(a,speed,cap,timing,station_bias,extra,audit_tag=None):
    ts=a.times[:-1]
    n=len(ts)
    d=7*n
    # First and last observed rays give a finite range bound given a bounded speed.
    # Reception position is observer + rho*u up to k*rho plus timestamp/site padding.
    k=(speed+35)/o.LIGHT_SPEED+(cap+station_bias)/o.ASEC+1e-7
    D=speed*(ts[-1]+2*timing)+13000
    sint=norm(np.cross(a.obs[0],a.obs[-1]))
    first=(D+k*D/(1-k))/(sint-k-k*(1+k)/(1-k))
    assert first>0 and sint>2*k/(1-k)
    ran=(first*(1+k)+D)/(1-k)+100
    assert ran*k+6500<20000
    bounds=( [(-ran-20000,ran+20000)]*3+[(-speed,speed)]*3 )*n+[(0,ran)]*n
    failures=0
    retain=False;certs=[]
    rows=[]
    rhs=[]
    def two(vec,c,w):
        rows.extend([vec,-vec]);rhs.extend([c+w,-c+w])
    def extrem(axis,idx,sign):
        nonlocal failures
        c=np.zeros(d)
        if isinstance(axis,int): c[6*idx+axis]=sign
        else: c[6*idx:6*idx+6]=sign*np.asarray(axis)
        cs=max(abs(c*scale))
        scaled_bounds=(np.array(bounds)-witness[:,None])/scale[:,None]
        res=linprog(-c*scale/cs,A_ub=smat,b_ub=sb,bounds=scaled_bounds,method="highs")
        if not res.success:
            failures+=1
            # Keep the pre-existing support of the entire variable box on numerical failure.
            return sign*(np.maximum(c*np.array(bounds)[:,0],c*np.array(bounds)[:,1]).sum()+1e-5)
        # Rigorous dual evaluation for the stored floating inequalities. Directed allowance for
        # the final dot products uses absolute product sums and finite variable bounds.
        active=np.flatnonzero(res.ineqlin.marginals<0)
        lam=-res.ineqlin.marginals[active].astype(np.longdouble)*cs/row_scales[active]
        A=mat[active].astype(np.longdouble); bb=b[active].astype(np.longdouble)
        residual=c-lam@A
        hi=np.array(bounds)[:,1];lo=np.array(bounds)[:,0]
        eps=np.finfo(np.longdouble).eps
        err=(len(active)+4)*eps*(abs(c)+lam@abs(A))
        support=np.maximum(residual*hi,residual*lo).sum()
        exact_upper=float(lam@bb+support+err@np.maximum(abs(hi),abs(lo))+(len(active)+d+8)*eps*(lam@abs(bb)+np.sum(abs(residual)*np.maximum(abs(hi),abs(lo)))))+.0000001
        if retain and audit_tag and not isinstance(axis,int):
            certs.append(dict(c=c.tolist(),active=active.tolist(),multipliers=[float(v) for v in lam],rounded_upper=exact_upper))
        return sign*exact_upper
    # Range-dependent light time and velocity terms retain the positive range variable.
    rlo=0.;rhi=ran
    for j in range(n):
        add_light(rows,rhs,a,j,rlo,rhi,np.full(3,-speed),np.full(3,speed),speed,cap,station_bias,timing)
    def dynamics(j,acc):
        h=ts[j+1]-ts[j]
        for ax in range(3):
            v=np.zeros(d);i=6*j+ax;v[i+6],v[i],v[i+3]=1,-1,-h
            two(v,0,acc*h*h/2+.001)
            v=np.zeros(d);v[i+9],v[i+3]=1,-1
            two(v,0,acc*h+1e-8)
    # Initial declared class has E >= 200,000 and Moon >= 100,000 km during reception span.
    coarse=E/2e5**2+M/1e5**2+M/340000**2+2*S*(ran+20000)/(1.44e8-ran-20000)**3+1e-9+extra
    for j in range(n-1): dynamics(j,coarse)
    base=len(rows)
    x=np.array(load("juice-forecast.json")["state"]); track=[]
    for j,t in enumerate(ts):
        if j: x=a.propagate(x,t-ts[j-1],ts[j-1])
        track.append(x.copy())
    rr=[norm(x[:3]-site)/(1+((x[3:]+ev)@u)/o.LIGHT_SPEED) for x,site,ev,u in zip(track,a.sites,a.earth_v,a.obs)]
    witness=np.r_[np.array(track).ravel(),rr]
    scale=np.r_[np.tile([1e4]*3+[.05]*3,n),np.full(n,1e4)]
    def shifted(mat,b):
        unscaled=mat*scale
        scales=np.maximum(1e-6,np.max(abs(unscaled),axis=1))
        return unscaled/scales[:,None],(b-mat@witness)/scales,scales
    for iteration in range(8):
        mat,b=np.array(rows),np.array(rhs)
        smat,sb,row_scales=shifted(mat,b)
        violation=mat@witness-b
        if np.max(violation)>.001: print("witness violated",iteration,np.argmax(violation),np.max(violation),flush=True)
        low=np.array([[extrem(ax,j,-1) for ax in range(3)] for j in range(n)])
        high=np.array([[extrem(ax,j,1) for ax in range(3)] for j in range(n)])
        def near(l,h):return norm(np.maximum(0,np.maximum(l,-h)))
        # Bound individual ranges from the state boxes; add tighter measurement inequalities at each pass.
        geo=[near(l,h) for l,h in zip(low,high)]
        newlo=max(0,(min(geo)-6500)/(1+k)-100)
        newhi=(max(norm(np.maximum(abs(l),abs(h))) for l,h in zip(low,high))+6500)/(1-k)+100
        rlo,rhi=max(rlo,newlo),min(rhi,newhi)
        del rows[base:]
        del rhs[base:]
        for j in range(n):
            vlo=np.array([extrem(ax,j,-1) for ax in (3,4,5)])
            vhi=np.array([extrem(ax,j,1) for ax in (3,4,5)])
            lo=max(0,(near(low[j]-a.sites[j],high[j]-a.sites[j])-100)/(1+k))
            hi=(norm(np.maximum(abs(low[j]-a.sites[j]),abs(high[j]-a.sites[j])))+100)/(1-k)
            for col,(lb,ub) in enumerate(zip(np.r_[low[j],vlo],np.r_[high[j],vhi])):
                old=bounds[6*j+col];bounds[6*j+col]=(max(lb,old[0]),min(ub,old[1]))
            bounds[6*n+j]=(max(lo,bounds[6*n+j][0]),min(hi,bounds[6*n+j][1]))
            add_light(rows,rhs,a,j,lo,hi,vlo,vhi,speed,cap,station_bias,timing)
        if iteration in (0,3,7,11): print("refine",speed,iteration,high[-1]-low[-1],flush=True)
        continuous_max=rhi+20000+speed*max(np.diff(ts))/2
        coarse=E/2e5**2+M/1e5**2+M/340000**2+2*S*continuous_max/(1.44e8-continuous_max)**3+1e-9+extra
        for j in range(n-1):
            add_forces(rows,rhs,a,j,low,high,coarse,extra)
            if iteration>=2: affine_gravity(rows,rhs,a,j,low,high,np.asarray(track),extra)
    mat,b=np.array(rows),np.array(rhs)
    smat,sb,row_scales=shifted(mat,b)
    low=np.array([extrem(i,n-1,-1) for i in range(6)])
    high=np.array([extrem(i,n-1,1) for i in range(6)])
    first_box=[[extrem(i,0,-1),extrem(i,0,1)] for i in range(6)]
    basis=np.array([a.obs[-1],*o.basis(a.obs[-1])])
    axes=[np.r_[z,np.zeros(3)] for z in basis]+[np.r_[np.zeros(3),z] for z in basis]
    retain=True
    rotated=[[extrem(z,n-1,-1),extrem(z,n-1,1)] for z in axes]
    if audit_tag:
        np.savez_compressed(HERE/(audit_tag+"-LP.npz"),A=mat,b=b,bounds=np.array(bounds))
        save(audit_tag+"-duals.json",certs)
    return low,high,dict(lp_support_fallbacks=failures,last_tangent_box=rotated,basis=basis.tolist(),initial_range_upper_km=ran,training_range_km=[rlo,rhi],first_state_box=first_box,last_state_lo=low.tolist(),last_state_hi=high.tolist())
def advance(a,low,high,speed,extra):
    t=a.times[-2];end=a.times[-1];x=(low+high)/2
    R=norm((high[:3]-low[:3])/2)+.01; V=norm((high[3:]-low[3:])/2)+1e-7
    best=dict(diameter_deg=180.,reason="upper hemisphere of visible sky")
    while t<=end:
        rho=norm(x[:3]-a.sites[-1])
        rad=R+speed*(end-t)+.02
        # Conservative emission shift from the final reception ball.
        rad+=(speed+35)*(rho+rad)/(o.LIGHT_SPEED-speed-35)
        if rad<rho:
            deg=2*degrees(asin(rad/rho))
            if deg<best["diameter_deg"]:
                best=dict(diameter_deg=deg,center=o.unit(x[:3]-a.sites[-1]).tolist(),hours_before=(end-t)/3600,
                          anchor_position_radius_km=R,terminal_position_radius_km=rad,anchor_speed_radius_km_s=V)
        if t>=end: break
        h=min(600.,end-t)
        # All true positions in this step lie inside the endpoint ball expanded by speed*h.
        env=o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)
        er=norm(x[:3])-R-speed*h-.01
        mr=norm(x[:3]-env["moon"])-R-(speed+2)*h-.01
        if er<10000 or mr<1800: break
        # Acceleration-distance upper bound for exact solutions; anchor uses OEL's actual displacement.
        sr=norm(env["sun"])-norm(x[:3])-R-(speed+35)*h
        tide=2*S*(norm(x[:3])+R+speed*h)/sr**3
        A=E/er**2+M/mr**2+M/340000**2+tide+12*.001083*E*6379**2/er**4+extra
        y=a.propagate(x,h,t)
        # Relative to the reference acceleration at the start, use Lipschitz position error,
        # bounded body motion and anchor drift. Both local reference defects measured explicitly.
        b=a.engine.acceleration_at(t_s=t,x_eci=x,command_accel_eci_km_s2=o.zero,env=a.env,ctx=o.ctx)
        delta=(2*E/er**3*(R+speed*h)+2*M/mr**3*(R+(speed+2)*h)
               +2*M/340000**3*2*h+2*S/sr**3*(R+(speed+35)*h)
               +2*S/(1.44e8)**3*35*h+24*.001083*E*6379**2/er**4+extra+1e-14)
        delta=min(delta,A+norm(b)+1e-13)
        R+=h*V+delta*h*h/2+norm(y[:3]-x[:3]-x[3:]*h-b*h*h/2)+1e-7
        V+=delta*h+norm(y[3:]-x[3:]-h*b)+1e-11
        x=y;t+=h
    return best
if __name__=="__main__":
    a=Arc("juice");ans=[]
    for v,cap,time,bias,force in [(6,2,0,0,0),(10,2,0,0,0),(6,2,1,1,1e-9),(10,2,1,1,1e-9)]:
        l,h,r=fit_cover(a,v,cap,time,bias,force)
        r.update(speed=v,cap=cap,timing=time,station_bias=bias,extra_force=force,cover=advance(a,l,h,v,force))
        ans.append(r); save("global.json",ans); print(v,time,r["cover"],flush=True)
