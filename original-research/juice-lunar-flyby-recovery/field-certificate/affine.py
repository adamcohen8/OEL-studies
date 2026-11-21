"""Affine enclosure of gravity about fixed chords on the optical training arc."""
from common import np,o
from tubes import E,M,S
from functools import lru_cache
def norm(x):return np.linalg.norm(x)
def near(l,u):return norm(np.maximum(0,np.maximum(l,-u)))
@lru_cache(maxsize=10000)
def body(a,t):
    return o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)
def add(rows,rhs,a,j,low,high,ref,extra):
    n=len(a.rows);ts=a.times;h=ts[j+1]-ts[j]; pieces=int(np.ceil(h/1800))
    sums=np.zeros((3,3));errors=np.zeros(3);w0=np.zeros((3,3,3));w1=w0.copy()
    l=np.minimum(low[j],low[j+1]);u=np.maximum(high[j],high[j+1])
    # Bound absolute acceleration across the true arc using the declared physical floors.
    br=body(a,ts[j]+h/2)
    pad0=.000022*h*h/8+.1
    er=max(2e5,near(l-pad0,u+pad0));mr=max(1e5,near(l-br['moon']-h-pad0,u-br['moon']+h+pad0))
    if norm(np.maximum(abs(l),abs(u)))>4e6: return
    accel=E/er**2+M/mr**2+M/340000**2+2*S*5e6/(1.44e8-5e6)**3+1e-9+extra
    chordpad=accel*h*h/8+.05; ds=h/pieces
    for p in range(pieces):
        t=(p+.5)/pieces;now=ts[j]+h*t
        r=(1-t)*ref[j,:3]+t*ref[j+1,:3];b=body(a,now)
        xyz=np.r_[r,[0.,0.,0.]]
        # Exact OEL model evaluation on the reference chord.
        acc=a.engine.acceleration_at(t_s=now,x_eci=xyz,command_accel_eci_km_s2=o.zero,env=a.env,ctx=o.ctx)
        z=(np.arange(2)+p)/pieces
        ll=np.minimum(*[(1-k)*low[j]+k*low[j+1] for k in z])-chordpad
        uu=np.maximum(*[(1-k)*high[j]+k*high[j+1] for k in z])+chordpad
        displacement=norm(np.maximum(abs(ll-r),abs(uu-r)))
        J=np.zeros((3,3));error=extra+2*12*.001083*E*6379**2/er**4+1e-14
        for mu,center,dr,indirect_floor in ((E,np.zeros(3),0.,1.),(M,b["moon"],ds+.02,340000.),(S,b["sun"],17.5*ds+.02,1.44e8)):
            sep=r-center;length=norm(sep)
            floor=near(np.minimum(ll,r)-center-dr,np.maximum(uu,r)-center+dr)
            if floor<1.:return
            jac=mu*(3*np.outer(sep,sep)/length**5-np.eye(3)/length**3)
            J+=jac
            error+=12*mu*(displacement+dr)**2/floor**4+2*mu*(chordpad+dr)/floor**3
            if dr: error+=2*mu*dr/indirect_floor**3
        # Integrals of polynomials of degree <=2 by two-point Gauss quadrature.
        for f in (-1,1):
            s=t+f/(np.sqrt(12)*pieces)
            weights=ds/2*np.array([1,h*(1-s),h*(.5-s)])
            sums+=weights[:,None]*(acc-J@r)
            w0+=weights[:,None,None]*J*(1-s);w1+=weights[:,None,None]*J*s
        # Absolute weight integral rounded upwards, including sign reversal inside a bin.
        errors+=ds*np.array([1,h*(1-t),h*(abs(.5-t)+.5/pieces)])*(error*1.000001+1e-15)
    for ax in range(3):
        i=6*j+ax
        for k in range(3):
            q=np.zeros(7*n)
            if k==0:q[i+9],q[i+3]=1,-1
            elif k==1:q[i+6],q[i],q[i+3]=1,-1,-h
            else:q[i+6],q[i],q[i+3],q[i+9]=1,-1,-h/2,-h/2
            q[6*j:6*j+3]-=w0[k,ax];q[6*j+6:6*j+9]-=w1[k,ax]
            e=errors[k]+(1e-8 if k==0 else .001)
            rows.extend([q,-q]);rhs.extend([sums[k,ax]+e,-sums[k,ax]+e])
