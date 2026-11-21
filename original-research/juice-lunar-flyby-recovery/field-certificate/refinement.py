"""Linear outer relaxations for bilinear light time and integrated interval accelerations."""
import numpy as np
from tubes import E,M,S
from common import o
def add_light(rows,rhs,a,j,rlo,rhi,vlo,vhi,speed,cap,bias,timing):
    n=len(a.rows)
    # Earth barycentric acceleration and object acceleration over emission lookback.
    rmax=rhi+20000
    tide=2*S*rmax/(1.44e8-rmax)**3
    amax=E/2e5**2+M/1e5**2+M/340000**2+tide+.00002
    for ax in range(3):
        for sign, corners in [(1,[(rhi,vlo[ax]),(rlo,vhi[ax])]),(-1,[(rlo,vlo[ax]),(rhi,vhi[ax])])]:
            for rc,vc in corners:
                q=np.zeros(7*n)
                q[6*j+ax]=sign
                q[6*j+3+ax]=-sign*rc/o.LIGHT_SPEED
                q[6*n+j]=-sign*(a.obs[j,ax]+(a.earth_v[j,ax]+vc)/o.LIGHT_SPEED)-(cap+bias)/o.ASEC
                pad=.05+(speed+.6)*timing+amax*(rhi/o.LIGHT_SPEED)**2
                rows.append(q);rhs.append(sign*(a.sites[j,ax]-rc*vc/o.LIGHT_SPEED)+pad)
def signed_gravity(low,high,mu,floor):
    nearest=np.maximum(0,np.maximum(low,-high))
    rmin=max(floor,np.linalg.norm(nearest)-.1)
    rmax=np.linalg.norm(np.maximum(abs(low),abs(high)))+.1
    if rmax<rmin: rmax=rmin
    vs=np.array([-mu*z/r**3 for z in [low,high] for r in [rmin,rmax]])
    return vs.min(axis=0),vs.max(axis=0)
def add_forces(rows,rhs,a,j,low,high,coarse,extra):
    ts=a.times
    h=ts[j+1]-ts[j]
    n=len(a.rows)
    pieces=max(1,int(np.ceil(h/3600)))
    delta=h/pieces
    s0=np.zeros((3,3));s1=s0.copy()
    for p in range(pieces):
        a0=p/pieces;a1=(p+1)/pieces
        left=(1-a0)*low[j]+a0*low[j+1];right=(1-a1)*low[j]+a1*low[j+1]
        hi0=(1-a0)*high[j]+a0*high[j+1];hi1=(1-a1)*high[j]+a1*high[j+1]
        pad=coarse*h*h/8+.1
        l=np.minimum(left,right)-pad;u=np.maximum(hi0,hi1)+pad
        mid=ts[j]+(a0+a1)*h/2
        moon=o.hpop_de440_positions_km(a.JD+mid/86400,a.ephem)["moon"]
        ml,mu=moon-delta-.1,moon+delta+.1
        el,eu=signed_gravity(l,u,E,2e5)
        dl,du=signed_gravity(l-mu,u-ml,M,1e5)
        il,iu=signed_gravity(ml,mu,M,340000)
        rmax=np.linalg.norm(np.maximum(abs(l),abs(u)))
        margin=2*S*rmax/(1.44e8-rmax)**3+12*.001083*E*6379**2/(2e5)**4+extra+1e-14
        al=el+dl+il-margin;au=eu+du+iu+margin
        for k,weight in enumerate([delta,delta*h*(1-(a0+a1)/2),delta*h*(.5-(a0+a1)/2)]):
            s0[k]+=np.minimum(weight*al,weight*au)
            s1[k]+=np.maximum(weight*al,weight*au)
        # If a subinterval straddles the zero of h/2-t, account for both signs.
        if a0<.5<a1:
            correction=(au-al)*delta**2/8
            s0[2]-=correction;s1[2]+=correction
    for ax in range(3):
        i=6*j+ax
        for k in range(3):
            q=np.zeros(7*n)
            if k==0:q[i+9],q[i+3]=1,-1
            elif k==1:q[i+6],q[i],q[i+3]=1,-1,-h
            else:q[i+6],q[i],q[i+3],q[i+9]=1,-1,-h/2,-h/2
            rows.extend([q,-q]);rhs.extend([s1[k,ax]+(1e-8 if k==0 else .001),-s0[k,ax]+(1e-8 if k==0 else .001)])
