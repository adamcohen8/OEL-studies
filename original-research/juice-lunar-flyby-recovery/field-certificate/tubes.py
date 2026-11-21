"""Bound continuous solutions about OEL reference knots using acceleration quadrature defects."""
from common import *
from sim.dynamics.orbit.environment import EARTH_MU_KM3_S2 as E, MOON_MU_KM3_S2 as M, SUN_MU_KM3_S2 as S
from math import asin, degrees
from time import perf_counter

def prepare(a,state,step_scale,extra_nodes=()):
    data=[]
    knots=np.unique(np.r_[extra_nodes,a.times[-1]])
    x=np.array(state)
    t=0.
    while True:
        bodies=o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)
        acc=a.engine.acceleration_at(t_s=t,x_eci=x,command_accel_eci_km_s2=o.zero,env=a.env,ctx=o.ctx)
        norms=[np.linalg.norm(z) for z in (x[:3],x[:3]-bodies["moon"],x[:3]-bodies["sun"],bodies["moon"],bodies["sun"])]
        data.append((t,x.copy(),acc,norms))
        if t>=a.times[-1]: break
        step=min(a.times[-1]-t,step_scale*min(240.,max(1.,norms[1]/1200.)))
        boundary=knots[np.searchsorted(knots,t,side="right")]
        next_t=min(t+step,boundary)
        x=a.propagate(x,next_t-t,t)
        t=next_t
    return data

def enclosure(data,a,radius,disturb,force_law=None,record=None):
    R,V=radius,radius/100000
    peak=0.
    for (t,x,b,d),(tn,y,c,_) in zip(data,data[1:]):
        h=tn-t
        # Bootstrap a tube about the first endpoint covering the whole step.
        speed=np.linalg.norm(x[3:])+V+1e-8
        disp=R+h*(speed+0.015*h/2)+1e-5
        e,m,s,bm,bs=np.array(d)-[disp,disp+2*h,disp+35*h,2*h+1e-5,35*h+1e-5]
        if min(e-6379,m-1737.4,s-1e8,bm-300000,bs-1e8)<=0:
            return dict(status="unresolved",at_hours=t/3600,reason="distance tube intersects a domain boundary")
        j2=12*.001083*E*6379**2/e**4
        local_extra=disturb+(force_law(m)+force_law(bm) if force_law else 0.)
        A=E/e**2+M/m**2+M/bm**2+S/s**2+S/bs**2+j2+local_extra
        if A>=.015: return dict(status="unresolved",at_hours=t/3600,reason="acceleration bootstrap")
        v=speed+A*h
        # Norm bound on the second time derivative of each inverse-square acceleration.
        snap=sum(mu*(24*speed_bound**2/floor**4+2*acc_bound/floor**3)
                 for mu,floor,speed_bound,acc_bound in ((E,e,v,A),(M,m,v+2,A+.00002),
                 (M,bm,2,.00002),(S,s,v+35,A+.00002),(S,bs,35,.00002)))
        L=2*(E/e**3+M/m**3+S/s**3)
        # 10 m envelope for numeric body-position evaluation and bounded endpoint force arithmetic.
        extra=2*j2+local_extra+.01*2*(M/m**3+M/bm**3+S/s**3+S/bs**3)+1e-16
        defect_r=np.linalg.norm((y[:3]-x[:3])-h*x[3:]-h*h*(2*b+c)/6)+1e-8
        defect_v=np.linalg.norm((y[3:]-x[3:])-h*(b+c)/2)+1e-13
        den=1-h*h*L/6
        if den<=0: return dict(status="unresolved",at_hours=t/3600,reason="positive denominator needed")
        new=(R+h*V+h*h*(2*L*R)/6+h*h*extra/2+h**4*snap/24+defect_r)/den
        new=new*(1+1e-12)+1e-12
        V=(V+h*(L*(R+new)/2+extra)+h**3*snap/12+defect_v)*(1+1e-12)+1e-16
        R=new
        peak=max(peak,R)
        if record is not None:
            record[tn]=(R,V,e,m,v)
    x=data[-1][1]
    distance=np.linalg.norm(x[:3]-a.sites[-1])
    # For the short light-time interval use barycentric spacecraft speed <=50 km/s.
    # Earlier knots cover this interval. This gives a cap centered on instantaneous geometric direction.
    if np.linalg.norm(x[3:])+V+35+1>50 or distance<=R: return dict(status="unresolved",reason="light time")
    lag=(distance+R)/(o.LIGHT_SPEED-50)
    q=(R+50*lag+.01)/distance
    if q>=1: return dict(status="unresolved",reason="projection")
    return dict(status="bounded",end_position_radius_km=R,end_velocity_radius_km_s=V,
                field_radius_deg=degrees(asin(q)),center=o.unit(x[:3]-a.sites[-1]).tolist(),peak_radius_km=peak)

if __name__=="__main__":
    a=Arc("juice")
    out=[]
    for j,w in enumerate(load("juice-clearance-support.json")["independent_replay"]):
        for step_scale in (1.,.5):
            data=prepare(a,w["state"],step_scale)
            for r in (.001,.01,.1,1.,10.):
                for force in (0.,1e-10,1e-9):
                    out.append(dict(witness=j,step_scale=step_scale,steps=len(data)-1,initial_radius_km=r,
                                  acceleration_km_s2=force,**enclosure(data,a,r,force)))
            print(j,step_scale,[(v["initial_radius_km"],v.get("field_radius_deg"),v.get("at_hours")) for v in out[-15::3]],flush=True)
    save("tubes.json",out)
