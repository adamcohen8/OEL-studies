"""Direction change across a single nonimpacting near-lunar passage with bounded perturbations."""
from math import sqrt,atan,pi
from functools import lru_cache
from tubes import E,M,S
def ext(r,extra):
    dm=360000.
    if r>=dm-6400:return 1.
    return (E*((dm-r)**-2-dm**-2)+5e-8+M/340000**2+
            2*S*(410000+r)/(1.44e8-410000-r)**3+12*.001083*E*6379**2/(dm-r)**4+extra+1e-12)
@lru_cache(maxsize=30000)
def turn_bound(duration,slow,fast,r0,speed,extra):
    """Prove convex squared radius, then bound the near-sphere hodograph and exterior impulse."""
    rs=1737.4; emin=slow**2/2-M/r0;emax=fast**2/2-M/r0
    radius=r0+.1;k=float("inf");outer_imp=0.
    h=120.;t=0.
    while t<duration:
        step=min(h,duration-t)
        ar=ext(radius+(speed+2)*step+1.,extra)
        if ar>=1:return None
        work=step*ar*(speed+2)
        emax+=min(work,step*ar*(sqrt(max(0,2*(emax+M/rs)))+step*ar/2))+1e-9
        emin-=min(work,step*ar*(sqrt(max(0,2*(emax+M/rs)))))+1e-9
        radius+=step*min(speed+2,sqrt(max(0,2*(emax+M/r0))))+.01
        k=min(k,2*emin+M/radius-radius*ar-.0001)
        if k<=0:return None
        outer_imp+=step*ar
        t+=step
    # A positive k over the entire horizon prevents a second inward crossing.
    ar=ext(r0+.1,extra); tin=min(duration,2*r0*fast/k+.1)
    e0=slow**2/2-M/r0-tin*ar*(speed+2)-.00001
    if e0<=0:return None
    vlocal=min(speed+2,sqrt(2*(fast**2/2-M/r0+tin*ar*(speed+2)+M/rs))+.001)
    torque=r0*ar*tin+.001
    hmin=sqrt(2*e0*rs**2+2*M*rs)-torque
    if hmin<=0:return None
    # e_vec = v x h / mu - r_hat; |mu e_dot| <= (|h|+r|v|)|f|.
    ecc=2*r0*vlocal*ar*tin+.001
    dv_near=sqrt((2*M+ecc+vlocal*torque)**2+(vlocal*torque)**2)/hmin
    dv_out=M/(r0*sqrt(k))*atan(sqrt(k)*duration/r0)+outer_imp
    return dv_near+dv_out+.0001
