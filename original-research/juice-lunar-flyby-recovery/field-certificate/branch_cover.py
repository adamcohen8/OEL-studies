"""Global post-arc enclosure by partitioning range and radial speed, then splitting at a lunar sphere."""
from common import *
from tubes import E,M,S
from gravity_model import envelope, potential_bound
from itertools import product
from turn import turn_bound
def norm(x):return np.linalg.norm(x,axis=-1)
def gravity(l,u,mu,floor):
    nearest=np.maximum(0,np.maximum(l,-u))
    rmin=np.maximum(floor,norm(nearest)-.01)
    rmax=np.maximum(rmin,norm(np.maximum(abs(l),abs(u)))+.01)
    candidates=np.stack([-mu*x/r[...,None]**3 for x in (l,u) for r in (rmin,rmax)])
    return candidates.min(axis=0),candidates.max(axis=0)
def moon(a,t):
    return o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)["moon"]
def travel(duration,entry_v,radius,speed,force,physical,moon_floor):
    """Upper distance from moving Moon after first entry, assuming survival exterior to the lunar surface."""
    p_surface=potential_bound(1737.4) if physical else 0.
    p_entry=potential_bound(radius) if physical else 0.
    rotating_work=4e-6*1737.4*envelope(1737.4) if physical else 0.
    e=.5*entry_v**2-M/radius+p_entry
    r=radius+.02
    steps=max(1,int(np.ceil(duration/240)));h=duration/steps
    for _ in range(steps):
        wide=r+(speed+2)*h+1.
        er=moon_floor-wide
        if er<6400: return min(r+(speed+2)*(duration-_*h)+2,1e8)
        # Moon geocentric acceleration is bounded by .00002. Both solar and indirect lunar
        # accelerations are bounded on the local sphere. Absolute perturbation suffices.
        tide=2*S*(410000+wide)/(1.44e8-410000-wide)**3
        ext=E*(1/er**2-1/moon_floor**2)+12*.001083*E*6379**2/er**4+5e-8+M/340000**2+tide+force
        if physical: ext+=envelope(340000.)
        # Ext over the entire current slice; bound lunar relative speed by prior for the work term.
        e+=h*rotating_work
        e+=min(h*ext*(speed+2),h*ext*(np.sqrt(max(0,2*(e+M/1737.4+p_surface)))+.5*h*ext))+1e-10
        # Outside the entry sphere the potential magnitude is at most mu/radius.
        r+=h*min(speed+2,np.sqrt(max(0,2*(e+M/radius+p_entry))))+.00001
    return r
def solve(a,row,N,radius,physical=False,directional=False,shift=0.):
    box=np.array(row["last_tangent_box"]);basis=np.array(row["basis"])
    edges=[np.linspace(*box[i],N+1) if i in (0,3) else box[i] for i in range(6)]
    cells=list(product(*[list(zip(q[:-1],q[1:])) for q in edges]))
    arr=np.array(cells)
    def rotate(z):
        mid=z.mean(axis=-1)@basis
        width=((z[:,:,1]-z[:,:,0])/2)@abs(basis)+.000001
        return mid-width,mid+width
    lo,hi=rotate(arr[:,:3]);vl,vu=rotate(arr[:,3:])
    t=a.times[-2];end=a.times[-1];speed=row["speed"];force=row["extra_force"]
    hit=False;earliest=None;max_hit=0.;earth_fail=False
    hit_angle_max=0.;turn_groups=0;total_groups=0;worst=None
    target=o.unit(moon(a,end)-a.sites[-1]);score_range=norm(moon(a,end)-a.sites[-1])
    near_pad=650. # Above 45 km/s times the largest allowed terminal light time on these arcs.
    mfinal=moon(a,end);site=a.sites[-1];target=o.unit(mfinal-site)
    if shift:
        tt=end-14*3600
        x=a.propagate(np.array(load("juice-forecast.json")["state"]),tt,0)
        drift=x[3:]-(moon(a,tt+1)-moon(a,tt-1))/2.
        direction=o.unit(mfinal-site+x[:3]-moon(a,tt)+drift*(end-tt))
        target=o.unit((1-shift)*target+shift*direction)
    offset=o.distance(target,o.unit(mfinal-site))/3600
    capfar=0.;old_key=None;cut=0
    while t<end:
        h=min(300.,end-t)
        mid=t+h/2;mp=moon(a,mid)
        mv=(moon(a,mid+1)-moon(a,mid-1))/2
        el=lo-speed*h-.1;eu=hi+speed*h+.1
        ml=el-mp-h-.01;mu=eu-mp+h+.01
        er=norm(np.maximum(0,np.maximum(el,-eu)))
        if np.min(er)<6400:
            earth_fail=True
            break
        # Exclude far cells that have certainly entered the sphere at this time.
        must=norm(np.maximum(abs(lo-moon(a,t)),abs(hi-moon(a,t))))<radius-.1
        if must.any() and t>a.times[-2]:
            cut+=sum(must)
            keep=~must
            lo,hi,vl,vu,el,eu,ml,mu,er=[q[keep] for q in (lo,hi,vl,vu,el,eu,ml,mu,er)]
            if len(lo)==0:break
        l0,u0=gravity(el,eu,E,6378.)
        l1,u1=gravity(ml,mu,M,radius)
        li,ui=gravity(mp-h-.01,mp+h+.01,M,340000.)
        rad=norm(np.maximum(abs(el),abs(eu)))
        pad=2*S*rad/(1.44e8-rad)**3+12*.001083*E*6379**2/np.maximum(er,6378.)**4+force+1e-13
        if physical:pad+=envelope(radius)+envelope(340000)
        al=l0+l1+li-pad[:,None];au=u0+u1+ui+pad[:,None]
        possible=norm(np.maximum(0,np.maximum(ml,-mu)))<=radius+.01
        if possible.any():
            hit=True
            if earliest is None:earliest=t
            vel=np.maximum(abs(vl[possible]-mv),abs(vu[possible]-mv))
            vmax=min(speed+2,max(norm(vel))+h*max(norm(np.maximum(abs(al[possible]),abs(au[possible]))))+.00002*h+.001)
            # Round upward and reuse bounds for conservative 15-minute / 0.1 km/s bins.
            duration=np.ceil((end-t)/900)*900;key=(duration,np.ceil(vmax*10)/10)
            if directional:
                fallback=travel(key[0],key[1],radius,speed,force,physical,360000.)+near_pad
                fallback_deg=np.rad2deg(np.arcsin(fallback/score_range)) if fallback<score_range else 180.
                centers=(lo[possible]+hi[possible])/2
                extent=(hi[possible]-lo[possible])/2+(speed+2)*h+.1
                # Center the entry position on the Moon when the radius shell beats the Cartesian box.
                q=centers-mp; qr=norm(extent)
                use_shell=qr>radius+.1
                q[use_shell]=0;qr=np.minimum(qr,radius+.1)
                vc=(vu[possible]+vl[possible])/2-mv
                verr=norm((vu[possible]-vl[possible])/2)+h*norm(np.maximum(abs(al[possible]),abs(au[possible])))+.00002*h+.001
                vmin=np.maximum(0,norm(vc)-verr);vmaxs=np.minimum(speed+2,norm(vc)+verr)
                centers=moon(a,end)-a.sites[-1]+q+vc*(end-t-h/2)
                errors=qr+verr*(end-t)+norm(vc)*h/2+near_pad
                grouped={}
                for i,(vlow,vhigh) in enumerate(zip(vmin,vmaxs)):
                    group=(np.floor(vlow*4)/4,np.ceil(vhigh*4)/4)
                    grouped.setdefault(group,[]).append(i)
                for group,ids in grouped.items():
                    total_groups+=1
                    dv=turn_bound(key[0],*group,radius,speed,force)
                    angle=fallback_deg+offset
                    if dv is not None:
                        turn_groups+=1
                        ratios=(errors[ids]+dv*(end-t))/norm(centers[ids])
                        if max(ratios)<1:
                            angle=min(angle,max(np.rad2deg(np.arcsin(ratios))+o.distance(o.unit(centers[ids]),target)/3600))
                            angle=max(angle,np.rad2deg(np.arcsin((radius+near_pad)/score_range))+offset)
                    if angle>hit_angle_max:
                        worst=dict(togo=(end-t)/3600,group=group,dv=dv,fallback=fallback_deg,velocity_error=float(max(verr[ids])),
                                   base_offset=float(max(o.distance(o.unit(centers[ids]),target)/3600)),
                                   approach_radius=float(max(errors[ids])),distance=float(min(norm(centers[ids]))))
                        hit_angle_max=angle
            if key!=old_key:
                max_hit=max(max_hit,travel(key[0],key[1],radius,speed,force,physical,360000.))
                old_key=key
        lo=lo+h*vl+h*h*al/2-1e-7;hi=hi+h*vu+h*h*au/2+1e-7
        vl=np.maximum(-speed,vl+h*al-1e-10);vu=np.minimum(speed,vu+h*au+1e-10)
        t+=h
    hit_rho=norm(mfinal-site)
    # Global barycentric light-time pad for the prior's speed and the terminal reachable region.
    last_radius=norm(np.maximum(np.abs(row["last_state_lo"]),np.abs(row["last_state_hi"])))+6500
    pad=(speed+35)*(last_radius+speed*(end-a.times[-2]))/(o.LIGHT_SPEED-speed-35)+.1
    if hit:max_hit+=pad
    if len(lo) and not earth_fail and t>=end:
        centers=(lo+hi)/2-site;radii=norm((hi-lo)/2)+pad
        ratio=radii/norm(centers)
        if np.max(ratio)>=1:capfar=180.
        else:capfar=max(o.distance(o.unit(centers),target)/3600+np.rad2deg(np.arcsin(ratio)))
    elif earth_fail:capfar=180.
    hitangle=0 if not hit else (np.rad2deg(np.arcsin(max_hit/hit_rho)) if max_hit<hit_rho else 180.)
    hitangle+=offset
    if directional:
        assert not physical
        hitangle=min(hitangle,hit_angle_max)
    diam=min(360.,2*max(capfar,hitangle))
    return dict(worst_turn=worst,turn_groups=turn_groups,total_groups=total_groups,cells=N*N,threshold_km=radius,physical=physical,near_max_radius_km=max_hit,hit_radius_deg=hitangle,
                far_radius_deg=capfar,diameter_deg=diam,cap_center=target.tolist(),
                earliest_possible_hit_hours_before_score=(end-earliest)/3600 if hit else None,
                earth_floor_failure=earth_fail,wholly_inside_removed=int(cut))
if __name__=="__main__":
    a=Arc("juice");configs=json.loads((HERE/"global.json").read_text())
    index=int(sys.argv[1]);res=[]
    for R in (12000.,25000.,50000.):
        sol=solve(a,configs[index],12,R)
        res.append(sol);save("branches-"+str(index)+".json",res);print(index,R,sol,flush=True)
