"""Fixed-budget covariance and profile systematic-ranging comparison on two historical arcs."""
import re,html
from datetime import datetime,timezone
from time import perf_counter
from itertools import product,combinations
from scipy.optimize import least_squares
from common import *
def second_tess(arc_key="september"):
    train,hold,site=("2018 09 09","2018 09 10","G96") if arc_key=="september" else ("2018 11 17","2018 11 19","703")
    eopfile="tess-eop.txt" if arc_key=="september" else "november-eop.txt"
    a=Arc("tess")
    rows=[]
    for m in re.finditer(r'<a name="o(\d+)"></a>([^\n]+)',(MISSIONS/"tess.html").read_text()):
        s=html.unescape(re.sub("<[^>]+>","",m[2]))
        if s[-3:] not in ("703","G96") or not s[15:32].strip().startswith((train+".",hold+".")): continue
        t=forecast.parse_time(s[15:32].strip())
        ra,dec=np.deg2rad([15*o.sexagesimal(s[32:44]),o.sexagesimal(s[44:56])])
        rows.append((t,s[-3:],np.array([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)]),int(m[1]),s[64]))
    rows.sort(key=lambda r:r[0]); a.rows=[r for r in rows if datetime.fromtimestamp(r[0],timezone.utc).day==int(train[-2:])]
    score=next(r for r in rows if r[1]==site and datetime.fromtimestamp(r[0],timezone.utc).day==int(hold[-2:]))
    epoch=a.rows[0][0];a.JD=2440587.5+epoch/86400;a.env=dict(o.env,jd_utc_start=a.JD)
    a.eop=np.array([[float(s[i:j]) for i,j in ((7,15),(18,27),(37,46),(58,68))]
                  for s in (HERE/eopfile).read_text().splitlines()])
    a.times=np.array([r[0]-epoch for r in a.rows]+[score[0]-epoch])
    a.obs=np.array([r[2] for r in a.rows]);a.east,a.north=o.basis(a.obs)
    a.sites=np.array([a.site(r[1],t) for r,t in zip(a.rows+[score],a.times)])
    a.earth_v=np.array([(o.hpop_de440_positions_km(a.JD+(t+60)/86400,a.ephem)["earth"]-
                         o.hpop_de440_positions_km(a.JD+(t-60)/86400,a.ephem)["earth"])/120 for t in a.times])
    use=[i for i,r in enumerate(a.rows) if r[1]==a.rows[0][1]]
    a.angular_slope=np.polyfit(a.times[use],a.obs[use],1)[0]
    a.site_v=(a.site(a.rows[0][1],1)-a.site(a.rows[0][1],-1))/2
    a.axes=np.array([a.obs[0],a.east[0],a.north[0]])
    return a,score
def run(name,arc_key="september"):
    if name=="tess": a,score=second_tess(arc_key)
    else: a=Arc("juice")
    n=len(a.rows)
    calls=0
    limit=4096
    def eval_sky(x,n):
        nonlocal calls
        calls+=1
        assert calls<=limit
        return a.sky(x,n)
    def res(x):
        pred=eval_sky(x,n)
        # Chords projected onto observed tangent frames, with antipodal points excluded later.
        return (o.ASEC*np.column_stack((np.sum(pred*a.east,axis=1),np.sum(pred*a.north,axis=1)))).ravel()
    def stats(x):
        u=eval_sky(x,n+1)
        return max(o.distance(u[:-1],a.obs)),u[-1]
    starts=[(r,v) for r in ((2e5,3e5,4e5) if name=="tess" else (1e6,2e6,3e6))
            for v in ((-1,1) if name=="tess" else (-4,-2))]
    fits=[];tic=perf_counter()
    for r,v in starts:
        seed=a.seed(r,v)
        opt=least_squares(lambda p:res(seed+a.delta(p)),np.zeros(6),jac="3-point",max_nfev=32)
        fits.append((np.linalg.norm(opt.fun),seed+a.delta(opt.x)))
    _,x=min(fits,key=lambda item:item[0])
    basemax,center=stats(x);common_calls=calls;common_time=perf_counter()-tic
    tic=perf_counter();eps=1e-3
    jac=np.column_stack([(res(x+a.delta(e))-res(x-a.delta(e)))/(2*eps) for e in eps*np.eye(6)])
    east,north=o.basis(center)
    def tangent(x):
        u=eval_sky(x,n+1)[-1];return np.array([u@east,u@north])*180/np.pi
    J=np.column_stack([(tangent(x+a.delta(e))-tangent(x-a.delta(e)))/(2*eps) for e in eps*np.eye(6)])
    # SVD inverse for unit standard deviation in each tangent coordinate.
    pinv=J@np.linalg.pinv(jac)
    axes=np.linalg.svd(pinv,compute_uv=False)
    covcalls=calls;cput=perf_counter()-tic
    tic=perf_counter();calls=common_calls
    # Profile over two fixed radial coordinates, fitting transverse position and velocity.
    relative=x[:3]-a.sites[0]; r0=relative@a.obs[0];v0=(x[3:]-a.site_v)@a.obs[0]
    if name=="juice":
        coarse=(np.linspace(0.5e6,3e6,7),np.linspace(-6,0,7))
        fine=(r0+np.linspace(-4000,4000,13),v0+np.linspace(-.024,.024,13))
    else:
        coarse=(np.linspace(200000,500000,7),np.linspace(-2,2,7))
        fine=(r0+np.linspace(-60000,60000,13),v0+np.linspace(-.45,.45,13))
    failed=0
    idx=[1,2,4,5];dirs={2:[],3:[]}
    for r,v in list(product(*coarse))+list(product(*fine)):
        y=x.copy();y[:3]+=(r-r0)*a.obs[0];y[3:]+=(v-v0)*a.obs[0]
        for step in range(2):
            base=res(y)
            D=np.column_stack([(res(y+a.delta(np.eye(6)[i]*.001))-base)/.001 for i in idx])
            p=np.zeros(6);p[idx]=np.linalg.lstsq(D,-base,rcond=1e-12)[0];y+=a.delta(p)
        if np.max(np.linalg.norm(res(y).reshape(-1,2),axis=1))>3.01:continue
        try:
            val,u=stats(y)
        except ValueError:
            failed+=1
            continue
        if val>3: continue
        geom=a.geometry(y)
        if not a.physical(geom):continue
        if name=="juice":
            enc=a.encounter(y)
            if not enc["resolved"] or enc["moon_radius_km"]<1737.4:continue
        moon=o.hpop_de440_positions_km(a.JD+a.times[-1]/86400,a.ephem)["moon"]-a.sites[-1]
        if u@a.sites[-1]<=0 or o.distance(u,o.unit(moon))/3600<=np.rad2deg(np.arcsin(1737.4/np.linalg.norm(moon))):continue
        for cap in dirs:
            if val<=cap: dirs[cap].append((u,y,val))
    out=dict(rows=[r[3] for r in a.rows],flags=[r[4] for r in a.rows],common_calls=common_calls,common_seconds=common_time,
             sigma1_deg=axes.tolist(),gaussian95_major_diameter_deg=float(2*np.sqrt(5.991465)*axes[0]),
             baseline_max_as=basemax,baseline_state=x.tolist(),center=center.tolist(),covariance_calls=covcalls,
             covariance_seconds=cput+common_time,ranging_calls=calls,ranging_seconds=perf_counter()-tic+common_time,max_calls=limit,
             failed_grid_forecasts=failed,grid_cells=218,grid={})
    for cap,items in dirs.items():
        pairs=list(combinations(items,2))
        if pairs:
            pair=max(pairs,key=lambda z:o.distance(z[0][0],z[1][0]))
            out["grid"][cap]=dict(count=len(items),lower_deg=o.distance(pair[0][0],pair[1][0])/3600,
               witnesses=[dict(state=y.tolist(),max_as=v,los=u.tolist()) for u,y,v in pair])
        else:out["grid"][cap]=dict(count=len(items))
    prefix=name if arc_key=="september" else name+"-"+arc_key
    if arc_key!="september":save(prefix+"-unscored.json",out)
    if name=="tess":
        out.update(score_source_row=score[3],holdout_error_deg=o.distance(center,score[2])/3600,
                   fit_geometry=a.geometry(x),epoch_jd=a.JD,forecast_jd=a.JD+a.times[-1]/86400)
    save(prefix+"-comparison.json",out);print(name,out["common_calls"],out["ranging_calls"],out["baseline_max_as"],out.get("holdout_error_deg"),[(k,v.get("lower_deg")) for k,v in out["grid"].items()],flush=True)
if __name__=="__main__":
    run(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else "september")
