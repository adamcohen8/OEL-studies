"""Optical forecasts using OEL, with each scoring direction excluded from the fit."""
from pathlib import Path
from datetime import datetime, timezone
import html
import json
import re
import sys
from itertools import product
import numpy as np
from scipy.optimize import least_squares, minimize, minimize_scalar
from scipy.spatial.distance import pdist, squareform

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent / "orion"))
import orion as o

def parse_time(s):
    if s.startswith("K"):
        digits=s[8:].strip()
        clock=datetime.strptime(s[:8]+digits[:6],"K%y%m%d:%H%M%S").replace(tzinfo=timezone.utc)
        return clock.timestamp()+float("0."+(digits[6:] or "0"))
    year,month,day=s.split()
    return datetime(int(year),int(month),int(float(day)),tzinfo=timezone.utc).timestamp()+(float(day)%1)*86400

def read_rows(target):
    records=[]
    for match in re.finditer(r'<a name="o(\d+)"></a>([^\n]+)',(HERE/(target+".html")).read_text()):
        line=html.unescape(re.sub("<[^>]+>","",match[2]))
        stamp=line[15:32].strip()
        if target=="tess":
            if line[-3:] not in ("703","G96") or not stamp.startswith(("2018 09 08","2018 09 09")):
                continue
        elif not re.match(r"2024 08 1[5-9]|K24081[5-9]",stamp):
            continue
        timestamp=parse_time(stamp)
        ra,dec=np.deg2rad([15*o.sexagesimal(line[32:44]),o.sexagesimal(line[44:56])])
        records.append((timestamp,line[-3:],np.array([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)]),
                        int(match[1]),line[64:65]))
    return sorted(records,key=lambda row:row[0])

class Arc:
    def __init__(self,target):
        records=read_rows(target)
        if target=="tess":
            self.rows=[r for r in records if datetime.fromtimestamp(r[0],timezone.utc).day==8]
            prediction=next(r[0] for r in records if r[1]=="703" and datetime.fromtimestamp(r[0],timezone.utc).day==9)
        else:
            self.rows=records
            prediction=datetime(2024,8,20,6,tzinfo=timezone.utc).timestamp()
        self.target=target
        self.JD=2440587.5+self.rows[0][0]/86400
        self.times=np.array([r[0]-self.rows[0][0] for r in self.rows]+[prediction-self.rows[0][0]])
        assert min(np.diff(self.times))>0
        self.obs=np.array([r[2] for r in self.rows])
        self.east,self.north=o.basis(self.obs)
        self.env=dict(o.env,jd_utc_start=self.JD)
        self.engine=o.OrbitPropagator(integrator="rkf78",numeric_backend="rust",
            plugins=o.engine.plugins,adaptive_atol=1e-11,adaptive_rtol=1e-12)
        self.ephem=o.ephem
        self.fixed={}
        for line in (HERE/"stations.txt").read_text().splitlines():
            lon=np.deg2rad(float(line[4:13]))
            self.fixed[line[:3]]=6378.14*np.array([np.cos(lon)*float(line[13:21]),np.sin(lon)*float(line[13:21]),float(line[21:30])])
        self.eop=np.array([[float(line[a:b]) for a,b in ((7,15),(18,27),(37,46),(58,68))]
                            for line in (HERE/"eop.txt").read_text().splitlines()])
        self.sites=np.array([self.site(code,t) for code,t in zip([r[1] for r in self.rows]+["703"],self.times)])
        self.earth_v=np.array([(o.hpop_de440_positions_km(self.JD+(t+60)/86400,self.ephem)["earth"]-
                                 o.hpop_de440_positions_km(self.JD+(t-60)/86400,self.ephem)["earth"])/120 for t in self.times])
        code=self.rows[0][1]
        idx=[i for i,r in enumerate(self.rows) if r[1]==code and self.times[i]<1800]
        self.angular_slope=np.polyfit(self.times[idx],self.obs[idx],1)[0]
        self.site_v=(self.site(code,1)-self.site(code,-1))/2
        self.axes=np.array([self.obs[0],self.east[0],self.north[0]])
        self.scale=np.array([1e5,10,10,1,.01,.01])
    def site(self,code,t):
        mjd=self.JD+t/86400-2400000.5
        assert self.eop[0,0]<=mjd<=self.eop[-1,0]
        xp,yp,dut=[np.interp(mjd,self.eop[:,0],s) for s in self.eop[:,1:].T]
        return o.eci_to_ecef_rotation_hpop_like(t,self.JD,dut1_s=dut,xp_arcsec=xp,yp_arcsec=yp,dat_s=37).T@self.fixed[code]
    def seed(self,r,rr):
        vel=self.site_v+rr*self.obs[0]+r*self.angular_slope
        return np.r_[self.sites[0]+r*self.obs[0]+(self.earth_v[0]+vel)*(r/o.LIGHT_SPEED),vel]
    def delta(self,p):
        q=p*self.scale
        return np.r_[q[:3]@self.axes,q[3:]@self.axes]
    def propagate(self,x,dt,t):
        return self.engine.propagate(x,dt,t,o.zero,self.env,o.ctx)
    def sky(self,x,n):
        states=[]
        prev=0.
        for t in self.times[:n]:
            if t>prev:
                x=self.propagate(x,t-prev,prev)
            states.append(x)
            prev=t
        states=np.array(states)
        rel=states[:,:3]-self.sites[:n]
        vel=states[:,3:]+self.earth_v[:n]
        ray=rel
        for _ in range(3):
            ray=rel-vel*np.linalg.norm(ray,axis=1)[:,None]/o.LIGHT_SPEED
        return o.unit(ray)
    def residual(self,x):
        pred=self.sky(x,len(self.obs))
        return o.ASEC*np.column_stack((np.sum(pred*self.east,axis=1),np.sum(pred*self.north,axis=1)))/np.sum(pred*self.obs,axis=1)[:,None]
    def geometry(self,x):
        finish=self.times[-1]
        # Sample the full coast, then minimize lunar distance around the closest sampled epoch.
        grid=np.unique(np.r_[np.arange(0,finish,600.),finish])
        prev=0.
        lunar=[]
        geo=[]
        speed=[]
        for t in grid:
            if t>prev:
                x=self.propagate(x,t-prev,prev)
            lunar.append(np.linalg.norm(x[:3]-o.hpop_de440_positions_km(self.JD+t/86400,self.ephem)["moon"]))
            geo.append(np.linalg.norm(x[:3]))
            speed.append(np.linalg.norm(x[3:]))
            prev=t
        return dict(min_earth_km=min(geo),max_earth_km=max(geo),min_moon_sample_km=min(lunar),max_speed_km_s=max(speed))
    def physical(self,g):
        if self.target=="tess":
            return g["min_earth_km"]>201000 and g["max_earth_km"]<499000 and g["min_moon_sample_km"]>101500 and g["max_speed_km_s"]<2.99
        return g["min_earth_km"]>6500
    def encounter(self,state):
        # Dense high-accuracy local search centered on ESA's published encounter date,
        # with a broad bracket located from a sampled minimum.
        def point(t):
            return self.propagate(state,t,0)
        def lunar(t):
            return np.linalg.norm(point(t)[:3]-o.hpop_de440_positions_km(self.JD+t/86400,self.ephem)["moon"])
        start=(datetime(2024,8,19,12,tzinfo=timezone.utc).timestamp()-self.rows[0][0])
        ts=np.linspace(start,self.times[-1],49)
        k=int(np.argmin([lunar(t) for t in ts]))
        if not 0<k<len(ts)-1: return dict(resolved=False)
        fit=minimize_scalar(lunar,bounds=(ts[k-1],ts[k+1]),method="bounded")
        return dict(resolved=True,moon_radius_km=fit.fun,utc=datetime.fromtimestamp(self.rows[0][0]+fit.x,timezone.utc).isoformat())

def forecast(target):
    arc=Arc(target)
    n=len(arc.obs)
    starts=product((225000,325000,425000),(-1,0,1)) if target=="tess" else product((1.5e6,2.5e6,3.5e6),(-5,-3))
    sols=[]
    for r,rr in starts:
        seed=arc.seed(r,rr)
        def res(p):
            return arc.residual(seed+arc.delta(p)).ravel()
        fit=least_squares(res,np.zeros(6),jac="3-point",max_nfev=65,
                          bounds=([-r*.75/1e5,-1000,-1000,-5,-100,-100],[r/1e5,1000,1000,5,100,100]))
        sols.append((np.linalg.norm(fit.fun),seed+arc.delta(fit.x)))
    best=min(sols,key=lambda pair:pair[0])[1]
    def stats(x):
        v=o.distance(arc.sky(x,n),arc.obs)
        return dict(rms_as=float(np.sqrt(np.mean(v**2))),max_as=max(v))
    baseline=stats(best)
    # Find an interior bounded-error fit even if least squares violates the cap.
    def prof(p):
        return arc.residual(best+arc.delta(p[:6]))
    sol=minimize(lambda p:p[6],np.r_[np.zeros(6),baseline["max_as"]+.01],
                 constraints=[dict(type="ineq",fun=lambda p:p[6]-np.linalg.norm(prof(p),axis=1))],
                 bounds=[(-4,4)]+[(-400,400)]*2+[(-5,5)]+[(-100,100)]*2+[(0,100)],
                 method="SLSQP",options=dict(maxiter=90,ftol=1e-8,eps=3e-6))
    interior=best+arc.delta(sol.x[:6])
    center=best
    u=arc.sky(center,n+1)[-1]
    geometry=arc.geometry(center)
    output=dict(target=target,rows=[r[3] for r in arc.rows],flags=[r[4] for r in arc.rows],
                ts=arc.times.tolist(),stations=[r[1] for r in arc.rows],state=best.tolist(),baseline=baseline,
                interior_state=interior.tolist(),interior=stats(interior),geometry=geometry,los=u.tolist(),
                test_site_alt_deg=float(np.rad2deg(np.arcsin(o.unit(arc.sites[-1])@u))))
    candidates=[]
    if stats(interior)["max_as"]<2:
        for x in [best,interior]+[x for _,x in sols]:
            if stats(x)["max_as"]<2: candidates.append(x)
        jac=least_squares(lambda p:arc.residual(interior+arc.delta(p)).ravel(),np.zeros(6),jac="3-point",max_nfev=1).jac
        _,sing,vt=np.linalg.svd(jac,full_matrices=False)
        for v,s in zip(vt[-2:],sing[-2:]):
            for sign in (-1,1):
                def slack(a):
                    return stats(interior+arc.delta(v*sign*a))["max_as"]-1.98
                hi=min(1/max(s,1e-9),.1)
                while hi<6 and slack(hi)<0: hi=min(2*hi,6.)
                if slack(0)<0:
                    a=o.brentq(slack,0,hi) if slack(hi)>0 else hi
                    candidates.append(interior+arc.delta(v*sign*a))
    # Filter duplicate candidates before the physical checks.
    saved=[]
    for x in candidates:
        if any(np.linalg.norm(x[:3]-np.array(v["state"][:3]))<.005 and np.linalg.norm(x[3:]-np.array(v["state"][3:]))<1e-7 for v in saved):
            continue
        g=arc.geometry(x)
        if not arc.physical(g): continue
        direction=arc.sky(x,n+1)[-1]
        alt=np.rad2deg(np.arcsin(o.unit(arc.sites[-1])@direction))
        if alt<0: continue
        record=dict(state=x.tolist(),los=direction.tolist(),**stats(x),**g,alt_deg=alt)
        if target=="juice":
            record["encounter"]=arc.encounter(x)
            if not record["encounter"]["resolved"] or record["encounter"]["moon_radius_km"]<1737.4: continue
        saved.append(record)
    if len(saved)>=2:
        dirs=np.array([v["los"] for v in saved])
        pair=np.unravel_index(np.argmax(squareform(pdist(dirs))), (len(saved),len(saved)))
        output.update(witnesses=[saved[i] for i in pair],lower_diameter_deg=o.distance(*dirs[list(pair)])/3600)
    if target=="juice": output["encounter"]=arc.encounter(best)
    output["feasible_visible_count"]=len(saved)
    (HERE/(target+"-forecast.json")).write_text(json.dumps(output,indent=2)+"\n")
    print(target,baseline,output["interior"],geometry,output.get("lower_diameter_deg"),flush=True)

if __name__=="__main__":
    forecast(sys.argv[1])
