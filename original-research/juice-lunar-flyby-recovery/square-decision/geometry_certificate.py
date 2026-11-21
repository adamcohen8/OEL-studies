"""Prove a detector exclusion for a finite set with angular errors using all roll angles."""
from itertools import combinations
from math import asin, atan, cos, sin, sqrt, tan, pi
from pathlib import Path
import json
import numpy as np
D=Path(__file__).resolve().parent
def decision(points,radii,side=4.4,samples=16384):
    points=np.array(points);radii=np.deg2rad(np.asarray(radii)/3600)
    assert len(points)>1 and np.max(abs(np.linalg.norm(points,axis=1)-1))<1e-12
    ids=np.array(list(combinations(range(len(points)),2)));i,j=ids.T
    diff=points[i]-points[j];k=np.argmax(np.linalg.norm(diff,axis=1))
    endpoints=points[i[k]]+points[j[k]]
    mid=endpoints/np.linalg.norm(endpoints)
    east=np.cross([0.,0.,1.],mid);east/=np.linalg.norm(east);north=np.cross(mid,east)
    limit=np.deg2rad(side/2)
    circ=atan(sqrt(2)*tan(limit))
    co=(cos(circ+radii[i[k]])+cos(circ+radii[j[k]]))/np.linalg.norm(endpoints)-1e-13
    if co>1:
        return dict(excluded=True,reason='pair exceeds spherical detector diagonal')
    delta=np.arccos(np.clip(co,-1,1))+1e-11
    dx,dy,dz=diff@east,diff@north,diff@mid
    tang=np.hypot(dx,dy);errors=(1-cos(delta))*tang+sin(delta)*abs(dz)+radii[i]+radii[j]
    assert delta<pi/2
    roll=np.arange(samples)*pi/(2*samples);best=2.
    for rr in np.array_split(roll,max(1,samples//2048)):
        C,S=np.cos(rr),np.sin(rr)
        w1=abs(dx[:,None]*C+dy[:,None]*S)-errors[:,None]
        w2=abs(-dx[:,None]*S+dy[:,None]*C)-errors[:,None]
        best=min(best,np.min(np.max(np.maximum(w1,w2),axis=0)))
    surplus=best-max(tang)*pi/(4*samples)-2*sin(limit)-2e-12
    return dict(excluded=bool(surplus>0),width_surplus=float(surplus),max_center_tilt_deg=float(np.rad2deg(delta)),
                tested_side_deg=side,roll_intervals=samples)
def summarize(rows):
    u=[q['validation']['forecast_center'] for q in rows];r=[q['validation']['forecast_radius_as'] for q in rows]
    ids=list(range(len(rows)))
    def check(ids,side=4.4,samples=16384):return decision([u[i] for i in ids],[r[i] for i in ids],side,samples)
    result=check(ids)
    lo,hi=3.,6.
    for _ in range(16):
        m=(lo+hi)/2
        if check(ids,m,8192)['excluded']:lo=m
        else:hi=m
    ans=dict(required_side_lower_deg=lo,all=result)
    if result['excluded']:
        # Greedy removal gives a small reproducible certificate.
        for i in ids[:]:
            subset=[j for j in ids if j!=i]
            if len(subset)>1 and check(subset)['excluded']:ids=subset
        q=np.array(u)[ids]
        max_pair=max(np.rad2deg(np.arctan2(np.linalg.norm(np.cross(x,y)),x@y))+(r[a]+r[b])/3600 for (x,y),(a,b) in zip(combinations(q,2),combinations(ids,2)))
        ans.update(indices=ids,minimal_by_single_deletion=check(ids),all_pair_distances_upper_deg=float(max_pair),
            max_training_error_as=max(rows[i]['validation']['training_max_upper_as'] for i in ids))
    return ans
if __name__=='__main__':
    out={}
    for name in ('nominal-certificate','threshold-certificate','transition-certificate'):
        p=D/(name+'.json')
        if p.is_file():
            for k,v in json.loads(p.read_text()).items():
                if k in out:continue
                out[k]=dict(source=name,**summarize(v['accepted']))
    (D/'all-rolls.json').write_text(json.dumps(out,indent=2)+'\n')
    print({k:(v['required_side_lower_deg'],v.get('indices')) for k,v in out.items()})
