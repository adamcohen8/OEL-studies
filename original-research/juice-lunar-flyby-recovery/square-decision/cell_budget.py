"""Measure the margin of the available continuous error-ball verifier at two limiting orbits."""
from explore import *
from tubes import prepare,enclosure
from time import perf_counter
from copy import copy
pth=D/'nominal-certificate.json'
w=json.loads(pth.read_text())['2']
q=w['geometry'];c=np.array(q['center']);e,n=o.basis(c);roll=np.deg2rad(q['roll_deg'])
axes=np.array([e*np.cos(roll)+n*np.sin(roll),-e*np.sin(roll)+n*np.cos(roll)])
face=np.vstack([axes,-axes])-np.tan(np.deg2rad(2.2))*c
face=face/np.linalg.norm(face,axis=1)[:,None]
from outer_test import path
segment=copy(a);segment.times=a.times[-2:]-a.times[-2];segment.sites=a.sites[-2:]
# prepare/enclosure use absolute epochs relative to a.JD. The propagator's environment
# epoch must move with it when starting at last reception.
segment.JD=a.JD+a.times[-2]/86400;segment.env=dict(a.env,jd_utc_start=segment.JD)
radii=[.001,.01,.03,.1,.3,1.,3.,10.,30.,100.]
out=[]
for j in (np.argmin([q['validation']['min_moon_floor_km'] for q in w['accepted']]),np.argmax([q['validation']['min_moon_floor_km'] for q in w['accepted']])):
    row=w['accepted'][j];last=path(np.array(row['state']))[6*len(a.obs)-6:6*len(a.obs)]
    start=perf_counter();knots=prepare(segment,last,.25);prep=perf_counter()-start;results=[]
    for R in radii:
        t=perf_counter();q=enclosure(knots,segment,R,0.)
        item=dict(radius_km=R,velocity_radius_km_s=R/100000,seconds=perf_counter()-t,status=q['status'])
        if q['status']=='bounded':
            center=np.array(q['center']);clear=np.rad2deg(np.arcsin(-max(face@center)))
            item.update(sky_ball_deg=q['field_radius_deg'],face_clearance_deg=clear,contained=bool(q['field_radius_deg']<clear))
        else:item.update(reason=q['reason'])
        results.append(item)
    out.append(dict(witness=int(j),reference_seconds=prep,steps=len(knots)-1,results=results))
save('cell-budget.json',out)
for v in out:print(v['witness'],[(q['radius_km'],q.get('contained')) for q in v['results']],flush=True)
