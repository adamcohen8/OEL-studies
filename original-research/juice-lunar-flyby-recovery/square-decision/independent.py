"""Replay the small exclusion and the tight nominal witnesses with the alternate OEL integrator."""
from explore import *
import ast
import forecast as flight
from itertools import combinations
fc=D.parent/'field-certificate'
src=(fc.parent/'followup/corrected/missions/emission.py').read_text()
fn=next(z for z in ast.parse(src).body if isinstance(z,ast.FunctionDef) and z.name=='validate')
ns=dict(f=flight,o=o,np=np);exec(compile(ast.Module(body=[fn],type_ignores=[]),'historical-emission-replay','exec'),ns)
data=json.loads((D/'all-rolls.json').read_text());out={}
for cap in ['2','2.15']:
    group=data[cap];rows=json.loads((D/(group['source']+'.json')).read_text())[cap]['accepted']
    select=group.get('indices',list(range(len(rows))))
    audits={}
    for k in range(0,len(select),2):
        ids=select[k:k+2]
        if len(ids)==1:ids=ids+[select[0]]
        rr=ns['validate']('juice',[rows[i] for i in ids])
        for i,r in zip(ids,rr['witnesses']):
            q=rows[i]['validation'];r['forecast_difference_as']=o.distance(np.array(r['score_los']),np.array(q['forecast_center']))
            assert min(r['min_limb_clearance_deg'],r['min_elevation_deg'])>0 and r['max_fit_as']<float(cap)
            assert r['forecast_difference_as']<q['forecast_radius_as']
            audits[i]=r
    out[cap]=audits;save('independent.json',out)
    print(cap,len(audits),max(r['max_fit_as'] for r in audits.values()),flush=True)
