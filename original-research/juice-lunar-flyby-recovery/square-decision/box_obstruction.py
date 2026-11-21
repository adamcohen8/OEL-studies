"""Find surviving trajectories inside each product box that witness loss of optical correlations."""
from refined_test import *
from geometry_certificate import decision
from certify_pairs import certify
from itertools import product
def examine(label,origin,J,bb):
    pts=[];seeds=[x for x in product((.02,.98),repeat=6)]
    for w in json.loads((D/'sky-support.json').read_text())['points'][:12]:
        zz=path(np.array(w['state']))[6*n-6:6*n];q=np.linalg.solve(J,zz-origin)
        unit=(q-bb[:,0])/np.ptp(bb,axis=1)
        if min(unit)>0 and max(unit)<1:seeds.append(unit.tolist())
    nominal=ref[6*n-6:6*n]
    for unit in seeds:
        q=bb[:,0]+np.array(unit)*np.ptp(bb,axis=1);z=origin+J@q
        rp=closest(z)
        if rp<1742:
            left,right=0.,1.
            for _ in range(17):
                f=(left+right)/2
                if closest(nominal+f*(z-nominal))>1743:left=f
                else:right=f
            z=nominal+left*(z-nominal);rp=closest(z)
            unit=(np.linalg.solve(J,z-origin)-bb[:,0])/np.ptp(bb,axis=1)

        obj=dict(last=z.tolist(),los=sky_last(z).tolist(),rp_km=rp,normalized_box_coordinates=list(unit))
        # First collect all cheap rays. Only repropagate decisive cases from training.
        pts.append(obj)
    pair=max(combinations(pts,2),key=lambda x:o.distance(*[np.array(w['los']) for w in x]))
    scan=o.distance(*[np.array(w['los']) for w in pair])/3600
    rec=[]
    for p in pair:
        z=np.array(p['last']);q=reverse(z);x=state(q);history=path(x)[6*n-6:6*n]
        err=np.max(abs(np.linalg.solve(J,history-origin)-(bb[:,0]+np.array(p['normalized_box_coordinates'])*np.ptp(bb,axis=1)))/np.ptp(bb,axis=1))
        verified=certify(a,x,2.)
        rec.append(dict(**p,state=x.tolist(),relative_box_mismatch=err,actual_max_error_as=maxerr(q),emission=verified))
    if all(r['emission'].get('status')=='bounded' for r in rec):
        resolved=(o.distance(*[np.array(r['emission']['forecast_center']) for r in rec])-sum(r['emission']['forecast_radius_as'] for r in rec))/3600
    else:resolved=None
    return dict(selection_count=len(pts),screened=len(seeds),pair_deg=scan,verified_separation_deg=resolved,rows=rec)
if __name__=='__main__':
    tab=json.loads((D/'outer-replay.json').read_text())
    b=np.array(json.loads((fc/'global.json').read_text())[0]['basis'])
    cfg=json.loads((D/'chart.json').read_text())
    out={}
    for label,origin,J,box in [('tangent',np.zeros(6),block_diag(b.T,b.T),np.array(tab['tightened_box'])),('correlated_chart',np.array(cfg['origin']),np.linalg.inv(cfg['transform']),np.array(cfg['box']))]:
        out[label]=examine(label,origin,J,box);save('box-obstruction.json',out)
        print(label,out[label]['verified_separation_deg'],flush=True)
