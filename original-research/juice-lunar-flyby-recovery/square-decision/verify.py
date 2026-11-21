"""Verify limiting OEL trajectories and convert them to an all-orientations side bound."""
from footprint import *
from certify_pairs import certify
from scipy.spatial.distance import pdist,squareform
def controls(pts):
    u=np.array([q['los'] for q in pts]);shape=fit(u);c=np.array(shape['center']);e,n=o.basis(c)
    roll=np.deg2rad(shape['roll_deg']);e,n=e*np.cos(roll)+n*np.sin(roll),-e*np.sin(roll)+n*np.cos(roll)
    # Corners of a detector-aligned bounding box and the farthest pair constrain
    # translation and roll. Include both endpoints of six nearby orientations.
    ids=set()
    for angle in np.deg2rad([-5,0,5,85,90,95,35,45]):
        z=(u@(e*np.cos(angle)+n*np.sin(angle)))/(u@c)
        ids.update([int(z.argmin()),int(z.argmax())])
    i,j=np.unravel_index(np.argmax(squareform(pdist(u))),(len(u),len(u)));ids.update([i,j])
    return [pts[k] for k in sorted(ids)]
def process(caps,savefile):
    results={};cache={}
    if (D/savefile).exists():results=json.loads((D/savefile).read_text())
    for cap,pts in caps.items():
        if cap in results:continue
        chosen=controls(pts);done=[];rejected=[]
        for q in chosen:
            key=tuple(q['state'])
            if key not in cache:cache[key]=certify(a,np.array(q['state']),float(cap))
            c=cache[key]
            if c.get('status')!='bounded':
                from types import FunctionType
                from tubes import prepare
                tighter=FunctionType(certify.__code__,dict(certify.__globals__,prepare=lambda a,x,scale,knots:prepare(a,x,scale/3,knots)),argdefs=certify.__defaults__)
                c=tighter(a,np.array(q['state']),float(cap))
                cache[key]=c
            if c.get('status')=='bounded' and c['training_max_upper_as']<=float(cap) and c['min_moon_floor_km']>1737.4 and c['min_earth_floor_km']>200000 and c['training_moon_floor_km']>1e5 and c['max_speed_upper']<6:
                u=np.array(c['forecast_center']);m=a.times[-1]
                lunar=o.hpop_de440_positions_km(a.JD+m/86400,a.ephem)['moon']-a.sites[-1]
                g=o.distance(u,o.unit(lunar))/3600-np.rad2deg(np.arcsin(1737.4/np.linalg.norm(lunar)))
                assert u@a.sites[-1]>0 and g>c['forecast_radius_as']/3600+.05
                done.append(dict(state=q['state'],validation=c))
            else:rejected.append(c)
        if len(done)<2:continue
        u=np.array([q['validation']['forecast_center'] for q in done]);rads=[q['validation']['forecast_radius_as'] for q in done]
        shape=fit(u,rads)
        pair=np.unravel_index(np.argmax(squareform(pdist(u))),(len(u),len(u)))
        pair_bound=o.distance(*u[list(pair)])/3600-sum(rads[i] for i in pair)/3600
        # Check explicit inclusion of the computed emission balls in the fitted tangent-plane square.
        center=np.array(shape['center']);e,n=o.basis(center);roll=np.deg2rad(shape['roll_deg'])
        axes=np.array([e*np.cos(roll)+n*np.sin(roll),-e*np.sin(roll)+n*np.cos(roll)])
        side=shape['fit_side_deg']+.01;face=np.vstack([axes,-axes])-np.tan(np.deg2rad(side/2))*center
        confined=bool(np.max(u@face.T+np.deg2rad(np.asarray(rads)/3600)[:,None]*np.linalg.norm(face,axis=1))<0)
        results[cap]=dict(accepted=done,rejections=rejected,minimum_pair_deg=pair_bound,geometry=shape,finite_enclosing_side_test_deg=side,finite_enclosing_passed=confined)
        save(savefile,results);print(cap,len(done),shape['lower_test_side_deg'],pair_bound,flush=True)
if __name__=='__main__':
    if sys.argv[1]=='nominal':
        pts=sum([json.loads((D/f).read_text())['points'] for f in ['exploration.json','sky-support.json']],[])
        process({'2':pts},'nominal-certificate.json')
    else:
        data=json.loads((D/'transition-search.json').read_text())
        process({k:v['extrema'] for k,v in data.items() if float(k)>2},'transition-certificate.json')
