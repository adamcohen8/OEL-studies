"""Approach the detector exclusion threshold from the large-error side."""
import sky_support as sk
from explore import *
from footprint import fit
from scipy.spatial.distance import pdist,squareform
from certify_pairs import certify
def compress(pts):
    u=np.array([q['los'] for q in pts]);i,j=np.unravel_index(np.argmax(squareform(pdist(u))), (len(u),len(u)))
    mid=o.unit(u[i]+u[j]);e,n=o.basis(mid);ids={i,j}
    for angle in np.linspace(0,2*np.pi,128,endpoint=False):
        ids.add(int(np.argmax(u@(e*np.cos(angle)+n*np.sin(angle)))))
    return [pts[i] for i in sorted(ids)]
if __name__=='__main__':
    start=json.loads((D/'sky-support.json').read_text())['points'];history={}
    for cap in (2.,2.2,2.4,2.1,2.15):
        sk.cap=cap-.012
        pts=list(start)
        for deg in (79.25,169.25,259.25,349.25,35,215):
            r=np.deg2rad(deg);axis=sk.east*np.cos(r)+sk.north*np.sin(r)
            seed=max(pts,key=lambda q:np.array(q['los'])@axis)
            row=sk.seek(r,seed)
            if row:pts.append(row)
        core=compress(pts);result=fit([r['los'] for r in core])
        history[str(cap)]=dict(extrema=core,fit=result)
        save('transition-search.json',history)
        print(cap,result,flush=True)
