from transition import *
from verify import process
from geometry_certificate import summarize
if __name__=='__main__':
    cap=2.15;sk.cap=cap-.007;sk.floor=1740.5
    initial=json.loads((D/'transition-search.json').read_text())[str(cap)]['extrema']
    pts=list(initial)
    for deg in (79.2,169.2,259.2,349.2):
        r=np.deg2rad(deg);axis=sk.east*np.cos(r)+sk.north*np.sin(r)
        q=sk.seek(r,max(pts,key=lambda p:np.array(p['los'])@axis))
        if q:pts.append(q)
    save('near-threshold.json',dict(points=pts))
    process({str(cap):pts},'threshold-certificate.json')
    out=summarize(json.loads((D/'threshold-certificate.json').read_text())[str(cap)]['accepted'])
    print(out,flush=True)
