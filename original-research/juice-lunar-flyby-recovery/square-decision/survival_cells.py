"""Joint optical and lunar-survival constraints on a subcell, including every detector face."""
from slab_test import *
from refined_test import closest
from collections import Counter
fullM,fullb=setup_full();limits=np.array(cfg['box']);origin=np.array(cfg['origin'])
G=np.zeros((6,N));G[:,6*n-6:6*n]=transform*scale[6*n-6:6*n]
gs=np.abs(G).max(axis=1);GG=G/gs[:,None]
lb,ub=((bounds-ref[:,None])/scale[:,None]).T
def cell_program(box):
    return vstack([fullM,GG,-GG],format='csr'),np.r_[fullb,box[:,1]/gs,-box[:,0]/gs]
def tightbox(box):
    m,b=cell_program(box);out=box.copy()
    for k,goal in enumerate(GG):
        hi=support(m,b,lb,ub,goal);lo=support(m,b,lb,ub,-goal)
        if lo is None or hi is None:return box
        out[k,0]=max(out[k,0],-lo*gs[k]-.000001);out[k,1]=min(out[k,1],hi*gs[k]+.000001)
    return out
def sc_module():
    return sys.modules[__name__]
def evaluate(box):
    mid=box.mean(axis=1);width=np.ptp(box,axis=1)/2
    last=origin+inverse@mid;rp=closest(last)
    if rp<80 and not getattr(sc_module(),'allow_deep_reference',False):return dict(status='unresolved',reason='reference near lunar center',ref_rp=rp)
    start=perf_counter();bundle=track(last);t,x,acc,d,j,ray,light=bundle
    per=t[d[:,1].argmin()]
    ids=np.unique(np.clip(np.searchsorted(t,per+np.array([-800,-480,-240,-120,-60,0,60,120,240,480,800])),0,len(t)-2))
    sampledP=np.zeros((len(ids),3,6));sampledR=np.full(len(ids),-1.)
    R,P=tube(t,x,acc,d,j,inverse@np.diag(width),ids,sampledP,sampledR)
    if R<0:return dict(status='collision',ref_rp=rp)
    m,b=cell_program(box)
    for k,PP,RR in zip(ids,sampledP,sampledR):
        if RR<0:continue
        moon=o.hpop_de440_positions_km(a.JD+(a0+t[k])/86400,a.ephem)['moon']
        delta=x[k,:3]-moon;length=np.linalg.norm(delta);unit=delta/length
        error=RR+.021
        span=maxball(PP)+error
        cut=(1737.4**2-length**2-span**2)/(2*length)-error-.00001
        if cut<=-np.sum(abs(unit@PP)):continue
        grad=(unit@PP)/width;goal=grad@G
        high=support(m,b,lb,ub,goal)
        if high is not None and high-grad@mid<cut-.0001:
            return dict(status='collision',ref_rp=rp,support_gap_km=cut-(high-grad@mid))
        size=max(1,np.max(abs(goal)));m=vstack([m,-goal/size],format='csr');b=np.r_[b,(-cut-grad@mid+.00001)/size]
    if R>1e8:return dict(status='unresolved',ref_rp=rp,reason='remainder')
    delay=(abs(light-np.linalg.norm(ray))+maxball(P)+R+.021)/(o.LIGHT_SPEED-50)
    tests=[]
    for face in ac.faces:
        grad=(face@P[:,:6])/width;goal=grad@G
        qbound=0.
        if P.shape[1]>6:
            from quadratic import upper_polynomial
            qbound=upper_polynomial(face@P[:,6:])
        high=support(m,b,lb,ub,goal)
        excess=face@ray+R+.021+50*delay+qbound
        excess+=min(sum(abs(face@P[:,:6])),high-grad@mid+.00001) if high is not None else sum(abs(face@P[:,:6]))
        tests.append(excess)
    return dict(status='contained' if max(tests)<-.001 else 'unresolved',remainder_km=R,
                max_outside_plane_km=max(tests),ref_rp=rp,seconds=perf_counter()-start)
if __name__=='__main__':
    filename='survival-cells.json'
    # Complete two refinement passes over the whole incoming set before interpreting the pilot.
    # Record every unresolved leaf. Classification never treats a numeric LP failure as empty.
    old=json.loads((D/'slabs.json').read_text())
    assert len(old)==32
    leaves=[dict(id=str(r['index']),box=r['box'],result=r['result'],depth=0) for r in old]
    for depth in range(2):
        nextlevel=[]
        for parent in leaves:
            if parent['result']['status']!='unresolved':nextlevel.append(parent);continue
            box=np.array(parent['box']);wid=np.ptp(box,axis=1)
            # Choose a split by the upstream encounter displacement. Include timing along track.
            priority=wid/np.array([200,65,40,40,5,5])
            k=np.argmax(priority);middle=box[k].mean()
            for nchild in range(2):
                child=box.copy();child[k,1-nchild]=middle;child=tightbox(child)
                test=evaluate(child)
                nextlevel.append(dict(id=parent['id']+'.'+str(nchild),depth=depth+1,box=child.tolist(),result=test))
        leaves=nextlevel
        save(filename,dict(depth=depth+1,counts=dict(Counter(r['result']['status'] for r in leaves)),leaves=leaves))
        print(depth+1,Counter(r['result']['status'] for r in leaves),flush=True)
