"""Adaptive all-orbit detector test with explicit bookkeeping of unresolved initial sets."""
from survival_cells import *
import survival_cells as sc
from fast_reference import fasttrack
from scipy.sparse import csr_matrix,hstack
def install():
    con=np.load(D/'contracted.npz');Q=con['A'];bs=con['b'];V=Q*scale
    s=np.max(abs(V),axis=1)
    mm=csr_matrix(V/s[:,None]);bb=(bs-Q@ref)/s
    sc.fullM=vstack([fullM,mm],format='csr');sc.fullb=np.r_[fullb,bb]+2e-9
    sc.lb=np.maximum(lb,(con['bounds'][:,0]-ref)/scale-1e-8)
    sc.ub=np.minimum(ub,(con['bounds'][:,1]-ref)/scale+1e-8)
    sc.track=fasttrack
    sc.allow_deep_reference=True
def impossible(box):
    # Feasibility of the cell in the incoming optical polytope. On a solver failure,
    # use a bounded phase-I slack with an outward dual bound before discarding it.
    m,b=cell_program(box)
    r=linprog(np.zeros(N),A_ub=m,b_ub=b,bounds=np.stack([sc.lb,sc.ub],axis=1),method='highs')
    if r.success:return False
    col=np.r_[np.zeros(len(sc.fullb)),-np.ones(12)]
    mat=hstack([m,csr_matrix(col[:,None])],format='csr')
    goal=np.r_[np.zeros(N),-1.]
    q=support(mat,b,np.r_[sc.lb,0.],np.r_[sc.ub,1e5],goal)
    return q is not None and q< -1e-6
def run_one(b):
    if impossible(b):return b,dict(status='optically empty')
    box=tightbox(b)
    if min(np.ptp(box,axis=1))<=0:return b,dict(status='unresolved',reason='degenerate numerical box')
    return box,evaluate(box)
def start():
    rows=json.loads((D/'slabs.json').read_text())
    # Reuse only containing proofs with the same camera and unchanged full-domain assumptions.
    return [dict(id=str(r['index']),depth=0,box=r['box'],result=r['result']) for r in rows]
def run(levels):
    install();pathout=D/'adaptive.json'
    if pathout.is_file():leaves=json.loads(pathout.read_text())['leaves']
    else:leaves=start()
    for _ in range(levels):
        todo=[p for p in leaves if p['result']['status']=='unresolved']
        if not todo:break
        depth=min(p['depth'] for p in todo)
        todo=[p for p in todo if p['depth']==depth]
        mapping={p['id']:p for p in leaves}
        starttime=perf_counter()
        def persist():
            leaves=list(mapping.values());save('adaptive.json',dict(counts=dict(Counter(p['result']['status'] for p in leaves)),leaves=leaves))
        for j,parent in enumerate(todo):
            b=np.array(parent['box']);wid=np.ptp(b,axis=1)
            k=np.argmax(wid*np.array([1.4,1.,.9,1.,.5,.5])) if depth>=5 else np.argmax(wid/np.array([170,55,32,26,2.5,2.5]));middle=b[k].mean()
            pair=[]
            for side in range(2):
                child=b.copy();child[k,1-side]=middle
                child,res=run_one(child)
                pair.append(dict(id=parent['id']+'.'+str(side),depth=depth+1,box=child.tolist(),result=res))
            del mapping[parent['id']]
            mapping.update({p['id']:p for p in pair})
            if j%12==0:persist()
        persist();leaves=list(mapping.values())
        print('depth',depth+1,Counter(p['result']['status'] for p in leaves),'seconds',perf_counter()-starttime,flush=True)
if __name__=='__main__':run(int(sys.argv[1]))
