"""Check the binary partition and independently bound any terminal state omitted by a leaf."""
from adaptive import *
from warm_support import CachedDual
from fractions import Fraction as F
from itertools import product
def exact(h,m,b,goal):
    if h.h.getModelStatus()!=__import__('scipy').optimize._highspy._core.HighsModelStatus.kOptimal:return None
    dual=-np.array(h.h.getSolution().row_dual);rr=[F(v) for v in goal];ub=F(0)
    for k in np.flatnonzero(dual>0):
        w=F(dual[k]);ub+=w*F(b[k]);row=m.getrow(k)
        for ix,aa in zip(row.indices,row.data):rr[ix]-=w*F(aa)
    ub+=sum(max(z*F(l),z*F(u)) for z,l,u in zip(rr,sc.lb,sc.ub))
    return float(ub)+1e-14
def audit(name):
    install()
    orig=json.loads((D/name).read_text());leaves={r['id']:r for r in orig['leaves']}
    ids=set(leaves);allnodes=set()
    for key in ids:
        fields=key.split('.')
        for n0 in range(1,len(fields)+1):allnodes.add('.'.join(fields[:n0]))
    assert {k for k in allnodes if '.' not in k}==set(map(str,range(32)))
    for k in allnodes:
        children={k+'.0',k+'.1'}&allnodes
        assert (k in ids and not children) or (k not in ids and len(children)==2)
    cubes={}
    for key in sorted(allnodes,key=lambda k:-len(k)):
        if key in ids:cubes[key]=np.array(leaves[key]['box'])
        else:
            l,r=cubes[key+'.0'],cubes[key+'.1']
            cubes[key]=np.column_stack([np.minimum(l[:,0],r[:,0]),np.maximum(l[:,1],r[:,1])])
    h=CachedDual();result=[];exact_count=0;audited={};routes={}
    def visit(key,domain):
        nonlocal exact_count
        tight=cubes[key];m,b=cell_program(domain);b=b+1e-10;current=domain.copy()
        for i,sign in product(range(6),(-1,1)):
            col=int(sign>0);requested=tight[i,col]*sign
            if current[i,col]*sign<=requested+1e-10:continue
            goal=sign*GG[i];upper=h(m,b,sc.lb,sc.ub,goal)
            if upper is not None and upper*gs[i]>requested-1e-9:
                q=exact(h,m,b,goal);exact_count+=1
                if q is not None:upper=min(upper,q)
            if upper is None:
                result.append(dict(id=key,axis=i,sign=sign,excess=None));continue
            gap=upper*gs[i]-requested
            if gap>1e-8:result.append(dict(id=key,axis=i,sign=sign,excess=gap))
            current[i,col]=sign*min(current[i,col]*sign,upper*gs[i]+.0000001)
        if key in ids:
            audited[key]=np.column_stack([np.minimum(current[:,0],tight[:,0])-1e-6,np.maximum(current[:,1],tight[:,1])+1e-6]).tolist()
            return
        l,r=cubes[key+'.0'],cubes[key+'.1']
        gap=r[:,0]-l[:,1]
        candidates=np.flatnonzero(gap>=-1e-8)
        if not len(candidates):raise ValueError((key,gap.tolist()))
        # Children are named in ascending order along their split coordinate.
        rel=(np.minimum(l[:,1]-current[:,0],current[:,1]-r[:,0]))/np.maximum(1e-8,np.ptp(current,axis=1))
        axis=candidates[np.argmax(rel[candidates])]
        cut=(r[axis,0]+l[axis,1])/2
        routes[key]=dict(axis=int(axis),cut=float(cut))
        for j in range(2):
            ch=current.copy()
            if j:ch[axis,0]=max(ch[axis,0],cut)
            else:ch[axis,1]=min(ch[axis,1],cut)
            if ch[axis,0]<=ch[axis,1]:visit(key+'.'+str(j),ch)
    edges=np.linspace(*limits[0],33)
    for i in range(32):
        box=limits.copy();box[0]=edges[i:i+2];visit(str(i),box)
    assert set(audited)==ids
    save('audited-boxes.json',audited);save('partition-routes.json',routes)
    out=dict(source=name,leaf_count=len(ids),prefix_complete=True,exact_support_rechecks=exact_count,
             issues=result,counts=orig['counts'])
    save('partition-audit.json',out);print(len(ids),len(result),result[:5],flush=True)
if __name__=='__main__':audit(sys.argv[1])
