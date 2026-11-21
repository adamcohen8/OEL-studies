"""Try detector inclusion on correlated optical subsets and retain every unresolved slab."""
from affine_cells import *
from refined_test import setup_full
from outer_test import N,scale,ref,bounds,n
from scipy.sparse import vstack,csr_matrix
# For the pilot, use the square optimized over relaxed terminal trajectories.
r=json.loads((D/'outer-edges.json').read_text())['fit'];c=np.array(r['center']);e,v=o.basis(c);t=np.deg2rad(r['roll_deg'])
sides=np.array([e*np.cos(t)+v*np.sin(t),-e*np.sin(t)+v*np.cos(t)])
import affine_cells as ac
ff=np.vstack([sides,-sides])-np.tan(np.deg2rad(2.2))*c;ac.faces=ff/np.linalg.norm(ff,axis=1)[:,None]
def support(m,b,l,u,goal):
    res=linprog(-goal,A_ub=m,b_ub=b,bounds=np.stack([l,u],axis=1),method='highs')
    if not res.success:return None
    y=np.maximum(0.,-res.ineqlin.marginals);idx=np.flatnonzero(y);w=y[idx].astype(np.longdouble)
    sub=m[idx].toarray().astype(np.longdouble);hb=b[idx].astype(np.longdouble)
    r=goal.astype(np.longdouble)-w@sub
    pad=(np.sum(abs(w*hb))+np.sum((w@abs(sub)+abs(goal))*np.maximum(abs(l),abs(u)))+1)*(len(y)+N)*1e-15
    return float(w@hb+sum(np.maximum(r*l,r*u))+pad)+1e-7
if __name__=='__main__':
    m,b=setup_full();base=np.array(cfg['box']);origin=np.array(cfg['origin'])
    goals=np.zeros((6,N));goals[:,6*n-6:6*n]=transform*scale[6*n-6:6*n]
    # Scale the chart constraints before giving them to the LP solver.
    scales=np.abs(goals).max(axis=1);gg=goals/scales[:,None]
    lower,upper=((bounds-ref[:,None])/scale[:,None]).T
    edges=np.linspace(*base[0],33);out=[]
    for i,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):
        m1=vstack([m,gg[0],-gg[0]],format='csr');b1=np.r_[b,hi/scales[0],-lo/scales[0]]
        local=np.array(base);local[0]=[lo,hi]
        for k in range(1,6):
            low=support(m1,b1,lower,upper,-gg[k]);high=support(m1,b1,lower,upper,gg[k])
            if low is not None and high is not None:
                local[k]=[max(base[k,0],-low*scales[k]),min(base[k,1],high*scales[k])]
        center=local.mean(axis=1);width=np.ptp(local,axis=1)/2
        last=origin+inverse@center
        # Approximate closest approach is a scheduling aid. Skipped or failed cells stay unresolved.
        from refined_test import closest
        rp=closest(last)
        row=dict(index=i,box=local.tolist(),nominal_peri_km=rp)
        if rp>1000:
            start=perf_counter();bundle=track(last);q=test(bundle,width)
            row.update(result=q,track_and_test_seconds=perf_counter()-start)
        else:row['result']=dict(status='unresolved',reason='near-center reference')
        out.append(row);save('slabs.json',out)
        if i%4==0: print(i,row['result'],flush=True)
