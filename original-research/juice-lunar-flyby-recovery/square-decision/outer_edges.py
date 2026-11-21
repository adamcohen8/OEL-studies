"""Maximize detector-edge displacement inside the tightened incoming relaxation."""
from refined_test import *
from scipy.sparse import vstack
from functools import lru_cache
mat,bv=setup_full()
endidx=slice(6*n-6,6*n);origin=ref[endidx]
def lastp(p):return path(state(p))[endidx]
J=np.column_stack([(lastp(h)-lastp(-h))/.02 for h in .01*np.eye(6)])
inv=np.linalg.inv(J)
sel=np.zeros((6,N));sel[:,endidx]=inv*scale[endidx];bounds0=(bounds-ref[:,None])/scale[:,None]
def last(z):return origin+J@z
@lru_cache(maxsize=4000)
def val(args):
    z=np.asarray(args);l=last(z);return np.r_[sky_last(l),closest(l)]
def run(axis,seed):
    z=inv@(np.asarray(seed)-origin);radius=3
    h=.015
    def obj(k):return k[:3]@axis*57.3
    for it in range(24):
        v=val(tuple(z))
        jac=np.column_stack([(val(tuple(z+e))-val(tuple(z-e)))/(2*h) for e in h*np.eye(6)])
        grad=axis@jac[:3];rgrad=jac[3]
        rows=vstack([mat,sel,-sel,csr_matrix((-rgrad@sel)[None,:])],format='csr')
        rhs=np.r_[bv,z+radius,radius-z,v[3]-1742-rgrad@z]
        c=grad@sel
        res=linprog(-c/max(abs(c)),A_ub=rows,b_ub=rhs,bounds=bounds0,method='highs')
        if not res.success:break
        new=sel@res.x
        for f in (1.,.75,.5,.25,.125,.05):
            p=z+f*(new-z);w=val(tuple(p))
            if w[3]>=1741.98 and obj(w)>obj(v)+1e-8:break
        else:
            radius/=3
            if radius<.01:break
            continue
        z=p;radius=min(5,radius*1.25)
    return dict(last=last(z).tolist(),los=sky_last(last(z)).tolist(),minimum_moon_km=closest(last(z)))
if __name__=='__main__':
    pts=json.loads((D/'refined-test.json').read_text())['points'];basepts=json.loads((D/'sky-support.json').read_text())['points']
    pts+=[dict(last=path(np.array(p['state']))[endidx].tolist(),los=p['los']) for p in basepts]
    mid=o.unit(np.mean([p['los'] for p in pts],axis=0));e,v=o.basis(mid);new=[]
    for deg in [79,169,259,349,35,215]:
        ang=np.deg2rad(deg);axis=np.cos(ang)*e+np.sin(ang)*v
        seed=max(pts,key=lambda r:np.array(r['los'])@axis)
        q=run(axis,seed['last']);pts.append(q);new.append(q)
        shape=fit([p['los'] for p in pts]);save('outer-edges.json',dict(points=new,fit=shape))
        print(deg,shape['fit_side_deg'],shape['lower_test_side_deg'],flush=True)
