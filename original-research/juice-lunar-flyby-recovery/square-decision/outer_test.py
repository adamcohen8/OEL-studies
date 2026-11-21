"""Pin complete OEL trajectories at last observation in the nominal and refined incoming relaxations."""
from explore import *
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from scipy.linalg import block_diag
from pathlib import Path
from itertools import product
fc=D.parent/'field-certificate'
dat=np.load(fc/'nominal-LP.npz')
A,b,bounds=[dat[k] for k in ('A','b','bounds')]
n=arc;N=7*n
def path(x):
    hist=[x]
    for j in range(1,n):hist.append(a.propagate(hist[-1],a.times[j]-a.times[j-1],a.times[j-1]))
    hist=np.array(hist);positions=hist[:,:3]-a.sites[:-1]
    ranges=np.linalg.norm(positions,axis=1)/(1+np.sum((hist[:,3:]+a.earth_v[:n])*a.obs,axis=1)/o.LIGHT_SPEED)
    return np.r_[hist.ravel(),ranges]
ref=path(base)
scale=np.r_[np.tile([1e4]*3+[.05]*3,n),np.full(n,1e4)]
def ray_rows(facets):
    rows=[];bs=[]
    for j in range(n):
        rr=bounds[6*n+j];vv=bounds[j*6+3:j*6+6];rmean=rr.mean();vmean=vv.mean(axis=1)
        # Tangent projection and exact centered-product remainder. Barycentric Earth
        # acceleration, finite station/ephemeris evaluation and delay share a position allowance.
        delay=max(rr)/o.LIGHT_SPEED
        pad=.05+.00004*delay**2
        for ang in np.linspace(0,2*np.pi,facets,endpoint=False):
            u=a.east[j]*np.cos(ang)+a.north[j]*np.sin(ang)
            q=np.zeros(N);q[j*6:j*6+3]=u;q[j*6+3:j*6+6]=-u*rmean/o.LIGHT_SPEED
            q[6*n+j]=-u@(a.earth_v[j]+vmean)/o.LIGHT_SPEED-2/o.ASEC
            bilinear=np.ptp(rr)/2*(abs(u)@(np.ptp(vv,axis=1)/2))/o.LIGHT_SPEED
            rows.append(q);bs.append(u@a.sites[j]-rmean*(u@vmean)/o.LIGHT_SPEED+pad+bilinear)
    return np.array(rows),np.array(bs)
def setup(facets):
    mat,bb=A,b
    if facets:
        Q,bb2=ray_rows(facets);mat=np.vstack([A,Q]);bb=np.r_[b,bb2]
    v=mat*scale;ss=np.max(abs(v),axis=1);return csr_matrix(v/ss[:,None]),(bb-mat@ref)/ss
Q0=setup(0);Q1=setup(32)
def feasible(x,sys):
    value=path(x);fixed=(value-ref)/scale
    lohi=(bounds-ref[:,None])/scale[:,None]
    lohi[6*n-6:6*n]=fixed[6*n-6:6*n,None];mat,bv=sys
    sol=linprog(np.zeros(N),A_ub=mat,b_ub=bv,bounds=lohi,method='highs')
    if sol.success:
        actual=ref+scale*sol.x
        return dict(feasible=True,max_scaled_violation=float(max(np.max(mat@sol.x-bv),np.max(lohi[:,0]-sol.x),np.max(sol.x-lohi[:,1]))),variables=actual.tolist())
    return dict(feasible=False,status=sol.message)
if __name__=='__main__':
    from footprint import fit
    old3=json.loads((fc/'expanded-witnesses.json').read_text())['3']['witnesses']
    pool=list(old3)
    more=json.loads((D/'sky-support.json').read_text())['points']
    for q in more:
        if all(np.linalg.norm(np.array(q['state'])[:3]-np.array(r['state'])[:3])>1 for r in pool):pool.append(q)
    # Sweep outward from the 2 arcsec extremal states toward feasible 3 arcsec trajectories.
    for w in old3:
        x=np.array(w['state'])
        seed=min(more,key=lambda q:abs(q['encounter']['moon_radius_km']-w['encounter']['moon_radius_km']))
        for frac in np.linspace(.08,.96,12):
            z=(1-frac)*np.array(seed['state'])+frac*x
            p=data(coords(z))
            if p:pool.append(p)
    groups={};records=[]
    for p in pool:
        x=np.array(p['state']);row=dict(state=p['state'],los=p['los'],actual_max_as=float(maxerr(coords(x))))
        for name,q in [('archived',Q0),('polygon32',Q1)]:row[name]=feasible(x,q)
        records.append(row)
    for key in ['archived','polygon32']:
        pts=[r for r in records if r[key]['feasible']]
        if pts:groups[key]=dict(count=len(pts),square=fit([q['los'] for q in pts]),max_error_as=max(q['actual_max_as'] for q in pts))
    save('outer-test.json',dict(groups=groups,records=records))
    print(groups,flush=True)
