"""Propagate correlated terminal states from the 32-facet optical relaxation in OEL."""
from outer_test import *
from footprint import fit
from scipy.spatial.distance import pdist,squareform
fc_time=a.times[n-1];end=a.times[-1]
def propagate(last,t):return a.propagate(last,t-fc_time,fc_time)
def reverse(last):
    # The OEL API integrates positive durations. Match the endpoint by shooting forward
    # using the weakly nonlinear state map on the historical observation arc.
    p=np.zeros(6)
    def last_of(p):return path(state(p))[6*n-6:6*n]
    for _ in range(5):
        err=last-last_of(p)
        if np.linalg.norm(err[:3])<1e-6 and np.linalg.norm(err[3:])<1e-11:break
        jac=np.column_stack([(last_of(p+h)-last_of(p-h))/.02 for h in .01*np.eye(6)])
        p+=np.linalg.solve(jac,err)
    return p
def feasible_supports():
    mat,bv=Q1;lohi=(bounds-ref[:,None])/scale[:,None]
    ind=slice(6*n-6,6*n)
    def terminal(p):return path(state(p))[ind]
    origin=terminal(np.zeros(6));jac=np.column_stack([(terminal(h)-terminal(-h))/.02 for h in .01*np.eye(6)])
    inv=np.linalg.inv(jac)
    axes=list(np.eye(6))+list(-np.eye(6))
    # Optical weak directions, along with directions of the downstream sky gradient.
    rng=np.random.default_rng(91301);axes+=list(rng.normal(size=(180,6)))
    for q in json.loads((D/'sky-support.json').read_text())['points'][:8]:
        p=coords(np.array(q['state']))
        skyjac=np.column_stack([(a.sky(state(p+h),n+1)[-1]-a.sky(state(p-h),n+1)[-1])/.04 for h in .02*np.eye(6)])
        center=np.array(q['los']);e,v=o.basis(center)
        for ang in np.linspace(0,2*np.pi,12,endpoint=False):
            axes.append((e*np.cos(ang)+v*np.sin(ang))@skyjac)
    chosen=[]
    for u in axes:
        c=np.zeros(N);c[ind]=(u@inv)*scale[ind]
        sol=linprog(-c/max(abs(c)),A_ub=mat,b_ub=bv,bounds=lohi,method='highs')
        if sol.success:
            yy=ref+scale*sol.x;z=yy[ind]
            if not any(np.linalg.norm((q['last']-z)@inv.T)<.01 for q in chosen):
                chosen.append(dict(last=z,lp=yy.tolist(),max_scaled_violation=float(max(0.,max(mat@sol.x-bv)))))
    return chosen
if __name__=='__main__':
    rows=feasible_supports();pool=[]
    # Join each relaxed extremum to a genuine, admitted grazing or distant solution.
    near=min(old,key=lambda q:q['encounter']['moon_radius_km']);z0=path(np.array(near['state']))[6*n-6:6*n]
    count=0
    for row in rows:
        z=np.array(row['last']);p=reverse(z);w=data(p)
        if w:
            w.update(last=z.tolist(),lp=row['lp'],max_scaled_violation=row['max_scaled_violation']);pool.append(w)
        else:
            # A surface intersection can hide a feasible angular extremum immediately
            # outside the sphere. Bisect along a convex combination of admitted terminals.
            lo=0.;hi=1.
            for _ in range(14):
                f=(lo+hi)/2;pp=reverse((1-f)*z0+f*z)
                if a.encounter(state(pp))['moon_radius_km']>1742:lo=f
                else:hi=f
            pp=reverse((1-lo)*z0+lo*z);w=data(pp)
            if w:
                w.update(last=((1-lo)*z0+lo*z).tolist(),lp_fraction=lo,extreme_last=z.tolist())
                pool.append(w)
        count+=1
        if count%15==0 or count==len(rows):
            shape=fit([q['los'] for q in pool]);print(count,len(rows),shape['fit_side_deg'],shape['lower_test_side_deg'],flush=True)
            save('outer-sky.json',dict(rows=count,total=len(rows),fit=shape,points=pool))
            if shape['excludes_4p4'] and shape['lower_test_side_deg']>4.43:break
