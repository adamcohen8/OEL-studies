"""Check locally reconstructed optical inequalities on the certified trajectories."""
from local_optics import *
from outer_test import path
if __name__=='__main__':
    rng=np.random.default_rng(1551);records=[]
    for k,row in enumerate(nom['accepted']):
        true=path(np.array(row['state']))[6*n-6:6*n]
        width=np.array([65.,30.,18.,10.,.3,.3]);z=rng.uniform(-.9,.9,6)
        midpoint=true-inverse@(width*z)
        A,b,diag=constraints(midpoint,width)
        excess=float(max(A@z-b));assert excess<0
        records.append(dict(witness=k,largest_constraint_excess_km=excess,**diag))
    save('local-optical-control.json',records);print(max(r['largest_constraint_excess_km'] for r in records))
