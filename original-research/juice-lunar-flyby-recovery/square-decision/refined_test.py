"""Combine tangent-plane optics with a reference-relative dynamics enclosure."""
import outer_sky as osky
import outer_test as ot
import differential as diff
from outer_sky import *
from scipy.optimize import minimize_scalar
from functools import lru_cache
def matrices():
    Q,bb=ray_rows(64)
    DD=np.load(D/'differential-constraints.npz')
    return np.vstack([A,Q,DD['A']]),np.r_[b,bb,DD['b']]
def setup_full():
    A,b=matrices()
    v=A*scale;ss=np.max(abs(v),axis=1);return csr_matrix(v/ss[:,None]),(b-A@ref)/ss
def closest(last):
    start=datetime(2024,8,19,19,tzinfo=timezone.utc).timestamp()-a.rows[0][0]
    def rad(t):
        pos=propagate(last,t)[:3];m=o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)['moon'];return np.linalg.norm(pos-m)
    q=minimize_scalar(rad,method='bounded',bounds=(start-3*3600,start+6*3600))
    return q.fun
def sky_last(last):
    x=propagate(last,end);obs=a.sites[-1]
    ray=x[:3]-obs;v=x[3:]+a.earth_v[-1]
    for _ in range(3):ray=x[:3]-obs-v*np.linalg.norm(ray)/o.LIGHT_SPEED
    return o.unit(ray)
if __name__=='__main__':
    qsys=setup_full();osky.Q1=qsys
    states=feasible_supports()
    near=path(np.array(min(old,key=lambda q:q['encounter']['moon_radius_km'])['state']))[6*n-6:6*n]
    assert feasible(np.array(min(old,key=lambda q:q['encounter']['moon_radius_km'])['state']),qsys)['feasible']
    pool=[]
    for row in states:
        z=row['last'];r=closest(z);frac=1.
        if r<1742:
            lo,hi=0.,1.
            for _ in range(18):
                f=(lo+hi)/2
                if closest(near+f*(z-near))>1742:lo=f
                else:hi=f
            frac=lo
        p=near+frac*(z-near)
        pool.append(dict(last=p.tolist(),los=sky_last(p).tolist(),minimum_moon_km=closest(p),fraction=frac))
    g=fit([p['los'] for p in pool]);save('refined-test.json',dict(fit=g,points=pool))
    print(len(pool),g,flush=True)
