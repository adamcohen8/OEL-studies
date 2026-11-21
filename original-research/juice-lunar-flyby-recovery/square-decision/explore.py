"""Search the bounded optical set from both sides of the nominal encounter."""
from pathlib import Path
import sys
D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'field-certificate'))
from common import Arc,load,o,np,json
from scipy.optimize import minimize
from datetime import datetime,timezone
from itertools import combinations
a=Arc('juice')
base=np.array(load('juice-forecast.json')['interior_state'])
arc=len(a.obs)
eps=1e-3
raw=np.column_stack([(a.residual(base+a.delta(e)).ravel()-a.residual(base-a.delta(e)).ravel())/(2*eps) for e in eps*np.eye(6)])
_,s,vt=np.linalg.svd(raw,full_matrices=False)
def state(p):return base+a.delta((p/s)@vt)
def coords(x):
    return s*(vt@(np.r_[a.axes@(x[:3]-base[:3]),a.axes@(x[3:]-base[3:])]/a.scale))
b=a.residual(base); jac=(raw@vt.T/s).reshape(arc,2,6)
def optical_linear(p):return b+jac@p
def maxerr(p):return np.linalg.norm(a.residual(state(p)),axis=1).max()
old=json.loads((D.parent/'field-certificate/expanded-witnesses.json').read_text())['2']['witnesses']
def data(p):
    x=state(p);r=a.encounter(x)
    if not r.get('resolved'):return None
    if r['moon_radius_km']<1738:return None
    sky=a.sky(x,arc+1)[-1]
    return dict(p=p.tolist(),state=x.tolist(),los=sky.tolist(),max_as=float(maxerr(p)),encounter=r)
def save(name,obj):(D/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':
    # Two weakest singular directions dominate incoming range and its rate. Probe their
    # whole circular support and all other coordinate extrema as independent starting states.
    seed=np.random.default_rng(56121)
    rays=list(np.eye(6))+list(-np.eye(6))
    for ang in np.linspace(0,2*np.pi,96,endpoint=False):
        v=np.zeros(6);v[-2:]=[np.sin(ang),np.cos(ang)];rays.append(v)
    rays.extend(seed.normal(size=(96,6)))
    pts=[q for q in old]
    info=[]
    for k,axis in enumerate(rays):
        opt=minimize(lambda p:-p@axis,np.zeros(6),jac=lambda p:-axis,method='SLSQP',
             constraints=dict(type='ineq',fun=lambda p:1.975**2-np.sum(optical_linear(p)**2,axis=1),
                jac=lambda p:-2*np.einsum('ij,ijk->ik',optical_linear(p),jac)),options=dict(ftol=1e-9,maxiter=80))
        p=opt.x
        trial=data(p)
        enc=trial['encounter'] if trial else a.encounter(state(p))
        info.append(dict(p=p.tolist(),min_r=enc.get('moon_radius_km'),max_as=maxerr(p)))
        if trial and trial['max_as']<1.988:pts.append(trial)
        if k%25==0:print(k,enc.get('moon_radius_km'),info[-1]['max_as'],flush=True)
    save('exploration.json',dict(s=s.tolist(),points=pts,supports=info))
    pair=max(combinations(pts,2),key=lambda z:o.distance(np.array(z[0]['los']),np.array(z[1]['los'])))
    print('maximum',o.distance(*[np.array(q['los']) for q in pair])/3600,'accepted',len(pts),'min radius',min(q['min_r'] for q in info if q['min_r']),flush=True)
