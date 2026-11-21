"""Recompute the limiting audited cell and export the four signed-face budgets."""
import sys,json
from pathlib import Path
from fractions import Fraction
D=Path(__file__).resolve().parent.parent/'square-decision'
sys.path.insert(0,str(D))
import joint_cells as jc
import numpy as np
from warm_support import CachedDual
key='20.1.0.1.1.0.0.1.0.1.1.1.1.1.0'
def read(n):return json.loads((D/(n+'.json')).read_text())
jc.install()
jc.sc.support=CachedDual()
box=np.array(read('audited-boxes')[key]);hist=[]
rows=read('partition-routes')
chain=key.split('.')
for i in range(1,len(chain)):
    parent='.'.join(chain[:i]);hist.append(dict(parent=parent,child=int(chain[i]),**rows[parent]))
budgets=[]
linear={}
goals=[];multipliers=[];scales=[]
def record(A,b,l,h,g,ub,y,f):
    if ub is None or 'face' not in f:return
    y=np.maximum(0.,y)
    residual=list(map(Fraction,g));s=Fraction(0)
    for j in np.flatnonzero(y):
        w=Fraction(y[j]);s+=w*Fraction(b[j]);row=A.getrow(j)
        for col,a in zip(row.indices,row.data):residual[col]-=w*Fraction(a)
    ex=s+sum(max(r*Fraction(lo),r*Fraction(hi)) for r,lo,hi in zip(residual,l,h))
    scale=f['scale_goal']
    assert float(ex)<=ub
    budgets.append(dict(face=f['face'].tolist(),reference_km=float(f['face']@f['ray']),
        remainder_km=f['R'],site_ephemeris_km=.021,light_time_km=50*f['delay'],
        dual_upper_km=ub*scale,exact_dual_km=float(ex)*scale,rounding_km=.0001,
        row_count=A.shape[0],variables=A.shape[1],n_knots=len(f['t'])))
    if not linear:linear.update(indptr=A.indptr,indices=A.indices,values=A.data,b=b,low=l,high=h)
    goals.append(g);multipliers.append(y);scales.append(scale)
class Dual:
    def __init__(self):self.lp=CachedDual()
    def __call__(self,A,b,l,h,g):
        ub=self.lp(A,b,l,h,g)
        if ub is not None:record(A,b,l,h,g,ub,-np.array(self.lp.h.getSolution().row_dual),sys._getframe(1).f_locals)
        return ub
jc.cert=Dual()
orig_support=jc.support
def fallback(A,b,l,h,g):
    ub=orig_support(A,b,l,h,g)
    if ub is not None:
        opt=jc.linprog(-g,A_ub=A,b_ub=b,bounds=np.stack([l,h],axis=1),method='highs')
        record(A,b,l,h,g,ub,-opt.ineqlin.marginals,sys._getframe(1).f_locals)
    return ub
jc.support=fallback
q=jc.eval_joint(box,local_optical=True)
assert len(budgets)==4,(q,len(budgets))
for b in budgets:
    b['total_km']=sum(b[k] for k in ['reference_km','remainder_km','site_ephemeris_km','light_time_km','dual_upper_km','rounding_km'])
assert abs(max(b['total_km'] for b in budgets)-q['max_outside_plane_km'])<1e-8
if len(sys.argv)>1:
    np.savez_compressed(sys.argv[1],**linear,goals=goals,multipliers=multipliers,scales=scales)
print(json.dumps(dict(id=key,route=hist,box=box.tolist(),last_state=(jc.origin+jc.inverse@box.mean(axis=1)).tolist(),faces=budgets,result=q),indent=2))
