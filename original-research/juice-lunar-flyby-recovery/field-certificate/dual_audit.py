"""Replay final LP supports using exact binary-rational arithmetic on the saved inequalities."""
from fractions import Fraction as F
from global_cover import *
tag="nominal"
if __name__=="__main__":
    if not (HERE/(tag+"-LP.npz")).is_file():
        fit_cover(Arc("juice"),6.,2.,0.,0.,0.,audit_tag=tag)
    data=np.load(HERE/(tag+"-LP.npz"));A,b,bounds=[data[k] for k in ("A","b","bounds")]
    cs=json.loads((HERE/(tag+"-duals.json")).read_text())
    assert len(cs)==12
    dif=[]
    for item in cs:
        residual=list(map(F,item["c"]));v=F(0)
        for i,q in zip(item["active"],item["multipliers"]):
            y=F(q);assert y>=0
            v+=y*F(b[i])
            for j in np.flatnonzero(A[i]):residual[j]-=y*F(A[i,j])
        v+=sum(max(f*F(lo),f*F(hi)) for f,(lo,hi) in zip(residual,bounds))
        diff=F(item["rounded_upper"])-v
        assert diff>=0,float(diff)
        dif.append(float(diff))
    save("dual-check.json",dict(checked=len(cs),outward_margins= dif,rows=A.shape[0],variables=A.shape[1]))
