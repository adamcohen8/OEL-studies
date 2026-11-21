"""Project the tightened incoming program with outward dual bounds and replay the complete flyby cover."""
from refined_test import *
from branch_cover import solve
def support(mat,bv,obj):
    v=mat*scale;ss=np.max(abs(v),axis=1)
    cs=max(abs(obj*scale));m=csr_matrix(v/ss[:,None]);bc=(bv-mat@ref)/ss
    res=linprog(-obj*scale/cs,A_ub=m,b_ub=bc,bounds=(bounds-ref[:,None])/scale[:,None],method='highs')
    if not res.success:raise RuntimeError(res.message)
    idx=np.flatnonzero(res.ineqlin.marginals<0)
    lam=-res.ineqlin.marginals[idx]*cs/ss[idx]
    # Binary-rational evaluation removes dual-feasibility and rounding ambiguity in this projection.
    from fractions import Fraction as F
    c=[F(q) for q in obj];val=F(0)
    for i,y in zip(idx,lam):
        w=F(float(y));assert w>=0
        val+=w*F(bv[i])
        for j in np.flatnonzero(mat[i]):c[j]-=w*F(mat[i,j])
    val+=sum(max(z*F(l),z*F(u)) for z,(l,u) in zip(c,bounds))
    return float(val)+1e-7
if __name__=='__main__':
    mat,bv=matrices();oldrow=json.loads((fc/'global.json').read_text())[0];row=oldrow.copy()
    basis=np.array(row['basis'])
    outs=[]
    for axes in (np.eye(6),block_diag(basis,basis)):
        cub=[]
        for axis in axes:
            q=np.zeros(N);q[6*n-6:6*n]=axis
            cub.append([-support(mat,bv,-q),support(mat,bv,q)])
        outs.append(cub)
    cub=np.array(outs[0]);row['last_state_lo']=cub[:,0].tolist();row['last_state_hi']=cub[:,1].tolist();row['last_tangent_box']=outs[1]
    full=solve(a,row,24,12000,directional=True,shift=1.)
    save('outer-replay.json',dict(former_box=oldrow['last_tangent_box'],tightened_box=outs[1],projection_rows=len(bv),cover=full))
    print('full cover',full['diameter_deg'],flush=True)
