"""Recompute optical and gravity difference bounds after contracting every training knot."""
from survival_cells import *
import differential as df
import outer_test as ot
if __name__=='__main__':
    new=bounds.copy()
    for j in range(N):
        e=np.zeros(N);e[j]=1
        high=support(fullM,fullb,lb,ub,e);low=support(fullM,fullb,lb,ub,-e)
        if low is None or high is None:continue
        new[j]=[max(bounds[j,0],ref[j]-low*scale[j]-1e-10),min(bounds[j,1],ref[j]+high*scale[j]+1e-10)]
    df.bounds=new;ot.bounds=new
    opt,bopt=ot.ray_rows(96);dyn,bdyn,diag=df.additional_rows()
    np.savez_compressed(D/'contracted.npz',A=np.vstack([opt,dyn]),b=np.r_[bopt,bdyn],bounds=new)
    save('contracted-budget.json',diag)
    print([(u['gap_s'],u['position_quadrature_allowance_km']) for u in diag if u['gap_s']>20000])
