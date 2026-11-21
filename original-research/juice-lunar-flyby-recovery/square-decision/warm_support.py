"""Reuse a HiGHS LP basis across objective directions. Evaluate returned duals with an outward budget."""
import numpy as np
from scipy.optimize._highspy import _core as hs
class CachedDual:
    def __init__(self):self.keys=None
    def __call__(self,m,b,l,u,goal):
        if self.keys is None or any(a is not z for a,z in zip(self.keys,(m,b,l,u))):
            self.keys=(m,b,l,u)
            lp=hs.HighsLp();lp.num_col_=m.shape[1];lp.num_row_=m.shape[0]
            lp.col_cost_=np.zeros(m.shape[1]);lp.col_lower_=l;lp.col_upper_=u
            lp.row_lower_=np.full(len(b),-np.inf);lp.row_upper_=b
            v=lp.a_matrix_;v.format_=hs.MatrixFormat.kRowwise
            v.start_=m.indptr;v.index_=m.indices;v.value_=m.data
            h=hs._Highs();h.setOptionValue('output_flag',False);h.setOptionValue('threads',1)
            h.passModel(lp);self.h=h;self.idx=np.arange(len(l),dtype=np.int32)
        self.h.changeColsCost(len(goal),self.idx,-goal)
        self.h.run()
        if self.h.getModelStatus()!=hs.HighsModelStatus.kOptimal:return None
        y=np.maximum(0.,-np.array(self.h.getSolution().row_dual));idx=np.flatnonzero(y)
        w=y[idx].astype(np.longdouble);sub=m[idx].toarray().astype(np.longdouble);hb=b[idx].astype(np.longdouble)
        r=goal.astype(np.longdouble)-w@sub
        pad=(np.sum(abs(w*hb))+np.sum((w@abs(sub)+abs(goal))*np.maximum(abs(l),abs(u)))+1)*(len(y)+len(goal))*1e-15
        return float(w@hb+sum(np.maximum(r*l,r*u))+pad)+1e-7
if __name__=='__main__':
    from adaptive import *
    box=np.array(json.loads((D/'adaptive.json').read_text())['leaves'][-1]['box']);m,b=cell_program(box);w=CachedDual()
    for k in range(6):
        t=perf_counter();v=w(m,b,lb,ub,GG[k]);s1=perf_counter()-t
        t=perf_counter();v0=support(m,b,lb,ub,GG[k]);s2=perf_counter()-t
        print(k,s1,s2,v,v0)
