"""Couple quadratic sky displacements and lunar survival in one lifted optical program."""
from adaptive import *
from fine_reference import fine_track,fine_tube,fine_grid
from scipy.sparse import hstack
from itertools import combinations_with_replacement
from warm_support import CachedDual
from math import sqrt
pairs=list(combinations_with_replacement(range(6),2));cert=CachedDual()
def eval_joint(box,local_optical=False):
    mid=box.mean(axis=1);width=np.ptp(box,axis=1)/2;last=origin+inverse@mid
    t,x,acc,d,j,ray,light=fine_track(last)
    if not np.all(np.isfinite(x)):return dict(status='unresolved',reason='nonfinite reference')
    per=t[d[:,1].argmin()]
    ids=np.unique(np.clip(np.searchsorted(t,per+np.arange(-720,721,30)),0,len(t)-2))
    sampled=np.zeros((len(ids),3,27));sr=np.full(len(ids),-1.)
    R,P=fine_tube(t,x,acc,d,j,inverse@np.diag(width),ids,sampled,sr)
    if R<0:return dict(status='collision',method='quadratic position')
    m,b=cell_program(box)
    nvar=N+27;sp=hstack([m,csr_matrix((len(b),27))],format='csr')
    low=np.r_[sc.lb,-np.ones(6),[0 if i==j else -1 for i,j in pairs]]
    high=np.r_[sc.ub,np.ones(27)]
    rows=[];rhs=[]
    def append(line,value):
        size=max(1.,max(abs(line)));rows.append(line/size);rhs.append((value+1e-7)/size)
    # p_i=(G_i y-mid_i)/width_i, enclosed by the normalized coordinate box.
    grad=G/width[:,None];bias=mid/width
    for i in range(6):
        row=np.r_[-grad[i],np.eye(6)[i],np.zeros(21)]
        append(row,-bias[i]);append(-row,bias[i])
    def form(p,q,bound):
        append(np.r_[np.zeros(N),p,q],bound)
    for k,(i,jj) in enumerate(pairs):
        q=np.eye(21)[k];ei=np.eye(6)[i];ej=np.eye(6)[jj]
        if i==jj:
            for h in [-1,-.5,0,.5,1]:form(2*h*ei,-q,h*h)
        else:
            form(ei+ej,-q,1);form(-ei-ej,-q,1)
            form(ei-ej,q,1);form(-ei+ej,q,1)
    phase_start=len(rows)
    if local_optical:
        from local_optics import constraints
        oa,ob,_=constraints(last,width)
        for optical_row,limit in zip(oa,ob):form(optical_row,np.zeros(21),limit)
    _,moons,_,_=fine_grid()
    for k,PP,rr in zip(ids,sampled,sr):
        if rr<0:continue
        pos=x[k,:3]-moons[k];length=np.linalg.norm(pos);u=pos/length
        L,Q=PP[:,:6],PP[:,6:];lin=maxball(L);quad=maxball(Q)
        error=rr+.021
        cut=(1737.4**2-length**2)/(2*length)-error-(2*lin*(quad+error)+(quad+error)**2)/(2*length)-.00002
        cL=u@L;cQ=u@Q+np.array([(L[:,i]@L[:,jj])*(1 if i==jj else 2)/(2*length) for i,jj in pairs])
        if cut>-sum(abs(cL))-sum(abs(cQ)):form(-cL,-cQ,-cut)
    sp=vstack([sp,csr_matrix(np.array(rows))],format='csr');bv=np.r_[b,rhs]
    # A proof that the lifted problem is empty also certifies that no surviving trajectory remains.
    res=linprog(np.zeros(nvar),A_ub=sp,b_ub=bv,bounds=np.stack([low,high],axis=1),method='highs')
    if not res.success:
        phase=hstack([sp,csr_matrix(np.r_[np.zeros(len(b)+phase_start),-np.ones(len(rhs)-phase_start)][:,None])],format='csr')
        bound=cert(phase,bv,np.r_[low,0],np.r_[high,1e6],np.r_[np.zeros(nvar),-1.])
        if bound is not None and bound< -1e-8:return dict(status='infeasible',method='empty lifted relaxation',phase_upper=bound)
        from fractions import Fraction as F
        target=np.r_[np.zeros(nvar),1.]
        phase_sol=linprog(target,A_ub=phase,b_ub=bv,bounds=np.stack([np.r_[low,0.],np.r_[high,1e6]],axis=1),method='highs')
        if phase_sol.success:
            y=np.maximum(0.,-phase_sol.ineqlin.marginals)
            resid=[F(0) for _ in low];upper=F(0)
            for ii in np.flatnonzero(y):
                lam=F(y[ii]);upper+=lam*F(bv[ii]);r=sp.getrow(ii)
                for kk,aa in zip(r.indices,r.data):resid[kk]-=lam*F(aa)
            upper+=sum(max(v*F(l),v*F(h)) for v,l,h in zip(resid,low,high))
            if upper<0:return dict(status='infeasible',method='rational infeasibility check',upper=float(upper))
        return dict(status='unresolved',reason='feasibility was not certified')
    if R>1e8:return dict(status='unresolved',reason='quadratic remainder')
    delay=(abs(light-np.linalg.norm(ray))+maxball(P)+R+.021)/(o.LIGHT_SPEED-50)
    upper=[]
    for face in ac.faces:
        lin=face@P[:,:6];q=face@P[:,6:];goal=np.r_[np.zeros(N),lin,q]
        scale_goal=max(1,max(abs(goal)))
        u=cert(sp,bv,low,high,goal/scale_goal)
        if u is None:
            u=support(sp,bv,low,high,goal/scale_goal)
        if u is not None:u*=scale_goal
        if u is None:return dict(status='unresolved',reason='objective')
        upper.append(face@ray+R+.021+50*delay+u+.0001)
    return dict(status='contained' if max(upper)<-.001 else 'unresolved',method='lifted projection',
                max_outside_plane_km=max(upper),remainder_km=R,reference_peri_km=min(d[:,1]))
if __name__=='__main__':
    install();q=json.loads((D/'quadratic-cover.json').read_text())['leaves']
    sample=[r for r in q if r['result']['status']=='unresolved']
    for r in sample[:6]:
        print(r['id'],eval_joint(np.array(r['box'])),flush=True)
