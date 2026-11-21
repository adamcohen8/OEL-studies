"""Tighten the incoming program by integrating differences from an enclosed OEL trajectory."""
from outer_test import *
from tubes import E,M,S,prepare,enclosure
from functools import lru_cache
MUS=np.array([E,M,S])
Vbody=np.array([0.,2.,35.])
def additional_rows():
    # Use emission-independent reception knots and subinterval midpoints in the reference run.
    ts=a.times[:n];pieces=np.maximum(1,np.ceil(np.diff(ts)/240).astype(int))
    midpoints=[ts[j]+(np.arange(k)+.5)*(ts[j+1]-ts[j])/k for j,k in enumerate(pieces)]
    # The validation only needs the training span. Supply a short view with the terminal site.
    import copy
    segment=copy.copy(a);segment.times=ts;segment.sites=a.sites[:n]
    refdata=prepare(segment,base,.25,np.unique(np.r_[ts[1:],np.concatenate(midpoints)]));records={}
    assert enclosure(refdata,segment,0.,0.,record=records)['status']=='bounded'
    records[0.]=(0.,0.);nodes={t:x for t,x,_,_ in refdata}
    means=[nodes[t] for t in ts]
    knotR=np.array([records[t][0] for t in ts]);knotV=np.array([records[t][1] for t in ts])
    box=bounds[:6*n].reshape(n,6,2);r0=np.max(abs(box[:,:3,:]-np.array(means)[:,:3,None]),axis=2)
    delta=np.linalg.norm(r0,axis=1)+knotR
    rows=[];rhs=[];widths=[]
    for j,mids in enumerate(midpoints):
        h=ts[j+1]-ts[j];ds=h/len(mids);dmax=max(delta[j:j+2])
        # Both trajectories satisfy the declared training floors. A uniform acceleration
        # cap applies to their difference before refinement by the actual ephemeris.
        wide=dmax+.000024*h*h/8+1
        pos=np.array([nodes[t][:3] for t in mids])
        bod=[o.hpop_de440_positions_km(a.JD+t/86400,a.ephem) for t in mids]
        body=np.array([[np.zeros(3),z['moon'],z['sun']] for z in bod])
        rel=pos[:,None,:]-body;length=np.linalg.norm(rel,axis=2)
        refR=max(records[t][0] for t in mids)
        travel=(6.+Vbody)*ds/2+refR+.02
        def constants(wide):
            floor=np.min(length-travel-wide,axis=0)
            assert min(floor)>2e5
            lips=np.sum(2*MUS/floor**3)
            j2=2*12*.001083*E*6379**2/floor[0]**4
            return floor,lips,j2
        for _ in range(3):
            floor,lips,j2=constants(wide)
            assert 1-h*h*lips/8>0
            wide=min(wide,(dmax+h*h*j2/8+.00001)/(1-h*h*lips/8))
        floor,lips,j2=constants(wide)
        chord=(lips*wide+j2)*h*h/8+1e-7
        curv=3.00001*MUS*wide**2/floor**4
        # Tensor of the derivative of the gravity Jacobian along a direction is
        # bounded by 24 mu/r^4. This term encloses its drift away from each midpoint.
        drift=24*MUS*travel*wide/floor**4
        err=sum(curv+drift)+j2+lips*chord+1e-15
        J=MUS[None,:,None,None]*(3*rel[:,:,:,None]*rel[:,:,None,:]/length[:,:,None,None]**5-np.eye(3)/length[:,:,None,None]**3)
        J=J.sum(axis=1);wl=np.zeros((3,3,3));wr=wl.copy();absolute=np.zeros(3)
        for k,mat in enumerate(J):
            for g in (-1.,1.):
                tau=(k+.5+g/np.sqrt(12))/len(mids)
                weights=ds/2*np.array([1.,h*(1-tau),h*(.5-tau)])
                wl+=weights[:,None,None]*(1-tau)*mat;wr+=weights[:,None,None]*tau*mat
        absolute=np.array([h,h*h/2,h*h/4])
        e=absolute*err
        for axis in range(3):
            for typ in range(3):
                q=np.zeros(N);i=6*j+axis
                if typ==0:q[i+9],q[i+3]=1,-1
                elif typ==1:q[i+6],q[i],q[i+3]=1,-1,-h
                else:q[i+6],q[i],q[i+3],q[i+9]=1,-1,-h/2,-h/2
                q[6*j:6*j+3]-=wl[typ,axis];q[6*j+6:6*j+9]-=wr[typ,axis]
                qlocal=q[6*j:6*j+12]
                center=qlocal@np.r_[means[j],means[j+1]]
                referr=np.abs(qlocal)@np.repeat([knotR[j],knotV[j],knotR[j+1],knotV[j+1]],3)
                allowance=e[typ]+referr+(2e-9 if typ==0 else .0002)
                rows.extend([q,-q]);rhs.extend([center+allowance,-center+allowance])
        widths.append(dict(gap_s=h,displacement=wide,chord=chord,position_quadrature_allowance_km=e[1]))
    return np.array(rows),np.array(rhs),widths
if __name__=='__main__':
    Q,bq,w=additional_rows()
    np.savez_compressed(D/'differential-constraints.npz',A=Q,b=bq)
    save('differential-budget.json',w)
    print('long gaps',sorted(w,key=lambda u:u['gap_s'],reverse=True)[:3],flush=True)
