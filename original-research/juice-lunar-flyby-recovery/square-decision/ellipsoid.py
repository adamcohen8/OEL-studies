"""Propagate an ellipsoidal bound on accumulated nonlinear error using the signed gravity gradient."""
from affine_cells import *
@njit
def ellipsoid(t,x,acc,d,J,B,sample_ids,sample_P,sample_R):
    P=B.copy();cov=np.diag(np.array([3e-14]*3+[3e-22]*3))
    idx=0;I=np.eye(3);mu=np.array([E,M,S])
    for k in range(len(t)-1):
        R=np.sqrt(np.trace(cov[:3,:3]));V=np.sqrt(np.trace(cov[3:,3:]))
        if idx<len(sample_ids) and k==sample_ids[idx]:
            sample_P[idx]=P[:3];sample_R[idx]=R;idx+=1
        h=t[k+1]-t[k];width=maxball(P[:3]);sz=width+R+.00001
        if d[k,1]+sz+.02<1737.4:return -1.,P[:3].copy()
        y=d[k]-np.array([sz+6*h,sz+8*h,sz+41*h,2*h,35*h])-.02
        earth=max(6378.,y[0]);moon=max(1737.4,y[1]);sun=max(1.e8,y[2]);bm=y[3];bs=y[4]
        if bm<3e5 or bs<1e8:return 1e9,P[:3].copy()
        j2=12*.001083*E*6379**2/earth**4
        A=E/earth**2+M/moon**2+M/bm**2+S/sun**2+S/bs**2+j2
        snap=(E*(24*6**2/earth**4+2*A/earth**3)+M*(24*8**2/moon**4+2*(A+.00002)/moon**3)+
              M*(24*2**2/bm**4+2*.00002/bm**3)+S*(24*41**2/sun**4+2*(A+.00002)/sun**3)+
              S*(24*35**2/bs**4+2*.00002/bs**3))
        K=np.linalg.inv(I-h*h/6*J[k+1])
        T=np.empty((6,6))
        T[:3,:3]=K@(I+h*h/3*J[k]);T[:3,3:]=h*K
        T[3:,:3]=h/2*(J[k]+J[k+1]@T[:3,:3]);T[3:,3:]=I+h/2*J[k+1]@T[:3,3:]
        P1=T@P
        # Bound the distance between true position and the reference at both endpoints
        # before applying Taylor's theorem to the actual gravitational acceleration.
        D1=sz+h*(maxball(P[3:])+V)+h*h/2*A+np.linalg.norm(x[k+1,:3]-x[k,:3]-h*x[k,3:])+.00001
        actual=np.array([earth,moon,sun])
        ref=np.minimum(d[k,:3],d[k+1,:3])-.02
        floor=np.maximum(ref-D1,np.sqrt(np.maximum(0.,np.minimum(actual,ref)**2-D1**2/4)))
        if min(floor)<100:return 1e9,P[:3].copy()
        extra=np.sum(3.00001*mu*D1**2/floor**4)+j2+12*.001083*E*6379**2/ref[0]**4+1e-14
        extra+=.011*(np.sum(2*mu/floor**3)+2*M/bm**3+2*S/bs**3)
        defect_r=np.linalg.norm(x[k+1,:3]-x[k,:3]-h*x[k,3:]-h*h*(2*acc[k]+acc[k+1])/6)+1e-8
        defect_v=np.linalg.norm(x[k+1,3:]-x[k,3:]-h*(acc[k]+acc[k+1])/2)+1e-12
        # Uniform operator bound on the implicit endpoint solve.
        L1=np.sum(2*mu/ref**3);den=1-h*h*L1/6
        if den<=0:return 1e9,P[:3].copy()
        er=(h*h*extra/2+h**4*snap/24+defect_r)/den+1e-8
        ev=h*extra+h**3*snap/12+defect_v+h*L1/2*er+1e-12
        q=T@cov@T.T
        # Young's inequality encloses the sum of ellipsoids. The scaling only chooses
        # the free positive parameter and does not change the enclosure.
        sc=1000.;tr=np.trace(q[:3,:3])+sc**2*np.trace(q[3:,3:])
        p=max(1e-16,np.sqrt((6*er**2+6*(sc*ev)**2)/max(1e-40,tr)))
        cov=q*(1+p)*(1+1e-12)
        for j in range(6):cov[j,j]+=(1+1/p)*2*(er**2 if j<3 else ev**2)
        P=P1
        if np.trace(cov[:3,:3])>1e10:return 1e9,P[:3].copy()
    return np.sqrt(np.trace(cov[:3,:3])),P[:3].copy()
if __name__=='__main__':
    import affine_cells as ac
    from fast_reference import fasttrack
    out=[]
    for j in (0,11):
        last=path(np.array(nom['accepted'][j]['state']))[6*len(a.obs)-6:6*len(a.obs)];bundle=fasttrack(last)
        for f in (0.,.1,.25,.5,1.):
            w=np.array([100.,20.,8.,8.,.46,.46])*f
            ac.tube=tube;old=ac.test(bundle,w)
            ac.tube=ellipsoid;new=ac.test(bundle,w)
            out.append(dict(witness=j,scale=f,ball=old,ellipsoid=new))
    save('ellipsoid-test.json',out)
    print([(q['witness'],q['scale'],round(q['ball'].get('remainder_km',-1)),round(q['ellipsoid'].get('remainder_km',-1))) for q in out])
