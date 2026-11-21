"""Transport bounded quadrature defects in separate signed linear-response blocks."""
from affine_cells import *
@njit
def block_tube(t,x,acc,d,J,B,sample_ids,sample_P,sample_R):
    count=0;starts=0.;Z=np.zeros((96,6,6));g=np.zeros((96,6));ident=np.eye(6)
    Z[0]=ident.copy();g[0,:3]=1e-7;g[0,3:]=1e-8
    P=B.copy();P[3:]*=1000.
    idx=0;I=np.eye(3);mus=np.array([E,M,S]);R=1e-6;V=1e-10
    for k in range(len(t)-1):
        if idx<len(sample_ids) and k==sample_ids[idx]:
            sample_P[idx]=P[:3];sample_R[idx]=R;idx+=1
        h=t[k+1]-t[k];width=maxball(P[:3]);sz=width+R+.00001
        if d[k,1]+sz+.02<1737.4:return -1.,P[:3].copy()
        y=d[k]-np.array([sz+6*h,sz+8*h,sz+41*h,2*h,35*h])-.02
        earth=max(6378.,y[0]);moon=max(1737.4,y[1]);sun=max(1.e8,y[2]);bm=y[3];bs=y[4]
        if bm<3e5 or bs<1e8:return 1e9,P[:3].copy()
        j2=12*.001083*E*6379**2/earth**4
        A=E/earth**2+M/moon**2+M/bm**2+S/sun**2+S/bs**2+j2
        snap=(E*(24*36/earth**4+2*A/earth**3)+M*(24*64/moon**4+2*(A+.00002)/moon**3)+
              M*(24*4/bm**4+2*.00002/bm**3)+S*(24*41**2/sun**4+2*(A+.00002)/sun**3)+
              S*(24*35**2/bs**4+2*.00002/bs**3))
        K=np.linalg.inv(I-h*h/6*J[k+1]);T=np.empty((6,6))
        T[:3,:3]=K@(I+h*h/3*J[k]);T[:3,3:]=(h/1000)*K
        T[3:,:3]=500*h*(J[k]+J[k+1]@T[:3,:3]);T[3:,3:]=I+500*h*J[k+1]@T[:3,3:]
        P1=T@P
        D1=sz+h*(maxball(P[3:])/1000+V)+h*h/2*A+np.linalg.norm(x[k+1,:3]-x[k,:3]-h*x[k,3:])+.00001
        actual=np.array([earth,moon,sun]);rr=np.minimum(d[k,:3],d[k+1,:3])-.02
        floor=np.maximum(rr-D1,np.sqrt(np.maximum(0.,np.minimum(actual,rr)**2-D1**2/4)))
        if min(floor)<100:return 1e9,P[:3].copy()
        extra=np.sum(3.00001*mus*D1**2/floor**4)+j2+12*.001083*E*6379**2/rr[0]**4+1e-14
        extra+=.011*(np.sum(2*mus/floor**3)+2*M/bm**3+2*S/bs**3)
        dr=np.linalg.norm(x[k+1,:3]-x[k,:3]-h*x[k,3:]-h*h*(2*acc[k]+acc[k+1])/6)+1e-8
        dv=np.linalg.norm(x[k+1,3:]-x[k,3:]-h*(acc[k]+acc[k+1])/2)+1e-12
        L1=np.sum(2*mus/rr**3);den=1-h*h*L1/6
        if den<=0:return 1e9,P[:3].copy()
        eps=1e-10*(1+R+width+h*(maxball(P[3:])/1000+V))
        er=(h*h*extra/2+h**4*snap/24+dr)/den+eps
        ev=1000*(h*extra+h**3*snap/12+dv+h*L1/2*er)+eps
        # Limit the duration of each disturbance block. Each older contribution
        # keeps its own transported axes and independent component widths.
        length=min(24000.,max(600.,d[k,1]/3.))
        if t[k]-starts>length and count<95:
            count+=1;starts=t[k];Z[count]=ident.copy()
        ext=np.zeros(6)
        for m in range(count+1):
            Z[m]=T@Z[m]
        invz=np.linalg.inv(Z[count]);res=np.max(np.sum(np.abs(ident-Z[count]@invz),axis=1))+1e-11
        if res>.001:return 1e9,P[:3].copy()
        invpad=np.max(np.sum(np.abs(invz),axis=1))*res/(1-res)
        for i in range(6):
            g[count,i]+=(np.linalg.norm(invz[i,:3])+2*invpad)*er+(np.linalg.norm(invz[i,3:])+2*invpad)*ev
        for m in range(count+1):
            for i in range(6):
                for j in range(6):
                    ext[i]+=abs(Z[m,i,j])*g[m,j]
        R=np.linalg.norm(ext[:3])+eps;V=(np.linalg.norm(ext[3:])+eps)/1000;P=P1
        if R>1e5:return 1e9,P[:3].copy()
    return R,P[:3].copy()
if __name__=='__main__':
    import affine_cells as ac
    from fast_reference import fasttrack
    out=[]
    for i in (0,11):
        last=path(np.array(nom['accepted'][i]['state']))[6*len(a.obs)-6:6*len(a.obs)];bundle=fasttrack(last)
        for f in (0,.1,.25,.5,1.):
            w=np.array([100.,20.,8.,8.,.46,.46])*f
            ac.tube=block_tube;q=ac.test(bundle,w)
            out.append(dict(witness=i,scale=f,**q))
    save('block-test.json',out);print([(r['witness'],r['scale'],r.get('remainder_km')) for r in out],flush=True)
