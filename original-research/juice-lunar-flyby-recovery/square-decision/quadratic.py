"""Second-order initial-state enclosure with a cubic gravitational remainder."""
from fast_reference import *
@njit
def hessian_forcing(pos,moon,sun,P):
    out=np.zeros((3,21))
    for m in range(3):
        v=pos.copy()
        if m==1:v-=moon
        if m==2:v-=sun
        r=np.sqrt(v@v);v/=r;w=v@P
        coeff=3*(E if m==0 else M if m==1 else S)/r**4
        k=0
        for i in range(6):
            for j in range(i,6):
                uv=P[0,i]*P[0,j]+P[1,i]*P[1,j]+P[2,i]*P[2,j]
                fac=coeff*(.5 if i==j else 1.)
                for idx in range(3):out[idx,k]+=fac*(w[i]*P[idx,j]+w[j]*P[idx,i]+(uv-5*w[i]*w[j])*v[idx])
                k+=1
    return out
@njit
def qr_tube(t,x,acc,d,J,B,ids,samp,err_rec,moons,suns):
    P=B[:3].copy();W=B[3:].copy();Q=np.zeros((3,21));U=Q.copy()
    R=1e-7;V=1e-11;idx=0;I=np.eye(3);mus=np.array([E,M,S])
    f0=hessian_forcing(x[0,:3],moons[0],suns[0],P)
    for k in range(len(t)-1):
        h=t[k+1]-t[k];lin=maxball(P);quad=maxball(Q);sz=lin+quad+R+.00001
        if idx<len(ids) and k==ids[idx]:
            if samp.shape[2]==27:
                samp[idx,:,:6]=P;samp[idx,:,6:]=Q;err_rec[idx]=R
            else:
                samp[idx]=P;err_rec[idx]=R+quad
            idx+=1
        if d[k,1]+sz+.02<1737.4:return -1.,np.hstack((P,Q))
        y=d[k]-np.array([sz+6*h,sz+8*h,sz+41*h,2*h,35*h])-.02
        earth=max(6378.,y[0]);moon=max(1737.4,y[1]);sun=max(1e8,y[2]);bm=y[3];bs=y[4]
        if bm<3e5 or bs<1e8:return 1e9,np.hstack((P,Q))
        j2=12*.001083*E*6379**2/earth**4
        A=E/earth**2+M/moon**2+M/bm**2+S/sun**2+S/bs**2+j2
        snap=(E*(6.00001*36/earth**4+2*A/earth**3)+M*(6.00001*64/moon**4+2*(A+.00002)/moon**3)+
              M*(6.00001*4/bm**4+2*.00002/bm**3)+S*(6.00001*41**2/sun**4+2*(A+.00002)/sun**3)+
              S*(6.00001*35**2/bs**4+2*.00002/bs**3))
        K=np.linalg.inv(I-h*h*J[k+1]/6);T=I+h*h*J[k]/3
        P1=K@(T@P+h*W);W1=W+h/2*(J[k]@P+J[k+1]@P1)
        f1=hessian_forcing(x[k+1,:3],moons[k+1],suns[k+1],P1)
        Q1=K@(T@Q+h*U+h*h/6*(2*f0+f1))
        U1=U+h/2*(J[k]@Q+J[k+1]@Q1+f0+f1)
        lin=max(lin,maxball(P1))+.00001;quad=max(quad,maxball(Q1))+.00001
        rr=np.minimum(d[k,:3],d[k+1,:3])-.02;floor=rr-lin-quad
        if min(floor)<100:return 1e9,np.hstack((P,Q))
        L=np.sum(2*mus/np.minimum(floor,np.array([earth,moon,sun]))**3)
        # For a unit direction with cosine z to r, 16 minus the squared norm of the
        # cubic coefficient equals 5(1-z^2)(35z^4+2z^2+11)/4 in units mu/r^5.
        extra=np.sum(4.00001*mus*(lin+quad)**3/floor**5+3.00001*mus*(2*lin*quad+quad**2)/rr**4)
        extra+=j2+12*.001083*E*6379**2/rr[0]**4+.011*(L+2*M/bm**3+2*S/bs**3)+1e-14
        den=1-h*h*L/6
        if den<=0:return 1e9,np.hstack((P,Q))
        eps=1e-11*(1+lin+quad+h*(maxball(W)+maxball(U)))
        dr=np.linalg.norm(x[k+1,:3]-x[k,:3]-h*x[k,3:]-h*h*(2*acc[k]+acc[k+1])/6)+1e-8+eps
        dv=np.linalg.norm(x[k+1,3:]-x[k,3:]-h*(acc[k]+acc[k+1])/2)+1e-12+eps/100
        new=(R+h*V+h*h*(L*R/3+extra/2)+h**4*snap/24+dr)/den
        V=(V+h*(L*(R+new)/2+extra)+h**3*snap/12+dv)*(1+1e-12)+1e-16
        R=new*(1+1e-12)+1e-12;P=P1;W=W1;Q=Q1;U=U1;f0=f1
        if R>1e5:return 1e9,np.hstack((P,Q))
    return R,np.hstack((P,Q))
def quadratic(t,x,acc,d,J,B,ids,samp,err):
    _,moons,suns,_=grid()
    return qr_tube(t,x,acc,d,J,B,ids,samp,err,moons,suns)
def upper_polynomial(q):
    diag=np.array([0,6,11,15,18,20]);return sum(abs(q))-sum(abs(q[diag]))+sum(np.maximum(0,q[diag]))
if __name__=='__main__':
    out=[]
    for j in (0,11):
        last=path(np.array(nom['accepted'][j]['state']))[6*len(a.obs)-6:6*len(a.obs)];bundle=fasttrack(last)
        for f in (0.,.1,.25,.5,1.):
            w=np.array([100.,20.,8.,8.,.46,.46])*f
            R,P=quadratic(*bundle[:5],inverse@np.diag(w),np.zeros(0,dtype=int),np.empty((0,3,6)),np.empty(0))
            out.append(dict(witness=j,scale=f,remainder_km=R,quad_contribution_km=maxball(P[:,6:])))
    save('quadratic-test.json',out)
    print([(q['witness'],q['scale'],round(q['remainder_km'],1),round(q['quad_contribution_km'],1)) for q in out],flush=True)
