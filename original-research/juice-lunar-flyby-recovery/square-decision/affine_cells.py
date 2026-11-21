"""Pilot affine-state enclosure about an OEL trajectory, conditioned on lunar survival."""
from explore import *
from tubes import prepare,E,M,S
from copy import copy
from time import perf_counter
from numba import njit
from scipy.optimize import linprog
from outer_test import path
fc=D.parent/'field-certificate'
nom=json.loads((D/'nominal-certificate.json').read_text())['2']
a0=a.times[-2];tr=a.times[-1];site=a.sites[-1]
cfg=json.loads((D/'chart.json').read_text());transform=np.array(cfg['transform']);inverse=np.linalg.inv(transform)
shape=nom['geometry'];center=np.array(shape['center']);east,north=o.basis(center);roll=np.deg2rad(shape['roll_deg'])
axes=np.array([east*np.cos(roll)+north*np.sin(roll),-east*np.sin(roll)+north*np.cos(roll)])
faces=np.vstack([axes,-axes])-np.tan(np.deg2rad(2.2))*center;faces/=np.linalg.norm(faces,axis=1)[:,None]
def track(last):
    nominal=a.propagate(last,tr-a0,a0)
    te=tr-np.linalg.norm(nominal[:3]-site)/o.LIGHT_SPEED
    seg=copy(a);seg.times=np.array([0,te-a0]);seg.sites=a.sites[-2:]
    seg.JD=a.JD+a0/86400;seg.env=dict(a.env,jd_utc_start=seg.JD)
    raw=prepare(seg,last,.35)
    t=np.array([q[0] for q in raw]);x=np.array([q[1] for q in raw]);acc=np.array([q[2] for q in raw])
    d=np.array([q[3] for q in raw]);J=[]
    for time,xx in zip(t,x):
        body=o.hpop_de440_positions_km(seg.JD+time/86400,seg.ephem);der=np.zeros((3,3))
        for mu,b in [(E,0),(M,body['moon']),(S,body['sun'])]:
            pos=xx[:3]-b;rad=np.linalg.norm(pos);der+=mu*(3*np.outer(pos,pos)/rad**5-np.eye(3)/rad**3)
        J.append(der)
    body=lambda time:o.hpop_de440_positions_km(a.JD+time/86400,a.ephem)['earth']
    ray=x[-1,:3]-site+(body(te)-body(tr))
    return t,x,acc,d,np.array(J),ray,(tr-te)*o.LIGHT_SPEED
@njit
def maxball(B):
    return np.linalg.norm(np.sum(np.abs(B),axis=1))
@njit
def tube(t,x,acc,d,J,B, sample_ids, sample_P, sample_R):
    P=B[:3].copy();W=B[3:].copy();R=1e-7;V=1e-11
    idx=0
    for k in range(len(t)-1):
        if idx<len(sample_ids) and k==sample_ids[idx]:
            sample_P[idx]=P;sample_R[idx]=R;idx+=1
        h=t[k+1]-t[k];width=maxball(P)
        if d[k,1]+width+R+.02<1737.4:return -1.,P
        y=d[k]-np.array([width+R+6*h,width+R+8*h,width+R+41*h,2*h,35*h])-.02
        earth=max(6378.,y[0]);moon=max(1737.4,y[1]);sun=max(1.0e8,y[2])
        bm=y[3];bs=y[4]
        if bm<3e5 or bs<1e8:return 1e9,P
        # Outward error for two-body spatial derivatives at each numerical reference.
        I=np.eye(3);P1=np.linalg.solve(I-h*h*J[k+1]/6,P+h*W+h*h/3*J[k]@P)
        W1=W+h/2*(J[k]@P+J[k+1]@P1)
        dp=maxball(P1-P-h*W-h*h/6*(2*J[k]@P+J[k+1]@P1))
        dv=maxball(W1-W-h/2*(J[k]@P+J[k+1]@P1))
        j2=12*.001083*E*6379**2/earth**4
        A=E/earth**2+M/moon**2+M/bm**2+S/sun**2+S/bs**2+j2
        snap=(E*(24*6**2/earth**4+2*A/earth**3)+M*(24*8**2/moon**4+2*(A+.00002)/moon**3)+
              M*(24*2**2/bm**4+2*.00002/bm**3)+S*(24*41**2/sun**4+2*(A+.00002)/sun**3)+
              S*(24*35**2/bs**4+2*.00002/bs**3))
        largest=max(width,maxball(P1))+.0001
        floors=np.minimum(d[k,:3],d[k+1,:3])-largest-.02
        if min(floors)<50:return 1e9,P
        L=2*(E/min(floors[0],earth)**3+M/min(floors[1],moon)**3+S/min(floors[2],sun)**3)
        # Taylor's second derivative along a repeated displacement has norm <= 6 mu/r^4.
        nonlinear=3.00001*largest**2*(E/floors[0]**4+M/floors[1]**4+S/floors[2]**4)
        refj2=12*.001083*E*6379**2/min(d[k,0],d[k+1,0])**4
        extra=nonlinear+j2+refj2+.011*(L+2*M/bm**3+2*S/bs**3)+1e-14
        defect_r=np.linalg.norm(x[k+1,:3]-x[k,:3]-h*x[k,3:]-h*h*(2*acc[k]+acc[k+1])/6)+dp+1e-8
        defect_v=np.linalg.norm(x[k+1,3:]-x[k,3:]-h*(acc[k]+acc[k+1])/2)+dv+1e-12
        denom=1-h*h*L/6
        if denom<=0:return 1e9,P
        new=(R+h*V+h*h*(L*R/3+extra/2)+h**4*snap/24+defect_r)/denom
        V=(V+h*(L*(R+new)/2+extra)+h**3*snap/12+defect_v)*(1+1e-12)+1e-16
        R=new*(1+1e-12)+1e-12;P=P1;W=W1
        if R>1e5:return 1e9,P
    return R,P
def test(bundle,widths):
    t,x,b,d,j,ray,light=bundle
    R,P=tube(t,x,b,d,j,inverse@np.diag(widths),np.zeros(0,dtype=np.int64),np.empty((0,3,6)),np.empty(0))
    if R<0:return dict(status='all intersect surface')
    if R>=1e9:return dict(status='unresolved')
    delay=(abs(light-np.linalg.norm(ray))+maxball(P)+R+.021)/(o.LIGHT_SPEED-50)
    projected=faces@ray+np.sum(abs(faces@P),axis=1)+R+.021+50*delay
    return dict(status='contained' if max(projected)<0 else 'unresolved',remainder_km=R,
                max_outside_plane_km=float(max(projected)))
if __name__=='__main__':
    widths=np.array([100.,20.,8.,8.,.46,.46])
    out=[]
    for idx in [np.argmax([q['validation']['min_moon_floor_km'] for q in nom['accepted']]),np.argmin([q['validation']['min_moon_floor_km'] for q in nom['accepted']])]:
        row=nom['accepted'][idx];last=path(np.array(row['state']))[6*len(a.obs)-6:6*len(a.obs)]
        t0=perf_counter();bundle=track(last);secs=perf_counter()-t0
        results=[]
        for f in (0.,.02,.05,.1,.25,.5,1.,2.):
            start=perf_counter();v=test(bundle,widths*f);results.append(dict(scale=f,seconds=perf_counter()-start,**v))
        out.append(dict(witness=int(idx),track_seconds=secs,knots=len(bundle[0]),base_half_widths=widths.tolist(),tests=results))
        print(idx,[(x['scale'],x['status'],round(x.get('max_outside_plane_km',0),1)) for x in results],flush=True)
        save('affine-cells.json',out)
