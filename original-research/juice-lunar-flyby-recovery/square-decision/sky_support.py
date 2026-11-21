"""Search actual square-axis supports with multiple starts and an explicit periapsis clearance."""
from explore import *
from footprint import fit
from functools import lru_cache
history=json.loads((D/'exploration.json').read_text())
initial=history['points'];mid=o.unit(np.mean([q['los'] for q in initial],axis=0));east,north=o.basis(mid)
t0=(datetime(2024,8,19,21,tzinfo=timezone.utc).timestamp()-a.rows[0][0])
def moon(t):return o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)['moon']
@lru_cache(maxsize=1000)
def evaluate(params):
    z=np.array(params);x=state(z[:6]);t=t0+z[6]*100
    xx=a.propagate(x,t,0);q=xx[:3]-moon(t);vv=xx[3:]-(moon(t+1)-moon(t-1))/2
    return np.r_[a.residual(x).ravel(),a.sky(x,arc+1)[-1],np.linalg.norm(q)/100,o.unit(q)@vv]
def get(z):return evaluate(tuple(z))
def deriv(z):
    # Central perturbations remove adaptive-step and ephemeris epoch quantization noise.
    h=.006
    return np.column_stack([(get(z+d)-get(z-d))/(2*h) for d in np.eye(7)*h])
cap=1.986;floor=1742.
def cons(z):
    q=get(z);return np.r_[cap**2-np.sum(q[:-5].reshape(-1,2)**2,axis=1),q[-2]-floor/100]
def cjac(z):
    q=get(z);j=deriv(z);return np.vstack([-2*np.einsum('ij,ijk->ik',q[:-5].reshape(-1,2),j[:-5].reshape(-1,2,7)),j[-2]])
def seek(angle,seed):
    axis=east*np.cos(angle)+north*np.sin(angle)
    p=coords(np.array(seed['state']));ts=datetime.fromisoformat(seed['encounter']['utc']).timestamp()-a.rows[0][0]
    x=np.r_[p,(ts-t0)/100]
    opt=minimize(lambda z:-get(z)[-5:-2]@axis*57.3,x,jac=lambda z:-axis@deriv(z)[-5:-2]*57.3,
         constraints=[dict(type='ineq',fun=cons,jac=cjac),dict(type='eq',fun=lambda z:get(z)[-1],jac=lambda z:deriv(z)[-1])],
         method='SLSQP',bounds=[(-32,32)]*6+[(-150,150)],options=dict(ftol=1e-8,maxiter=32))
    if min(cons(opt.x))<-.015:return None
    q=data(opt.x[:6])
    return q if q and q['max_as']<cap+.005 else None
if __name__=='__main__':
    pts=initial[:];out=[]
    # Test each suggested roll first then the rest of the boundary.
    angles=[79,-11,169,259,35,215,0,90,180,270,55,235]+list(range(0,360,15))
    for deg in dict.fromkeys(angles):
        r=np.deg2rad(deg);axis=east*np.cos(r)+north*np.sin(r)
        seeds=sorted(pts,key=lambda q:np.array(q['los'])@axis,reverse=True)
        choose=[seeds[0]]
        # Start on both sides of the narrow family even for directions transverse to it.
        if abs(np.cos(r-np.deg2rad(35)))<.6:
            choose.append(max(initial,key=lambda q:q['encounter']['moon_radius_km']) if seeds[0]['encounter']['moon_radius_km']<3300 else min(initial,key=lambda q:q['encounter']['moon_radius_km']))
        for seed in choose:
            q=seek(r,seed)
            if q:out.append(q);pts.append(q)
        result=fit([q['los'] for q in pts])
        save('sky-support.json',dict(points=out,fit=result,angle=deg))
        print(deg,len(out),result['fit_side_deg'],result['lower_test_side_deg'],flush=True)
        if result['lower_test_side_deg']>4.418:break
