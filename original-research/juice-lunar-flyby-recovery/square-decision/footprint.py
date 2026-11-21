"""Fit the ideal square and bound the smallest possible side for a finite set."""
from explore import *
from scipy.optimize import minimize
def fit(points,rad_as=None):
    from scipy.spatial.distance import pdist,squareform
    points=np.array(points)
    i,j=np.unravel_index(np.argmax(squareform(pdist(points))), (len(points),len(points)))
    mid=o.unit(points[i]+points[j]);east,north=o.basis(mid)
    def faces(p):
        center=o.unit(mid+p[0]*east+p[1]*north)
        e,n=o.basis(center)
        roll=p[2]
        axes=np.array([np.cos(roll)*e+np.sin(roll)*n,-np.sin(roll)*e+np.cos(roll)*n])
        return np.arctan2(points@axes.T,(points@center)[:,None]),center
    def slack(p):
        f,_=faces(p)
        return np.r_[p[3]-f.ravel(),p[3]+f.ravel()]
    opt=[]
    for roll in np.linspace(0,np.pi/2,8,endpoint=False):
        z=minimize(lambda p:p[3],[0.,0.,roll,.05],constraints=dict(type='ineq',fun=slack),method='SLSQP',options=dict(ftol=1e-12,maxiter=100))
        f,c=faces(z.x);opt.append((abs(f).max(),z.x.tolist(),c.tolist()))
    half,p,c=min(opt);ra=np.zeros(len(points)) if rad_as is None else np.deg2rad(np.array(rad_as)/3600)
    # For the 4.4 deg decision, every possible boresight lies within the detector
    # circumradius of every included direction. The nearest individual point to mid is enough.
    limit=np.deg2rad(2.2)
    tested_half=max(limit,half+2*max(ra)+1e-4)
    circ=np.arctan(np.sqrt(2)*np.tan(tested_half))
    # The sum of the endpoint unit rays points along mid. A feasible center has
    # dot product at least cos(circ + error) with each nominal endpoint.
    delta=np.arccos(np.clip((np.cos(circ+ra[i])+np.cos(circ+ra[j]))/np.linalg.norm(points[i]+points[j]),-1,1))
    delta+=1e-9
    # Parallel transport from mid to the boresight maps any pair of orthogonal
    # tangent axes to the true axes. The in-plane error on a ray difference is
    # bounded by (1-cos(delta))*|tangent difference| + sin(delta)*|normal difference|.
    e=points@east;n=points@north;z=points@mid
    diameter=np.linalg.norm([np.ptp(e),np.ptp(n)])
    correction=diameter*(1-np.cos(delta))+np.ptp(z)*np.sin(delta)+2*max(ra)
    rolls=np.linspace(0,np.pi/2,8192,endpoint=False);best=1
    for batch in np.array_split(rolls,32):
        C=np.cos(batch);S=np.sin(batch)
        b=np.maximum(np.ptp(e[:,None]*C+n[:,None]*S,axis=0),np.ptp(-e[:,None]*S+n[:,None]*C,axis=0))
        best=min(best,min(b))
    lower=best-correction-diameter*(np.pi/2/8192)/2-1e-12
    # Covering sides above the tested width satisfy the resulting lower bound directly.
    # Below the tested width the preceding orientation and pointing bounds apply.
    lower_side=min(2*tested_half,2*np.arcsin(lower/2))
    return dict(fit_side_deg=float(2*np.rad2deg(half)),center=c,roll_deg=float(np.rad2deg(p[2])%90),
                lower_test_side_deg=float(np.rad2deg(lower_side)),excludes_4p4=bool(lower_side>2*limit),
                max_axis_correction=correction,mesh_angles=len(rolls))
if __name__=='__main__':
    pts=json.loads((D/'exploration.json').read_text())['points']
    print(fit([q['los'] for q in pts]))
