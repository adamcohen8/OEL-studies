"""Construct encounter-oriented coordinates and bound all six widths with the optical program."""
from outer_replay import *
from scipy.linalg import orthogonal_procrustes
def build():
    pts=json.loads((D/'refined-test.json').read_text())['points']
    z=np.array([p['last'] for p in pts]);origin=ref[6*n-6:6*n];T=221000.
    # A linear change of variables only. This drift centers the coordinate system
    # near the encounter and does not replace the OEL propagation.
    drift=np.eye(6);drift[:3,3:]=np.eye(3)*T
    drift[3:,3:]*=1e4
    samples=(z-origin)@drift.T
    _,_,v=np.linalg.svd(samples,full_matrices=False)
    transform=v@drift
    A,bv=matrices();box=[]
    for row in transform:
        u=np.zeros(N);u[6*n-6:6*n]=row
        box.append([-support(A,bv,-u)-row@origin,support(A,bv,u)-row@origin])
    box=np.array(box)
    return dict(origin=origin.tolist(),transform=transform.tolist(),box=box.tolist(),drift_time_s=T)
if __name__=='__main__':
    info=build();save('chart.json',info);print('widths',np.ptp(info['box'],axis=1))
