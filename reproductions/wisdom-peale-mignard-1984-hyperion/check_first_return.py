"""Check the first native return of the production seeds with an independent planar model."""
import numpy as np
from scipy.integrate import solve_ivp

def reference(theta,spin,periods=1):
    # Dimensionless t=n*tphysical. Integrate eccentric anomaly with Kepler's time equation.
    e=.1
    def derivatives(t,y):
        E,theta,spin=y
        r=1-e*np.cos(E)
        f=np.arctan2(np.sqrt(1-e*e)*np.sin(E),np.cos(E)-e)
        return [1/r, spin, -0.5*.89**2*np.sin(2*(theta-f))/r**3]
    t=np.arange(periods+1)*2*np.pi
    y=solve_ivp(derivatives,(0,t[-1]),[0,theta,spin],t_eval=t,method='DOP853',rtol=2e-12,atol=2e-13)
    assert y.success
    return (y.y[1] % np.pi, y.y[2])

if __name__=="__main__":
    import json
    from pathlib import Path
    folder=Path(__file__).resolve().parent / "outputs" / "production-c16"
    errors=[]
    for row in json.loads((folder/"summary.json").read_text())["tracks"]:
        angles, spins=reference(*row["seed"])
        with np.load(folder/row["data"]) as data:
            first=data["peri"][1]
        errors.append([abs((first[0]-angles[1]+np.pi/2)%np.pi-np.pi/2),abs(first[1]-spins[1])])
    maximum=np.max(errors,axis=0)
    print(json.dumps({"tracks":len(errors),"max_angle_error_rad":maximum[0],"max_normalized_spin_error":maximum[1]}))
    if maximum[0]>6e-6 or maximum[1]>3e-6:
        raise SystemExit("First-return cross-check exceeded the recorded tolerance")
