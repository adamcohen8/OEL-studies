"""Compare the baseline methods at the same 24 additional curve evaluations around the common fit."""
from common import *
from comparison import second_tess
from itertools import product
from time import perf_counter
out=[]
for case in ("juice","tess","tess-november"):
    a=Arc("juice") if case=="juice" else second_tess("november" if case.endswith("november") else "september")[0]
    f=json.loads((HERE/(case+"-comparison.json")).read_text())
    x=np.array(f["baseline_state"]);n=len(a.rows)
    unit=np.eye(6);rate=1e-3
    start=perf_counter()
    def residual(y):
        return a.residual(y).ravel()
    jac=np.column_stack([(residual(x+a.delta(e))-residual(x-a.delta(e)))/(2*rate) for e in rate*unit])
    center=np.array(f["center"]);E,N=o.basis(center)
    def pointing(y):
        ray=a.sky(y,n+1)[-1];return np.rad2deg([ray@E,ray@N])
    pj=np.column_stack([(pointing(x+a.delta(e))-pointing(x-a.delta(e)))/(2*rate) for e in rate*unit])
    width=2*np.sqrt(5.991465)*np.linalg.svd(pj@np.linalg.pinv(jac),compute_uv=False)[0]
    cov_seconds=perf_counter()-start
    start=perf_counter();stars=[]
    # Two opposite radial neighbors one fine-grid spacing from the training fit. Retain the
    # short-arc angular slope by adjusting the transverse velocity with the range change.
    for sign in (-1,1):
        dr=sign*(4000/6 if case=="juice" else 60000/6)
        y=x+np.r_[dr*a.obs[0],dr*a.angular_slope]
        for it in range(2):
            r=residual(y)
            j=np.column_stack([(residual(y+a.delta(e*rate))-r)/rate for e in unit[[1,2,4,5]]])
            z=np.zeros(6);z[[1,2,4,5]]=np.linalg.lstsq(j,-r,rcond=1e-12)[0];y+=a.delta(z)
        rr=max(o.distance(a.sky(y,n),a.obs))
        u=a.sky(y,n+1)[-1]
        stars.append(dict(max_as=rr,los=u.tolist()))
    grid_seconds=perf_counter()-start
    out.append(dict(case=case,new_curve_calls_each=24,cov_seconds=cov_seconds,grid_seconds=grid_seconds,
                    gaussian95_diameter_deg=width,trial_max_residuals_as=[q["max_as"] for q in stars],
                    trial_pair_deg=o.distance(*[np.array(q["los"]) for q in stars])/3600,below_2as=int(sum(q["max_as"]<=2 for q in stars))))
save("equal-cost.json",out)
