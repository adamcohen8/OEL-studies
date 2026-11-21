"""Continuous derivative/radius budgets and lunar central-acceleration residual on the Juice arc."""
from numpy.polynomial import chebyshev as ch
from common import *
from tubes import E
from sim.dynamics.orbit.de440_hpop import jd_utc_to_jd_tdb
a=Arc("juice")
p=np.load(a.ephem["de440_coeff_path"])
start=jd_utc_to_jd_tdb(a.JD-.02,tai_utc_s=37);end=jd_utc_to_jd_tdb(a.JD+a.times[-1]/86400+.02,tai_utc_s=37)
tables={}
for body in ("moon","sun","earthmoon"):
    records=[]
    for k,(lo,hi) in enumerate(zip(p[body+"_start_jd_tdb"],p[body+"_end_jd_tdb"])):
        if hi<start or lo>end:continue
        coeff=np.stack([p[body+"_"+ax][k] for ax in "xyz"],axis=-1)
        pol=[coeff]
        for _ in range(3):
            pol.append(ch.chebder(pol[-1],axis=0)*2/((hi-lo)*86400))
        records.append((lo,hi,pol))
    tables[body]=records
def get(body,jd,order):
    lo,hi,pol=next(q for q in tables[body] if q[0]<=jd<=q[1])
    tau=2*(jd-lo)/(hi-lo)-1
    return ch.chebval(tau,pol[order])
def rate(body,k):
    return max(np.linalg.norm(np.sum(abs(segment[2][k]),axis=0)) for segment in tables[body])
def sample(jd):
    m=get("moon",jd,0);v=get("moon",jd,1);b=get("moon",jd,2)
    vals=[get("sun",jd,k)-get("earthmoon",jd,k)+get("moon",jd,k)/82.300568221497 for k in range(3)]
    return [np.linalg.norm(m),np.linalg.norm(v),np.linalg.norm(b),np.linalg.norm(vals[0]),np.linalg.norm(vals[1]),np.linalg.norm(vals[2]),np.linalg.norm(b+E*m/np.linalg.norm(m)**3)]
grid=np.linspace(start,end,int(np.ceil((end-start)*86400/60))+1)
r=np.array([sample(t) for t in grid])
mx=np.max(r,axis=0);mn=np.min(r,axis=0);dv=31
sun_rates=[sum([rate("sun",k),rate("earthmoon",k),rate("moon",k)/82.3]) for k in (1,2,3)]
moon_r=mn[0]-dv*rate("moon",1)-1
# Bound offset of any time to sampled time, allowing both sides of record joins.
out=dict(moon_floor_km=moon_r,moon_ceiling_km=mx[0]+dv*rate("moon",1)+1,
         moon_speed_kms=mx[1]+dv*rate("moon",2)+1e-6,moon_acc_kms2=mx[2]+dv*rate("moon",3)+1e-10,
         sun_floor_km=mn[3]-dv*sun_rates[0]-1,sun_speed_kms=mx[4]+dv*sun_rates[1]+1e-5,
         sun_acc_kms2=mx[5]+dv*sun_rates[2]+1e-9,
         moon_central_residual_kms2=mx[6]+dv*(rate("moon",3)+2*E*rate("moon",1)/moon_r**3)+1e-9)
out["max_record_join_differences"]={key:[max([float(np.linalg.norm(ch.chebval(1,left[2][k])-ch.chebval(-1,right[2][k]))) for left,right in zip(rows,rows[1:])]+[0]) for k in range(3)] for key,rows in tables.items()}
save("ephemeris.json",out);print(out)
