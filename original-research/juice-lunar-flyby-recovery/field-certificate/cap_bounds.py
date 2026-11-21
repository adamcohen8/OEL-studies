from global_cover import *
from branch_cover import solve
a=Arc("juice");out=[]
for cap in (1.,3.):
    l,h,r=fit_cover(a,6,cap,0.,0.,0.)
    r.update(speed=6,cap=cap,timing=0,station_bias=0,extra_force=0.)
    r["cover"]=solve(a,r,24,12000.,directional=True,shift=1.)
    out.append(r);save("cap-bounds.json",out)
    print(cap,r["cover"]["diameter_deg"],flush=True)
