from global_cover import *
from branch_cover import solve
a=Arc("juice")
out=[]
for timing,bias in [(1,0),(0,1)]:
    lo,hi,row=fit_cover(a,6,2,timing,bias,0.)
    row.update(speed=6,cap=2,timing=timing,station_bias=bias,extra_force=0.)
    row["post"]=solve(a,row,32,12000.)
    out.append(row);save("isolated-errors.json",out)
    print(timing,bias,row["post"]["diameter_deg"],flush=True)
