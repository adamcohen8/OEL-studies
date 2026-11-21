"""Replay selected global certificates against the final saved incoming boxes."""
from branch_cover import *
from time import perf_counter
a=Arc("juice");rows=json.loads((HERE/"global.json").read_text());out=[]
for i,(N,radius,shift,direct) in enumerate([(24,12000,1.,True),(24,12000,1.,True),(32,12000,0.,False),(32,12000,0.,False)]):
    start=perf_counter()
    q=solve(a,rows[i],N,radius,physical=i>=2,directional=direct,shift=shift)
    out.append(dict(global_row=i,elapsed_seconds=perf_counter()-start,**q))
    save("final-bounds.json",out);print(i,q["diameter_deg"],flush=True)
