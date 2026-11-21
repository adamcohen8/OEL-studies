from tubes import *
from gravity_model import envelope
a=Arc("juice");out=[]
for i,w in enumerate(load("juice-clearance-support.json")["independent_replay"]):
    data=prepare(a,w["state"],.5)
    for rad in (.001,.1,1.,10.):
        res=enclosure(data,a,rad,5e-10,envelope)
        out.append(dict(witness=i,initial_position_radius_km=rad,**res))
        print(i,rad,res.get("field_radius_deg"),res.get("at_hours"),flush=True)
save("physical-tubes.json",out)
