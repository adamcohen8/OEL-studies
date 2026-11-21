from common import *
from comparison import second_tess
# Use the archived independent emission-time implementation, redirecting its arc factory in memory.
from types import SimpleNamespace
import ast
source=(MISSIONS/"emission.py").read_text()
module=ast.parse(source);body=[n for n in module.body if isinstance(n,ast.FunctionDef) and n.name=="validate"]
fn=compile(ast.Module(body=body,type_ignores=[]),str(MISSIONS/"emission.py"),"exec")
out={}
for case in ("juice","tess","tess-november"):
    target="juice" if case=="juice" else "tess"
    arc=Arc("juice") if target=="juice" else second_tess("november" if case=="tess-november" else "september")[0]
    ns=dict(np=np,o=o,f=SimpleNamespace(Arc=lambda _:arc))
    exec(fn,ns)
    f=json.loads((HERE/(case+"-comparison.json")).read_text())
    result={}
    for cap,rec in f["grid"].items():
        rows=rec["witnesses"]
        check=ns["validate"](target,rows)
        geo=[];bounds=[]
        for w in rows:
            x=np.array(w["state"]);g=arc.geometry(x)
            if target=="juice":g["encounter"]=arc.encounter(x)
            geo.append(g)
            assert arc.physical(g)
        assert all(v["max_fit_as"]<float(cap) and min(v["min_elevation_deg"],v["min_limb_clearance_deg"])>0 for v in check["witnesses"])
        angle=check["lower_diameter_deg"]
        diagonal=lambda side:2*np.rad2deg(np.arctan(np.sqrt(2)*np.tan(np.deg2rad(side/2))))
        side=4.4 if case!="tess" else 10560*.77/3600
        check.update(geometry=geo,minimum_square_side_for_pair_deg=2*np.rad2deg(np.arctan(np.tan(np.deg2rad(angle/2))/np.sqrt(2))),
                    detector_side_deg=side,detector_diagonal_deg=diagonal(side),pair_excludes_square=bool(angle>diagonal(side)))
        result[cap]=check
    out[case]=result
save("independent-pairs.json",out)
