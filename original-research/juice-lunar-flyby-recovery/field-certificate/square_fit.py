"""Place an ideal rectangular gnomonic field over the saved nominal Juice witness directions."""
from common import *
from scipy.optimize import minimize
a=Arc("juice")
def read(s):return json.loads((HERE/(s+".json")).read_text())
points=[load("juice-forecast.json")["los"]]
points.extend(w["los"] for w in read("expanded-witnesses")["2"]["witnesses"])
for q in read("tight-cap").values():points.extend(w["los"] for w in q["witnesses"])
for f in [load("juice-forecast.json")["witnesses"],load("juice-clearance-support.json")["independent_replay"],read("juice-comparison")["grid"]["2"]["witnesses"]]:
    points.extend(w["los"] for w in f)
points=np.unique(points,axis=0)
mid=o.unit(np.mean(points,axis=0));east,north=o.basis(mid)
def faces(p):
    center=o.unit(mid+p[0]*east+p[1]*north)
    e,n=o.basis(center)
    ar=np.deg2rad(p[2]);u=np.cos(ar)*e+np.sin(ar)*n;v=-np.sin(ar)*e+np.cos(ar)*n
    return np.rad2deg(np.arctan((points@np.array([u,v]).T)/(points@center)[:,None])),center
def obj(p):return np.max(np.abs(faces(p)[0]))
res=min([minimize(obj,[0.,0.,roll],method="Nelder-Mead",options=dict(maxiter=800,xatol=1e-11,fatol=1e-10)) for roll in range(0,91,15)],key=lambda z:z.fun)
coords,center=faces(res.x)
save("square-fit.json",dict(sampled_unique_directions=len(points),verified_covering_side_deg=2*np.max(abs(coords)),center=center.tolist(),
                   roll_deg_from_east_toward_north=res.x[2]%90,points=points.tolist(),coordinates_deg=coords.tolist(),
                   nominal_703_side_deg=4.4,diagonal_deg=2*np.rad2deg(np.arctan(np.sqrt(2)*np.tan(np.deg2rad(2.2))))))
print(len(points),2*np.max(abs(coords)))
