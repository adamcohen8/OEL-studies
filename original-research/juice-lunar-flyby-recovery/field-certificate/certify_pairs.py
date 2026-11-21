"""Enclose individual emission directions for every witness used in the final diameter test."""
from tubes import *
from comparison import second_tess
def certify(a,x,cap,clock=None,bias=None,stations=None):
    n=len(a.rows)
    if clock is None: clock=np.zeros(n);bias=np.zeros((n,2))
    else:
        ids=[stations.index(r[1]) for r in a.rows]
        clock=np.asarray(clock)[ids];bias=np.asarray(bias)[ids]
    recv=a.times+np.r_[clock,0.]
    sites=np.array([a.site(r[1],t) for r,t in zip(a.rows,recv[:-1])]+[a.sites[-1]])
    x=np.array(x)
    emit=[]
    for t,site in zip(recv,sites):
        # A trial epoch suffices. The light-cone defect below bounds the error in the root.
        pos=a.propagate(x,t,0)[:3] if t>0 else x[:3]+t*x[3:]
        emit.append(t-np.linalg.norm(pos-site)/o.LIGHT_SPEED)
    data=prepare(a,x,.25,[t for t in emit if t>0])
    history={entry[0]:entry[1] for entry in data}
    boxes={}
    status=enclosure(data,a,0.,0.,record=boxes)
    if status["status"]!="bounded":return dict(status="unresolved",reason=status)
    def earth(t):return o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)["earth"]
    biasray=(bias[:,0,None]*a.east+bias[:,1,None]*a.north)/o.ASEC
    mag=np.linalg.norm(biasray,axis=1)[:,None]
    obs=np.cos(mag)*a.obs+np.sinc(mag/np.pi)*biasray
    errs=[];radii=[];last_ray=None
    for j,(t,te,site) in enumerate(zip(recv,emit,sites)):
        if te<=0:
            pos=x[:3]+te*x[3:];R=.5*.00002*te**2+1e-7
        else:
            pos=history[te][:3];R=boxes[te][0]
        ray=pos-site+(earth(te)-earth(t))
        dist=np.linalg.norm(ray)
        # Positive interpolation/rounding budget for both evaluations of barycentric Earth,
        # finite reference site evaluation and a 50 km/s barycentric speed envelope.
        eps=R+.0202
        delay=(abs((t-te)*o.LIGHT_SPEED-dist)+eps)/(o.LIGHT_SPEED-50)
        angle=np.rad2deg(np.arcsin((eps+50*delay)/dist))*3600+1e-6
        radii.append(angle);last_ray=o.unit(ray)
        if j<n:errs.append(o.distance(last_ray,obs[j])+angle)
    # Maximal geocentric speed over all step tubes checks both the physical guard and light-time speed bound.
    mx=max(v[4] for v in boxes.values());assert mx+35<50
    return dict(status="bounded",training_max_upper_as=max(errs),accepted=bool(max(errs)<=cap),
                forecast_center=last_ray.tolist(),forecast_radius_as=radii[-1],max_speed_upper=mx,
                min_moon_floor_km=min(v[3] for v in boxes.values()),min_earth_floor_km=min(v[2] for v in boxes.values()),
                training_moon_floor_km=min(v[3] for t,v in boxes.items() if t<=a.times[-2]+2),steps=len(data)-1)
def suite(name,a,cap,ws,stations=None):
    out=[]
    for w in ws:
        out.append(certify(a,w["state"],cap,w.get("station_clock_s"),w.get("station_bias_as"),stations))
    result=dict(rows=out)
    if all(q.get("accepted") for q in out):
        result["diameter_lower_deg"]=(o.distance(*[np.array(q["forecast_center"]) for q in out])-sum(q["forecast_radius_as"] for q in out))/3600
    print(name,result.get("diameter_lower_deg"),[r.get("training_max_upper_as") for r in out],flush=True)
    return result
if __name__=="__main__":
    a=Arc("juice")
    out={};full=json.loads((HERE/"expanded-witnesses.json").read_text())
    for k,r in full.items():
        out[k]=suite(k,a,2 if k=="shared" else float(k),r["witnesses"],r.get("stations"))
        save("certified-pairs.json",out)
    for k,r in json.loads((HERE/"tight-cap.json").read_text()).items():
        out["tight-"+k]=suite(k,a,float(k),r["witnesses"])
        save("certified-pairs.json",out)
    for name in ("tess","tess-november"):
        arc=second_tess("september" if name=="tess" else "november")[0]
        for k,r in json.loads((HERE/(name+"-comparison.json")).read_text())["grid"].items():
            out[name+"-"+k]=suite(name+k,arc,float(k),r["witnesses"]);save("certified-pairs.json",out)
