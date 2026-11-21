"""Score frozen forecasts, compare mission-state replay, and verify witness geometry."""
import csv
import json
import numpy as np
from datetime import datetime, timezone
import forecast as f
from forecast import o, HERE

def table(name):
    data=(HERE/(name+".txt")).read_text().split("$$SOE")[1].split("$$EOE")[0].strip()
    return list(csv.reader(data.splitlines()))
def direction(ra,dec):
    ra,dec=np.deg2rad([ra,dec])
    return np.array([np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])
def score(target):
    a=f.Arc(target)
    frozen=json.loads((HERE/(target+"-forecast.json")).read_text())
    pred=np.array(frozen["los"])
    sky=table(target+"-reference-sky")[0]
    reference=direction(float(sky[3]),float(sky[4]))
    vectors=table(target+"-vectors")
    traj=np.array([list(map(float,row[2:8])) for row in vectors])
    times=np.array([(float(row[0])-a.JD)*86400 for row in vectors])
    rec=traj[0].copy()
    def integrate_early(x):
        return a.sky(x,len(a.rows)+1)
    geocentric=[np.linalg.norm(x[:3]) for x in traj]
    lunar=[np.linalg.norm(x[:3]-o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)["moon"]) for t,x in zip(times,traj)]
    control=integrate_early(rec)
    control_errors=o.distance(control[:-1],a.obs)
    a.engine=o.OrbitPropagator(integrator="dopri5",plugins=a.engine.plugins,adaptive_rtol=5e-13,adaptive_atol=3e-12)
    replays=[]
    for row in [dict(state=frozen['state'])]+frozen.get("witnesses",[]):
        x=np.array(row["state"])
        u=a.sky(x,len(a.rows)+1)
        body=o.hpop_de440_positions_km(a.JD+a.times[-1]/86400,a.ephem)
        moon=body["moon"]-a.sites[-1]
        replays.append(dict(training_max_as=max(o.distance(u[:-1],a.obs)),los=u[-1].tolist(),
                       limb_clearance_deg=o.distance(u[-1],o.unit(moon))/3600-np.rad2deg(np.arcsin(1737.4/np.linalg.norm(moon))),
                       elevation_deg=np.rad2deg(np.arcsin(u[-1]@o.unit(a.sites[-1]))),
                       forecast_error_deg=o.distance(u[-1],reference)/3600))
    result=dict(horizons_ra_deg=float(sky[3]),horizons_dec_deg=float(sky[4]),
               forecast_error_deg=o.distance(pred,reference)/3600,mission_initial_position_difference_km=np.linalg.norm(rec[:3]-np.array(frozen["state"][:3])),
               mission_initial_velocity_difference_m_s=1000*np.linalg.norm(rec[3:]-np.array(frozen["state"][3:])),
               same_model_mission_seed_forecast_error_as=o.distance(control[-1],reference),
               same_model_mission_seed_training_max_as=max(control_errors),
               reference_sampled_coast=dict(min_earth_km=min(geocentric),max_earth_km=max(geocentric),
                                            min_moon_km=min(lunar),max_speed_km_s=max(np.linalg.norm(traj[:,3:],axis=1))),
               verified_pair_diameter_deg=o.distance(*[np.array(v["los"]) for v in replays[-2:]])/3600,
               replay=replays)
    if target=="tess":
        first=next(r for r in f.read_rows(target) if r[1]=="703" and datetime.fromtimestamp(r[0],timezone.utc).day==9)
        result.update(optical_score_epoch=datetime.fromtimestamp(first[0],timezone.utc).isoformat(),
                      optical_forecast_error_deg=o.distance(pred,first[2])/3600,
                      optical_vs_horizons_as=o.distance(reference,first[2]))
    (HERE/(target+"-score.json")).write_text(json.dumps(result,indent=2)+"\n")
    print(target,{k:v for k,v in result.items() if k!="replay"},flush=True)
if __name__=="__main__":
    for name in ("tess","juice"): score(name)
