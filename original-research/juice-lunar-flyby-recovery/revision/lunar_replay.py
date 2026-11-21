"""Fixed-initial-state terrain and gravity sensitivity, using OEL dynamics throughout."""
import sys, json
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import spiceypy as sp
from scipy.interpolate import CubicHermiteSpline
D=Path(__file__).resolve().parent
F=D/'inputs'
sys.path.insert(0,str(D.parent/'field-certificate'))
from common import Arc,o
from sim.dynamics.orbit.propagator import accel_third_body, MOON_MU_KM3_S2
from sim.dynamics.orbit.spherical_harmonics import (SphericalHarmonicTerm,
    compile_spherical_harmonic_terms,_load_normalized_spherical_harmonic_kernel)

def read(folder,name):
    return json.loads((D.parent/folder/name).read_text())
def run(label,degree,knot,integrator='rkf78'):
    a=Arc('juice')
    for name in ['naif0012.tls','moon_080317.tf','moon_pa_de421_1900-2050.bpc']:
        sp.furnsh(str(F/name))
    et0=sp.str2et('2024-08-15T07:27:53.856')
    lib=np.load(F/'libration.npz')
    start=float(lib['start']);span=(float(lib['end'])-start)/4
    coeff=lib['angles']
    def pa(t):
        sub=(2451545.-start+(et0+t)/86400)/span
        k=int(np.floor(sub))
        assert 0<=k<4
        phi,theta,psi=np.polynomial.chebyshev.chebval(2*(sub-k)-1,coeff[k].T)
        return sp.eul2m(psi,theta,phi,3,1,3)
    def bodies(t):
        return o.hpop_de440_positions_km(a.JD+t/86400,a.ephem)
    rot_delta=sp.raxisa(pa(0)@sp.pxform('J2000','MOON_PA_DE421',et0).T)[1]*o.ASEC
    assert rot_delta<3,rot_delta
    if degree:
        terms=[]
        with (F/('grgm1200a.tab' if degree>256 else 'grail256.tab')).open() as f:
            radius,gm=map(float,next(f).split(',')[:2])
            for line in f:
                n,m,c,s,*_=line.split(',')
                if int(n)>degree:break
                if int(n)>=2:terms.append(SphericalHarmonicTerm(int(n),int(m),float(c),float(s),True))
        ct=compile_spherical_harmonic_terms(terms)
        kernel=_load_normalized_spherical_harmonic_kernel()
        args=(gm,radius,ct.c_nm,ct.s_nm,ct.legendre_diag_scale,ct.legendre_subdiag_scale,
              ct.legendre_recur_a,ct.legendre_recur_b,ct.legendre_recur_c)
        def harmonics(q,rot):
            n=degree
            if degree>256:
                ratio=np.linalg.norm(q)/radius
                n=32 if ratio>=2.5 else 128 if ratio>=1.3 else 256 if ratio>=1.12 else degree
            return kernel(q,rot,*args,n,n)
        def lunar(t,x,env,ctx):
            m=bodies(t)['moon']
            rot=pa(t)
            return (harmonics(x[:3]-m,rot)-harmonics(-m,rot)
                    +accel_third_body(x[:3],m,gm-MOON_MU_KM3_S2))
        a.engine=o.OrbitPropagator(integrator=integrator,numeric_backend='python',plugins=[*a.engine.plugins,lunar],
                                   adaptive_rtol=3e-13,adaptive_atol=3e-12)
    row=read('square-decision','threshold-certificate.json')['2.15']['accepted'][dict(A=1,B=3,C=4)[label[0]]]
    x0=np.array(row['state'])
    nominal=next(r for r in read('square-decision','near-threshold.json')['points'] if r['state']==row['state'])
    peri=datetime.fromisoformat(nominal['encounter']['utc']).timestamp()-a.rows[0][0]
    if label.startswith('A+'):
        alpha=float(label[2:])
        best=np.array(read('followup/corrected/missions','juice-forecast.json')['state'])
        x0=x0+alpha*(best-x0)
    knots=np.unique(np.r_[0.,a.times,np.maximum(a.times-10.,0.),np.arange(0,a.times[-1],1800.),
                          np.arange(peri-1800,peri+1800+knot/2,knot)])
    states=[x0]
    for left,right in zip(knots,knots[1:]):
        states.append(a.propagate(states[-1],right-left,left))
    states=np.array(states)
    curve=CubicHermiteSpline(knots,states[:,:3],states[:,3:])
    rays=[]
    for i,t in enumerate(a.times):
        lag=0.
        for _ in range(5):
            emit=t-lag
            if emit<0:
                acc=a.engine.acceleration_at(t_s=0.,x_eci=x0,command_accel_eci_km_s2=o.zero,env=a.env,ctx=o.ctx)
                pos=x0[:3]+emit*x0[3:]+.5*emit**2*acc
            else:pos=curve(emit)
            ray=(pos-a.sites[i])+(bodies(emit)['earth']-bodies(t)['earth'])
            lag=np.linalg.norm(ray)/o.LIGHT_SPEED
        rays.append(o.unit(ray))
    residual=o.distance(np.array(rays[:-1]),a.obs)
    dem=np.memmap(F/'ldem_16.img',mode='r',dtype='<i2',shape=(2880,5760))
    patch=np.load(F/'terrain128.npz');high=patch['grid'];r0,c0=int(patch['row0']),int(patch['col0'])
    def shape(t):
        q=sp.pxform('J2000','MOON_ME_DE421',et0+t)@(curve(t)-bodies(t)['moon'])
        r,lon,lat=sp.reclat(q)
        lon=np.rad2deg(lon)%360;lat=np.rad2deg(lat)
        col=lon*16-.5;line=(90-lat)*16-.5
        i,j=int(np.floor(line)),int(np.floor(col))
        dx,dy=col-j,line-i
        grid=dem[np.clip(np.array([i,i+1]),0,2879)[:,None],np.array([j,j+1])%5760]*.0005
        h=np.array([1-dy,dy])@grid@np.array([1-dx,dx])
        entry=[r-1737.4-h,r-1737.4,h,lon,lat,r-1737.4-grid.max(),t]
        if r-1737.4<14:
            row=(90-lat)*128-.5-r0;col=lon*128-.5-c0
            i,j=int(np.floor(row)),int(np.floor(col))
            assert 0<=i<high.shape[0]-1 and 0<=j<high.shape[1]-1,(lon,lat)
            g=high[i:i+2,j:j+2]*.0005
            rh=np.array([1-(row-i),row-i])@g@np.array([1-(col-j),col-j])
            return entry+[r-1737.4-rh,rh,r-1737.4-g.min(),r-1737.4-g.max()]
        return entry+[np.inf,np.nan,np.inf,np.inf]
    # The high terrain-risk section is sampled at 0.25 seconds, with 5 s shoulders.
    samples=np.unique(np.r_[np.arange(peri-1800,peri+1801,5.),np.arange(peri-220,peri+220.25,.25),np.arange(peri-85,peri+85.025,.025)])
    data=np.array([shape(t) for t in samples]);v=data[data[:,0].argmin()]
    # Coarser guard against a lower passage elsewhere in the full arc.
    outside=[np.linalg.norm(x[:3]-bodies(t)['moon']) for t,x in zip(knots,states) if abs(t-peri)>1799]
    high_result=None
    if np.isfinite(data[:,7]).any():
        u=data[data[:,7].argmin()]
        high_result=dict(min_clearance_km=u[7],terrain_km=u[8],sphere_height_km=u[1],
            lon_east_deg=u[3],lat_deg=u[4],t_s=u[6],clearance_above_lowest_corner_km=u[9],
            min_clearance_above_highest_corner_km=min(data[:,10]))
    return dict(initial_state=x0.tolist(),highres=high_result,witness=label,degree=degree,knot_s=knot,integrator=integrator,angle_PA430_to_PA421_epoch_arcsec=rot_delta,
                min_DEM_clearance_km=v[0],sphere_height_at_min_km=v[1],terrain_at_min_km=v[2],
                longitude_east_deg=v[3],latitude_deg=v[4],min_four_pixel_clearance_km=min(data[:,5]),
                t_min_s=v[6],min_sphere_height_km=min(data[:,1]),min_outer_knot_radius_km=min(outside),
                forecast_los=rays[-1].tolist(),shift_from_archived_center_as=o.distance(rays[-1],row['validation']['forecast_center']),
                max_sample_speed_km_s=max(np.linalg.norm(states[:,3:],axis=1)),min_sample_earth_km=min(np.linalg.norm(states[:,:3],axis=1)),
                min_training_moon_km=min(np.linalg.norm(curve(t)-bodies(t)['moon']) for t in knots if t<=a.times[-2]),
                max_residual_as=max(residual),residuals_as=residual.tolist())
if __name__=='__main__':
    lab,deg,step=sys.argv[1:4]
    result=run(lab,int(deg),float(step),*sys.argv[4:])
    print(json.dumps(result,indent=2))
