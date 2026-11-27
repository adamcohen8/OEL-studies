"""Independent supervisor checks against the published true-anomaly equation."""
import math
import numpy as np
from scipy.integrate import solve_ivp
from spin_bridge import SpinBridge

def check():
    maxima = dict(tide=0., gravity_gradient=0., phase=0.)
    for a, inertia in [(57909050., 1.), (149597870.7, 1e25)]:
        model = SpinBridge(a_km=a, C_kg_m2=inertia)
        e, n = model.eccentricity, model.mean_motion_rad_s
        for E in np.linspace(0., 2*math.pi, 51, endpoint=False):
            M = E-e*np.sin(E)
            f = math.atan2(math.sqrt(1-e**2)*math.sin(E), math.cos(E)-e)
            rho = 1-e*math.cos(E)
            fd = math.sqrt(1-e**2)/rho**2
            s = model.orbit_sample(M/n)
            maxima['phase'] = max(maxima['phase'], abs(math.remainder(s.true_anomaly_rad-f,math.tau)))
            np.testing.assert_allclose(s.f_dot_rad_s/n,fd,rtol=3e-13)
            for spin in [fd-0.02,fd+0.07]:
                # Eq. 34 in true anomaly. Transform the bridge acceleration back for comparison.
                value = model.tide_accel_rad_s2(s,n*spin)/(fd*n)**2
                expected = -model.D*np.sign(spin/fd-1)/(rho**2*(1-e**2))
                maxima['tide'] = max(maxima['tide'],abs(value-expected))
            assert model.tide_accel_rad_s2(s,s.f_dot_rad_s)==0
            for theta in [-0.32, 0.67, 2.91]:
                analytic = -1.5*model.triaxiality*math.sin(2*(theta-f))/rho**3
                value = model.conservative_accel_rad_s2(s,theta)/n**2
                maxima['gravity_gradient'] = max(maxima['gravity_gradient'],abs(value-analytic))
            assert SpinBridge(tide_enabled=False).tide_accel_rad_s2(s,2*n)==0
    assert maxima['tide']<1e-15 and maxima['gravity_gradient']<3e-14 and maxima['phase']<1e-12
    # Independently step the full anomaly-domain Eq. 34 and the time-domain bridge.
    model=SpinBridge()
    n=model.mean_motion_rad_s
    e=model.eccentricity
    fd0=math.sqrt(1-e**2)/(1-e)**2
    def by_anomaly(f,y):
        psi,prime=y
        rho=(1-e**2)/(1+e*math.cos(f))
        torque=-model.D*math.copysign(1.,prime) if prime else 0.
        return prime, (2*rho*e*math.sin(f)*(1+prime)-1.5*model.triaxiality*rho*math.sin(2*psi)+torque/rho**2)/(1-e**2)
    def by_time(M,y):
        return y[1],model.angular_accel_rad_s2(M/n,y[0],y[1]*n)/n**2
    kw=dict(method='DOP853',rtol=2e-12,atol=2e-14,max_step=.04)
    theta,spin=.32,1.509
    end=math.tau
    ya=solve_ivp(by_anomaly,(0,end),(theta,spin/fd0-1),**kw)
    yt=solve_ivp(by_time,(0,end),(theta,spin),**kw)
    assert ya.success and yt.success
    delta=yt.y[:,-1]-[ya.y[0,-1]+end,(1+ya.y[1,-1])*fd0]
    assert max(abs(delta))<3e-8
    return dict(max_kernel_errors=maxima,one_orbit_angle_rad=float(delta[0]),one_orbit_spin_over_n=float(delta[1]))

if __name__=='__main__':
    import json
    print(json.dumps(check(),indent=2))
