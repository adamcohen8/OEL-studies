"""Higher-accuracy OEL reference for cells remaining near the survival boundary."""
from fast_reference import *
import fast_reference as fr
from quadratic import qr_tube
from types import FunctionType
@lru_cache(maxsize=1)
def fine_grid():
    tt=[a0]
    for stop,step in ((tr-13*3600,90.),(tr-10.75*3600,5.),(tr-7.0*3600,1.25),(tr-5*3600,5.),(tr-1,18.)):
        tt.extend(np.linspace(tt[-1],stop,int(np.ceil((stop-tt[-1])/step))+1)[1:])
    t=np.array(tt);body=[o.hpop_de440_positions_km(a.JD+s/86400,a.ephem) for s in t]
    return t,np.array([b['moon'] for b in body]),np.array([b['sun'] for b in body]),body[-1]['earth']-o.hpop_de440_positions_km(a.JD+tr/86400,a.ephem)['earth']
fine_track=FunctionType(fr.fasttrack.__code__,dict(fr.fasttrack.__globals__,grid=fine_grid))
def fine_tube(t,x,acc,d,J,B,ids,samp,err):
    _,m,s,_=fine_grid()
    return qr_tube(t,x,acc,d,J,B,ids,samp,err,m,s)
