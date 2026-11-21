"""Verify all four exported linear supports by exact arithmetic on the stored doubles."""
from fractions import Fraction as F
from pathlib import Path
import json
import numpy as np
D=Path(__file__).resolve().parent/'results'
p=np.load(D/'worked-linear.npz')
faces=json.loads((D/'worked-cell.json').read_text())['faces']
for k,(g,y,scale,row) in enumerate(zip(p['goals'],p['multipliers'],p['scales'],faces,strict=True)):
    r=list(map(F,g));upper=F(0)
    assert np.all(y>=0)
    # Load CSR arrays once, including bounds used to account for every dual residual.
    b,l,h,ptr,idx,val=[p[n] for n in ('b','low','high','indptr','indices','values')]
    for i in np.flatnonzero(y):
        w=F(y[i]);upper+=w*F(b[i])
        for j in range(ptr[i],ptr[i+1]):r[idx[j]]-=w*F(val[j])
    upper+=sum(max(a*F(lo),a*F(hi)) for a,lo,hi in zip(r,l,h,strict=True))
    assert upper*F(scale)<=F(row['dual_upper_km'])
    total=sum(row[i] for i in ('reference_km','remainder_km','site_ephemeris_km','light_time_km','rounding_km','dual_upper_km'))
    assert abs(total-row['total_km'])<1e-8 and total<0
    print(f'Face {k+1}: rational support {float(upper*F(scale)):.8f} km; signed bound {total:.8f} km')
