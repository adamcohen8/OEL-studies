"""What pure angular geometry plus a speed limit alone can establish at the first reception."""
from common import *
from scipy.optimize import linprog
a=Arc("juice")
s=a.sites[:-1];u=a.obs;t=a.times[:-1];n=len(u)
out=[]
for speed in (6.,10.):
    k=(speed+35)/o.LIGHT_SPEED+2/o.ASEC+1e-7
    D=speed*t[-1]+13000
    sine=np.linalg.norm(np.cross(u[0],u[-1]))
    first=(D+k*D/(1-k))/(sine-k-k*(1+k)/(1-k))
    universal=(first*(1+k)+D)/(1-k)+100
    rows=[];rhs=[]
    for i in range(n):
        for j in range(i):
            h=t[i]-t[j];directions=[*np.eye(3),*o.basis(u[i]),u[i],*o.basis(u[j]),u[j]]
            for axis in directions:
                for sign in (-1,1):
                    z=sign*axis
                    row=np.zeros(n);row[i]=z@u[i]-k;row[j]=-z@u[j]-k
                    rows.append(row);rhs.append(speed*h-z@(s[i]-s[j])+.05)
    bounds=[(0,universal)]*n;c=np.eye(n)[0]
    fits=[linprog(sign*c,A_ub=rows,b_ub=rhs,bounds=bounds,method='highs') for sign in (1,-1)]
    assert all(f.success for f in fits)
    out.append(dict(speed=speed,first_range_extrema_of_relaxation_km=[q.x[0] for q in fits],analytic_max_km=universal))
save("range-only.json",out)
