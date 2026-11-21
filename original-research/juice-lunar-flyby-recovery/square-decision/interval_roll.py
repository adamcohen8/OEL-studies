"""Independent outward interval replay of all possible detector orientations for the compact witnesses."""
from pathlib import Path
from itertools import combinations
import json
from mpmath import iv
iv.dps=36
D=Path(__file__).resolve().parent
def dot(a,b):return sum(i*j for i,j in zip(a,b))
def norm(a):return iv.sqrt(dot(a,a))
def unit(a):return [v/norm(a) for v in a]
def cross(a,b):return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def run(group,ids,N=4096,side=4.4):
    # Normalize exact representations of the stored double centers. Include an extra
    # 1e-10 radian in their radii for representation and normalization.
    q=[group[i]['validation'] for i in ids]
    u=[unit([iv.mpf(v) for v in x['forecast_center']]) for x in q]
    r=[iv.mpf(x['forecast_radius_as'])*iv.pi/(180*3600)+iv.mpf('1e-10') for x in q]
    pairs=list(combinations(range(len(q)),2))
    i,j=min(pairs,key=lambda ij:float(dot(u[ij[0]],u[ij[1]]).a))
    add=[v+w for v,w in zip(u[i],u[j])];mid=unit(add)
    east=unit(cross([0,0,1],mid));north=cross(mid,east)
    t=iv.mpf(str(side))*iv.pi/360
    # cos(circ)=1/sqrt(1+2*tan(t)^2). Dot product with mid bounds the boresight.
    cc=1/iv.sqrt(1+2*iv.tan(t)**2);sc=iv.sqrt(1-cc**2)
    co=(cc*(iv.cos(r[i])+iv.cos(r[j]))-sc*(iv.sin(r[i])+iv.sin(r[j])))/norm(add)
    if float(co.a)>1:return dict(pair_excludes=True)
    st=iv.sqrt(1-co**2)
    vals=[]
    for i,j in pairs:
        d=[v-w for v,w in zip(u[i],u[j])]
        x,y,z=dot(d,east),dot(d,north),dot(d,mid)
        error=(1-co)*iv.sqrt(x**2+y**2)+st*abs(z)+r[i]+r[j]
        vals.append((x,y,error))
    worst=1
    for k in range(N):
        angle=iv.mpf([k,k+1])*iv.pi/(2*N)
        c,s=iv.cos(angle),iv.sin(angle)
        margins=[abs(v)-err-2*iv.sin(t) for x,y,err in vals for v in [x*c+y*s,-x*s+y*c]]
        bound=max(float(v.a) for v in margins)
        if bound<=0:return dict(success=False,unresolved_bin=k,margin=bound)
        worst=min(worst,bound)
    return dict(success=True,roll_intervals=N,precision_decimal_digits=iv.dps,minimum_width_margin=worst)
if __name__=='__main__':
    allroll=json.loads((D/'all-rolls.json').read_text());out={}
    for cap,row in allroll.items():
        if 'indices' not in row:continue
        g=json.loads((D/(row['source']+'.json')).read_text())[cap]['accepted']
        out[cap]=run(g,row['indices'])
    (D/'interval-roll.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out)
