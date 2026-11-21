"""Plot camera coordinates and the small set ruling out every detector roll."""
from explore import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from geometry_certificate import decision
from itertools import combinations
def read(f):return json.loads((D/f).read_text())
def coords(u,center,roll):
    e,n=o.basis(center);t=np.deg2rad(roll)
    return np.rad2deg(np.arctan2(u@np.array([np.cos(t)*e+np.sin(t)*n,-np.sin(t)*e+np.cos(t)*n]).T,(u@center)[:,None]))
fig,(left,right)=plt.subplots(1,2,figsize=(10.7,4.65),gridspec_kw={'width_ratios':[.92,1.16]})
fig.subplots_adjust(left=.075,right=.98,bottom=.18,top=.82,wspace=.34)
cam=read('outer-edges.json')['fit'];center=np.array(cam['center'])
u=np.array([r['validation']['forecast_center'] for r in read('nominal-certificate.json')['2']['accepted']])
z=coords(u,center,cam['roll_deg'])
left.add_patch(Rectangle((-2.2,-2.2),4.4,4.4,facecolor='#edf5f2',edgecolor='#277668',lw=1.5))
sampled=coords(np.array([k['los'] for k in read('exploration.json')['points']]),center,cam['roll_deg'])
left.scatter(sampled[:,0],sampled[:,1],s=9,color='#90aeb8',alpha=.72,label='Sampled orbits',zorder=2)
left.scatter(z[:,0],z[:,1],s=23,color='#214e76',linewidths=.5,edgecolors='white',label='Verified boundary orbits',zorder=3)
ls=np.array(load('juice-forecast.json')['los']);v=coords(ls[None,:],center,cam['roll_deg'])[0]
left.scatter(*v,s=66,marker='*',color='#b15a28',label='Least-squares prediction',zorder=4)
left.set(xlim=(-2.6,2.6),ylim=(-2.6,2.6),aspect='equal',xlabel='Camera horizontal angle (degrees)',ylabel='Camera vertical angle (degrees)')
left.set_title('2.00″ per observation',loc='left',fontweight='bold')
left.legend(loc='upper right',fontsize=8,framealpha=.97)
case=read('all-rolls.json')['2.15'];rows=read(case['source']+'.json')['2.15']['accepted']
u=np.array([rows[i]['validation']['forecast_center'] for i in case['indices']])
r=np.deg2rad([rows[i]['validation']['forecast_radius_as']/3600 for i in case['indices']]);pairs=np.array(list(combinations(range(3),2)));i,j=pairs.T
dif=u[i]-u[j];k=np.argmax(np.linalg.norm(dif,axis=1));mid=o.unit(u[i[k]]+u[j[k]]);e,n=o.basis(mid)
x,y,z=dif@e,dif@n,dif@mid
circ=np.arctan(np.sqrt(2)*np.tan(np.deg2rad(2.2)))
co=(np.cos(circ+r[i[k]])+np.cos(circ+r[j[k]]))/np.linalg.norm(u[i[k]]+u[j[k]])
de=np.arccos(co);err=(1-co)*np.hypot(x,y)+np.sin(de)*abs(z)+r[i]+r[j]
angles=np.linspace(0,90,20001);cc=np.cos(np.deg2rad(angles));ss=np.sin(np.deg2rad(angles))
w=np.max(np.maximum(abs(x[:,None]*cc+y[:,None]*ss),abs(-x[:,None]*ss+y[:,None]*cc))-err[:,None],axis=0)
limit=2*np.sin(np.deg2rad(2.2));delta=(w-limit)*1000
k=delta.argmin();ra=angles[k]
right.plot(angles,delta,color='#174c72',lw=2)
right.axhline(0,color='#a34b30',ls='--',lw=1.3,label='Largest width the detector allows')
right.scatter(ra,delta[k],s=32,color='#174c72')
right.annotate(f'Minimum +{delta[k]:.3f}',(ra,delta[k]),xytext=(ra-.37,delta[k]+.65),fontsize=9,arrowprops=dict(arrowstyle='-',color='#174c72'))
right.set(xlim=(ra-1.2,ra+1.2),ylim=(-.17,1.6),xlabel='Roll around midpoint direction (degrees)',ylabel='Excess projected width (×10⁻³)')
right.set_title('2.15″ per observation',loc='left',fontweight='bold')
right.legend(loc='upper center',fontsize=8,frameon=False);right.grid(axis='y',alpha=.2)
for ax in (left,right):ax.spines[['top','right']].set_visible(False)
fig.suptitle('Could one Catalina exposure cover Juice after the lunar flyby?',x=.075,ha='left',fontsize=15,fontweight='bold')
for ext in ('png','pdf'):fig.savefig(D/('detector-decision.'+ext),dpi=155)
