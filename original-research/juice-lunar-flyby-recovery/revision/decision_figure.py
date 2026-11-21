"""Draw the pointing decision directly from the frozen unit vectors."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
D=Path(__file__).resolve().parent.parent
def read(folder,name):return json.loads((D/folder/name).read_text())
def unit(u):return u/np.linalg.norm(u)
def basis(c):
    e=unit(np.cross([0,0,1],c))
    return e,np.cross(c,e)
def angle(a,b):return np.rad2deg(np.arctan2(np.linalg.norm(np.cross(a,b)),a@b))
res=read('square-decision','RESULTS.json')
rows=read('square-decision','nominal-certificate.json')['2']['accepted']
ls=np.array(read('followup/corrected/missions','juice-score.json')['replay'][0]['los'])
camera=np.array(res['boresight']);roll=np.deg2rad(res['roll_east_toward_north_deg'])
e,n=basis(camera);axes=np.array([np.cos(roll)*e+np.sin(roll)*n,-np.sin(roll)*e+np.cos(roll)*n])
def project(d):
    d=np.atleast_2d(d)
    return np.rad2deg(np.arctan2(d@axes.T,(d@camera)[:,None])).T
u=np.array([r['validation']['forecast_center'] for r in rows])
beta=np.arctan(np.sqrt(2)*np.tan(np.deg2rad(2.2)))
lo=angle(ls,u[5])-rows[5]['validation']['forecast_radius_as']/3600
phi=np.linspace(0,2*np.pi,501)
ax1,ax2=basis(ls)
ring=ls*np.cos(beta)+np.sin(beta)*(np.cos(phi[:,None])*ax1+np.sin(phi[:,None])*ax2)
blue='#214f78';rust='#a84d2d';green='#247463'
fig,(a,b)=plt.subplots(1,2,figsize=(9.15,4.1),gridspec_kw={'width_ratios':[1.12,1]})
fig.subplots_adjust(left=.08,right=.98,top=.86,bottom=.16,wspace=.36)
a.add_patch(Rectangle((-2.2,-2.2),4.4,4.4,color=green,alpha=.095))
a.add_patch(Rectangle((-2.2,-2.2),4.4,4.4,fill=False,color=green,lw=1.5))
a.plot(*project(ring),ls='--',lw=1,color=rust)
a.scatter(*project(u),s=14,c=blue,zorder=3)
a.scatter(*project(u[5]),s=52,facecolors='white',edgecolors=blue,zorder=4)
a.scatter(*project(ls),marker='*',s=90,c=rust,zorder=4)
a.scatter(0,0,marker='+',s=80,c=green,zorder=4)
a.set(xlim=(-3.42,3.42),ylim=(-3.42,3.42),aspect='equal',xlabel='Horizontal camera angle (deg)',ylabel='Vertical camera angle (deg)')
a.set_title('A   The two pointings',loc='left',weight='bold',fontsize=12)
b.set(xlim=(3.06,3.205),ylim=(-.1,1.08),xlabel='Angle from least-squares boresight (deg)',yticks=[])
b.set_title('B   A direction beyond every corner',loc='left',weight='bold',fontsize=12)
for x,y,col,lab in [(np.rad2deg(beta),.30,rust,'Farthest detector corner'),(lo,.78,blue,'Compatible orbit (lower bound)')]:
    b.hlines(y,3.06,x,color=col,lw=2);b.plot(x,y,'o',color=col)
    b.text(3.063,y+.085,lab,fontsize=10,color=col)
    b.text(x,y-.12,f'{x:.5f}°',ha='center',fontsize=10,color=col)
b.annotate('',(lo,.50),(np.rad2deg(beta),.50),arrowprops=dict(arrowstyle='<->',color='#42484f'))
b.text((lo+np.rad2deg(beta))/2,.52,'4.23 arcmin',ha='center',va='bottom',fontsize=10)
b.spines[['top','right','left']].set_visible(False)
a.spines[['top','right']].set_visible(False)
fig.savefig(D/'manuscript/detector-decision.pdf',bbox_inches='tight',pad_inches=.06)
