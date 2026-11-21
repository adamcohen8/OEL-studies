from common import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
def read(name):return json.loads((HERE/(name+".json")).read_text())
cert=read("certified-pairs");outer=read("final-bounds");caps=read("cap-bounds")
fig,(left,right)=plt.subplots(1,2,figsize=(9.5,4.15),layout="constrained")
x=np.array([.65,1,1.25,1.5,2,3])
y=np.maximum.accumulate([cert[q]["diameter_lower_deg"] for q in ("tight-0.65","tight-1.0","tight-1.25","tight-1.5","2","3")])
left.step(x,y,where="post",color="#194b6d",linewidth=1.6)
left.scatter(x,y,color="#194b6d",s=27,zorder=3)
diag=np.rad2deg(2*np.arctan(np.sqrt(2)*np.tan(np.deg2rad(2.2))))
for yline,label in [(4.4,"Inscribed circle  4.40°"),(diag,"Square diagonal  6.22°")]:
    left.axhline(yline,ls="--",lw=1,color="#915127")
    left.text(2.04 if yline==4.4 else .61,yline+.06,label,color="#75401e",fontsize=8)
left.set(xlabel="Maximum angular error per observation (arcsec)",ylabel="Lower bound on field diameter (degrees)",xlim=(.58,3.16),ylim=(3.45,6.8))
left.set_title("Which fields can be ruled out?",loc="left",fontsize=12,fontweight="bold")
xs=[1,2,3];ys=[caps[0]["cover"]["diameter_deg"],outer[0]["diameter_deg"],caps[1]["cover"]["diameter_deg"]]
ls=[y[1],y[4],y[5]]
for xx,lo,hi in zip(xs,ls,ys):
    right.plot([xx,xx],[lo,hi],c="#8ba0ac",lw=3)
    right.scatter(xx,hi,c="#236959",s=36,zorder=3)
    right.scatter(xx,lo,c="#194b6d",s=36,zorder=3)
    right.annotate(f"{hi:.1f}°",(xx,hi),xytext=(0,6),textcoords="offset points",ha="center",fontsize=9)
    right.annotate(f"{lo:.2f}°",(xx,lo),xytext=(0,-12),textcoords="offset points",ha="center",fontsize=9)
right.scatter([],[],c="#236959",label="Sufficient diameter");right.scatter([],[],c="#194b6d",label="Required diameter ≥ this value")
right.set(xlabel="Maximum angular error per observation (arcsec)",ylabel="Field diameter (degrees)",xlim=(.45,3.55),ylim=(-1,78),xticks=xs)
right.legend(loc="upper left",bbox_to_anchor=(0.01,.99),fontsize=8,framealpha=.95)
right.set_title("How large is sufficient?",loc="left",fontsize=12,fontweight="bold")
for ax in (left,right):
    ax.spines[["right","top"]].set_visible(False);ax.grid(axis="y",alpha=.12)
fig.suptitle("Juice after the lunar flyby",x=.075,ha="left",fontsize=17,fontweight="bold")
for suffix in ("png","pdf"):
    fig.savefig(HERE/("juice-field-bounds."+suffix),dpi=160)
