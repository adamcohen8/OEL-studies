"""Reviewer's post-submission check of alternative attraction-error definitions."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from cube_study import LEVELS, errors

# Lower positive-y crossing of the 1% right-hand lobe, read from Figure 15.
# x, y are in cube side lengths. Digitization detail is in the README.
READINGS = [(0.8,0.684), (1.0,0.695), (1.2,0.640), (1.4,0.501)]


def crossing(x, index):
    low, high = 0.0, 0.8
    if not errors(x, low)[index] > 1 > errors(x, high)[index]:
        raise ValueError("Crossing is outside the selected interval")
    for _ in range(45):
        mid = (low + high) / 2
        if errors(x, mid)[index] > 1:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def main():
    out = Path("results")
    out.mkdir(exist_ok=True)
    readings = [dict(x=x, paper_y=y, magnitude_difference_y=crossing(x,0),
                     vector_difference_y=crossing(x,1)) for x,y in READINGS]
    (out/"figure_check.json").write_text(json.dumps(readings, indent=2)+"\n")
    axis = np.linspace(0, 2, 201)
    quadrant = np.full((len(axis),len(axis)), np.nan)
    for i,x in enumerate(axis):
        if x >= 0.5:
            for j in range(i+1):
                quadrant[i,j] = quadrant[j,i] = errors(x,axis[j])[1]
    idx = np.abs(np.arange(-200,201))
    xy = np.linspace(-2,2,401)
    colors = ("#313695","#4575b4","#d95f02","#a50026")
    fig,ax = plt.subplots(figsize=(7,6),layout="constrained")
    ax.contour(xy,xy,quadrant[np.ix_(idx,idx)],levels=LEVELS,colors=colors,corner_mask=False)
    ax.add_patch(plt.Rectangle((-.5,-.5),1,1,facecolor="0.88",edgecolor="black"))
    ax.set(xlim=(-2,2),ylim=(-2,2),aspect="equal",xlabel="x / side",ylabel="y / side",
           title="Relative acceleration-vector error")
    ax.legend(handles=[Line2D([],[],color=c,label=f"{v:g}%") for v,c in zip(LEVELS,colors)],
              title="Vector error",loc="upper left",bbox_to_anchor=(1.02,1))
    fig.savefig(out/"vector_comparison.png",dpi=150)
    plt.close(fig)
    print(json.dumps(readings,indent=2))


if __name__=="__main__":
    main()
