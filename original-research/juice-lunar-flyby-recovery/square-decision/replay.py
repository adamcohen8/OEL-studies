"""Replay stored terminal boxes, independently of the original subdivision schedule."""
from joint_cells import *
from quadratic import quadratic
from affine_cells import track as original_track, tube as ball_tube
from fast_reference import fasttrack
import affine_cells as acmod
from collections import Counter
def classify(row):
    box=np.array(row['box']);mid=box.mean(axis=1);width=np.ptp(box,axis=1)/2
    if row.get('joint'):
        out=eval_joint(box,local_optical=True)
        if out['status']!='unresolved':return out
        if not row.get('local_optical'):return eval_joint(box,local_optical=False)
        return out
    if row['depth']==0:
        bk=original_track(origin+inverse@mid)
        acmod.tube=ball_tube
        return acmod.test(bk,width)
    sc.tube=quadratic if row.get('quadratic') else ball_tube
    sc.track=fasttrack;return evaluate(box)
if __name__=='__main__':
    install();sc.support=CachedDual()
    group=json.loads((D/'optical-cover.json').read_text())['leaves']
    # Deterministic ordering allows an interrupted replay to resume at the next box.
    group=sorted(group,key=lambda r:r['id'])
    name='replayed-cells.json';output=json.loads((D/name).read_text()) if (D/name).is_file() else {}
    for k,row in enumerate(group):
        key=row['id']
        if key in output:continue
        q=classify(row)
        if q['status']=='unresolved' and not row.get('joint'):
            q=eval_joint(np.array(row['box']),local_optical=True)
        output[key]=q
        if k%15==0 or len(output)==len(group):
            save(name,output);print(k+1,Counter(v['status'] for v in output.values()),flush=True)
