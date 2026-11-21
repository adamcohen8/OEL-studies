"""Replay independently enlarged terminal boxes. Subdivide any inconclusive box without contraction."""
from replay import *
from hashlib import sha256
def certify_box(row,box,level):
    leaf=dict(row,box=box.tolist())
    if level:
        q=eval_joint(box,local_optical=True)
    else:
        q=classify(leaf)
        if q['status']=='unresolved' or (q['status']=='contained' and q.get('max_outside_plane_km',-100.)> -2):
            tighter=eval_joint(box,local_optical=True)
            if q['status']=='unresolved' or (tighter['status']!='unresolved' and
                   tighter.get('max_outside_plane_km',-1e99)<q.get('max_outside_plane_km',0)):q=tighter
    if q['status']!='unresolved' or level>=6:return dict(box=box.tolist(),result=q)
    axis=int(np.argmax(np.ptp(box,axis=1)*np.array([1.4,1.,.9,1.,.5,.5])))
    cut=float(box[axis].mean());left=box.copy();right=box.copy();left[axis,1]=cut;right[axis,0]=cut
    return dict(axis=axis,cut=cut,children=[certify_box(row,left,level+1),certify_box(row,right,level+1)])
def states(record):
    if 'children' in record:
        for child in record['children']:yield from states(child)
    else:yield record['result']
def main(shard,total):
    install();sc.support=CachedDual()
    bfile=D/'audited-boxes.json';digest=sha256(bfile.read_bytes()).hexdigest();boxes=json.loads(bfile.read_text())
    src={r['id']:r for r in json.loads((D/'optical-cover.json').read_text())['leaves']}
    name='final-shard-'+str(shard)+'.json'
    if (D/name).exists():
        packet=json.loads((D/name).read_text());assert packet['boxes_sha256']==digest and packet['shards']==total
    else:packet=dict(boxes_sha256=digest,shards=total,records={})
    rows=packet['records'];work=sorted(boxes)[shard::total]
    for key in work:
        if key in rows:continue
        rows[key]=certify_box(src[key],np.array(boxes[key]),0)
        if len(rows)%10==0 or len(rows)==len(work):
            tmp=D/(name+'.tmp');tmp.write_text(json.dumps(packet,indent=2,allow_nan=False)+'\n');tmp.replace(D/name)
            print(shard,len(rows),Counter(q['status'] for z in rows.values() for q in states(z)),flush=True)
if __name__=='__main__':main(int(sys.argv[1]),int(sys.argv[2]))
