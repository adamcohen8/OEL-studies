"""Finish all outstanding branches, keeping each accepted, empty and unresolved leaf."""
from joint_cells import *
def main():
    install();sc.support=CachedDual()
    local='--optical' in sys.argv;file=D/('optical-cover.json' if local else 'joint-cover.json')
    if file.is_file():leaves=json.loads(file.read_text())['leaves']
    else:leaves=json.loads((D/('joint-cover.json' if local else 'quadratic-cover.json')).read_text())['leaves']
    mapping={r['id']:r for r in leaves}
    def record():
        data=list(mapping.values());save(file.name,dict(counts=dict(Counter(r['result']['status'] for r in data)),leaves=data))
    for level in range(int(sys.argv[1])):
        pending=[r for r in mapping.values() if r['result']['status']=='unresolved']
        if not pending:break
        first=any(not r.get('local_optical' if local else 'joint') for r in pending)
        depth=min(v['depth'] for v in pending)
        pending=[r for r in pending if (not r.get('local_optical' if local else 'joint') if first else r['depth']==depth)]
        for j,parent in enumerate(pending):
            bx=np.array(parent['box']);children=[]
            if first:cases=[(parent['id'],parent['depth'],bx)]
            else:
                width=np.ptp(bx,axis=1);k=np.argmax(width*np.array([1.4,1.,.9,1.,.5,.5]))
                cases=[]
                for s in (0,1):
                    box=bx.copy();box[k,1-s]=bx[k].mean()
                    cases.append((parent['id']+'.'+str(s),parent['depth']+1,box))
            for key,dep,b in cases:
                if not first and impossible(b):res=dict(status='optically empty')
                else:
                    if not first:b=tightbox(b)
                    res=eval_joint(b,local_optical=local)
                children.append(dict(id=key,depth=dep,box=b.tolist(),result=res,joint=True,local_optical=local))
            del mapping[parent['id']];mapping.update({r['id']:r for r in children})
            if j%8==0:record()
        record();print(level,Counter(r['result']['status'] for r in mapping.values()),flush=True)
if __name__=='__main__':main()
