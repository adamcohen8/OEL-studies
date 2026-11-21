"""Resolve outstanding terminal cells with quadratic transport and surface constraints."""
from adaptive import *
from quadratic import quadratic
from warm_support import CachedDual
def main():
    install();sc.tube=quadratic;sc.support=CachedDual()
    name='quadratic-cover.json';pathout=D/name
    if pathout.is_file():initial=json.loads(pathout.read_text())['leaves']
    else:initial=json.loads((D/'adaptive.json').read_text())['leaves']
    mapping={p['id']:p for p in initial}
    def persist():
        leaves=list(mapping.values())
        save(name,dict(algorithm='affine or quadratic with signed detector-face and lunar-survival constraints',counts=dict(Counter(r['result']['status'] for r in leaves)),leaves=leaves))
    for dep in range(int(sys.argv[1])):
        todo=[p for p in mapping.values() if p['result']['status']=='unresolved']
        if not todo:break
        # First retest every inherited unresolved cell without splitting. Then visit a breadth layer.
        first=dep==0 and any(not p.get('quadratic') for p in todo)
        if not first:
            depth=min(p['depth'] for p in todo);todo=[p for p in todo if p['depth']==depth]
        for j,p in enumerate(todo):
            b=np.array(p['box'])
            if first:
                if p.get('quadratic'):continue
                tests=[(p['id'],p['depth'],b)]
            else:
                wid=np.ptp(b,axis=1);k=np.argmax(wid*np.array([1.4,1.,.9,1.,.5,.5]));mid=b[k].mean()
                tests=[]
                for s in range(2):
                    bx=b.copy();bx[k,1-s]=mid
                    tests.append((p['id']+'.'+str(s),p['depth']+1,bx))
            new=[]
            for key,depth,bx in tests:
                if first:box=bx;v=evaluate(box)
                else:box,v=run_one(bx)
                new.append(dict(id=key,depth=depth,box=box.tolist(),result=v,quadratic=True))
            del mapping[p['id']];mapping.update({v['id']:v for v in new})
            if j%8==0:persist()
        persist();print('cycle',dep,Counter(q['result']['status'] for q in mapping.values()),flush=True)
if __name__=='__main__':main()
