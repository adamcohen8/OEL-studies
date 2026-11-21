from geometry_certificate import *
from interval_roll import run
group=json.loads((D/'nominal-certificate.json').read_text())['2']['accepted']
ids=list(range(len(group)));width=4.313
def check(idx):
    q=[group[k]['validation'] for k in idx]
    return decision([r['forecast_center'] for r in q],[r['forecast_radius_as'] for r in q],width,32768)['excluded']
for i in ids[:]:
    sub=[j for j in ids if j!=i]
    if len(sub)>=2 and check(sub):ids=sub
out=run(group,ids,32768,width);out.update(side_deg=width,indices=ids)
(D/'nominal-lower-interval.json').write_text(json.dumps(out,indent=2)+'\n')
print(out)
