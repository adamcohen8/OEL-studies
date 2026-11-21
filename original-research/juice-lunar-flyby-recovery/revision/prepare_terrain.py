"""Extract one georeferenced LOLA subimage using an HTTP range and preserve its provenance."""
import hashlib,json,urllib.request
from pathlib import Path
import numpy as np
out=Path(__file__).resolve().parent/'inputs'
url='https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/ldem_128.img'
r0,r1=95*128,113*128
c0,c1=90*128,120*128
stride=46080*2
first,last=r0*stride,r1*stride-1
req=urllib.request.Request(url,headers={'Range':f'bytes={first}-{last}'})
with urllib.request.urlopen(req,timeout=120) as resp:
    assert resp.headers['Content-Range']==f'bytes {first}-{last}/2123366400'
    data=resp.read()
assert len(data)==(r1-r0)*stride
crop=np.frombuffer(data,dtype='<i2').reshape(-1,46080)[:,c0:c1].copy()
np.savez_compressed(out/'terrain128.npz',grid=crop,row0=r0,col0=c0)
(out/'terrain128-source.json').write_text(json.dumps(dict(url=url,first_byte=first,last_byte=last,
    band_sha256=hashlib.sha256(data).hexdigest(),row_range=[r0,r1],column_range=[c0,c1]),indent=2)+'\n')
