"""Extract the used coefficients verbatim from NASA products. Check hashes in inputs/sources.json."""
from pathlib import Path
import hashlib, shutil, sys, json
import numpy as np
from itertools import islice
raw=Path(sys.argv[1])
out=Path(__file__).resolve().parent/'inputs'
out.mkdir(exist_ok=True)
sources={}
def source(name,url):
    p=raw/name
    sources[name]=dict(url=url,sha256=hashlib.file_digest(p.open('rb'),'sha256').hexdigest())
    return p
pds='https://pds-geosciences.wustl.edu/lro/lro-l-lola-3-rdr-v1/lrolol_1xxx/data/lola_gdr/cylindrical/img/'
naif='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/'
jpl='https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de430/'
grail='https://pds-geosciences.wustl.edu/grail/grail-l-lgrs-5-rdr-v1/grail_1001/shadr/'
for name,url in [('ldem_16.img',pds),('ldem_16.lbl',pds),('ldem_128.lbl',pds),('moon_080317.tf',naif+'fk/satellites/'),('moon_pa_de421_1900-2050.bpc',naif+'pck/'),('naif0012.tls',naif+'lsk/'),('gggrx_1200a_sha.lbl',grail),('header.430_229',jpl)]:
    shutil.copyfile(source(name,url+name),out/name)
shutil.copyfile(raw/'gggrx_1200a_sha.tab',out/'grgm1200a.tab')
with source('gggrx_1200a_sha.tab',grail+'gggrx_1200a_sha.tab').open('rb') as f, (out/'grail256.tab').open('wb') as dest:
    dest.write(next(f))
    for line in f:
        if int(line.split(b',')[0])>256:break
        dest.write(line)
with source('ascp1950.430',jpl+'ascp1950.430').open() as f:
    for header in f:
        _,count=map(int,header.split())
        chunk=''.join(islice(f,(count+2)//3))
        vals=np.fromstring(chunk.replace('D','E'),sep=' ')[:count]
        if vals[0]<=2460537 and vals[1]>=2460544:
            np.savez(out/'libration.npz',start=vals[0],end=vals[1],angles=vals[898:1018].reshape(4,3,10));break
    else:raise ValueError('missing dates')
(out/'sources.json').write_text(json.dumps(sources,indent=2)+'\n')
