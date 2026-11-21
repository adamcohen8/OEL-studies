"""Degree power and an orientation-independent bound for the fixed PDS gravity field."""
from hashlib import sha256
from common import *
def extract(path):
    power=np.zeros(1201);count=np.zeros(1201,dtype=int)
    with open(path,'rb') as fd:
        digest=sha256()
        header=fd.readline()
        digest.update(header)
        for line in fd:
            digest.update(line)
            deg,order,c,s,*_=line.split(b',')
            deg=int(deg)
            power[deg]+=float(c)**2+float(s)**2
            count[deg]+=1
    assert np.array_equal(count[1:],np.arange(2,1202))
    save("grail-power.json",dict(source_sha256=digest.hexdigest(),header=header.decode().strip(),power=power.tolist()))
from functools import lru_cache
@lru_cache(maxsize=512)
def envelope(r):
    if r<1737.4:return float("inf")
    n=np.arange(1,1201)
    return 4902.801/r**2*np.sum(factors*(1738/r)**n)+.0001/r**2
if __name__=="__main__":
    extract(sys.argv[1])
if (HERE/"grail-power.json").exists():
    power=np.array(json.loads((HERE/"grail-power.json").read_text())["power"])
    n=np.arange(1,1201)
    factors=1.000001*np.sqrt(power[1:]*(2*n+1)*(n+1)*(2*n+1))

@lru_cache(maxsize=64)
def potential_bound(r):
    n=np.arange(1,1201)
    return 4902.801/r*np.sum(1.000001*np.sqrt(power[1:]*(2*n+1))*(1738/r)**n)+.0001/r
