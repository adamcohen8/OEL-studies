"""Post-fit check of the TESS September 9/10 optical errors against Horizons."""
import csv
from common import *
from comparison import second_tess
sys.path.insert(0,str(HERE.parent/"followup"))
import reference
reference.HERE=HERE
arc_key=sys.argv[1] if len(sys.argv)>1 else "september"
prefix="tess" if arc_key=="september" else "tess-"+arc_key
a,score=second_tess(arc_key)
assert (HERE/(prefix+"-comparison.json")).exists()
out=[]
for site in ("703","G96"):
    rows=[r for r in a.rows+[score] if r[1]==site]
    name=prefix+"-check-"+site
    path=HERE/(name+".txt")
    if not path.exists():
        meta=reference.fetch(name,COMMAND=-95,CENTER=site,EPHEM_TYPE="OBSERVER",CSV_FORMAT="YES",CAL_FORMAT="JD",TIME_DIGITS="FRACSEC",ANG_FORMAT="DEG",EXTRA_PREC="YES",QUANTITIES=1,OBJ_DATA="NO",TLIST="','".join(f"{2440587.5+r[0]/86400:.12f}" for r in rows))
        save(name+"-request.json",meta)
    table=path.read_text().split("$$SOE")[1].split("$$EOE")[0].strip()
    for r,fields in zip(rows,csv.reader(table.splitlines()),strict=True):
        ra,dec=np.deg2rad([float(v) for v in fields[3:5]])
        assert abs(float(fields[0])-(2440587.5+r[0]/86400))*86400<.001
        out.append(dict(row=r[3],site=site,score=r[3]==score[3],error_as=o.distance(r[2],[np.cos(ra)*np.cos(dec),np.sin(ra)*np.cos(dec),np.sin(dec)])))
save(prefix+"-horizons.json",out)
