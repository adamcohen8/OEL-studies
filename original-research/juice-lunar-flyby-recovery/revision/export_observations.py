"""Print the 25 inputs with archived identifiers and received line-of-sight vectors."""
import csv,sys
from pathlib import Path
from datetime import datetime,timezone
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'field-certificate'))
from common import Arc
a=Arc('juice');writer=csv.writer(sys.stdout,lineterminator='\n')
writer.writerow(('archive_entry','observatory','reception_utc','ux','uy','uz','site_x_km','site_y_km','site_z_km'))
for (t,code,u,idx,flag),pos in zip(a.rows,a.sites):
    writer.writerow((idx,code,datetime.fromtimestamp(t,timezone.utc).isoformat(),*u,*pos))
