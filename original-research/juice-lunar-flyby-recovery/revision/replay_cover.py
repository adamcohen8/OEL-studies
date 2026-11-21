"""Replay the global partition and each terminal cell in a fresh destination."""
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
src=Path(__file__).resolve().parent.parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('destination',type=Path,help='new directory for derived results')
p.add_argument('--frozen-optics',action='store_true',help='start with the archived optical matrices')
p.add_argument('--oel',type=Path,required=True,help='installed OEL checkout root')
args=p.parse_args()
dst=args.destination.resolve()
if dst.exists():raise SystemExit(f'Destination already exists: {dst}')
dst.mkdir(parents=True)
for name in ('followup','field-certificate','square-decision'):
    shutil.copytree(src/name,dst/name,ignore=shutil.ignore_patterns('__pycache__','final-shard-*.json'))
work=dst/'square-decision'
env=dict(os.environ,PYTHONPATH=str(args.oel.resolve())+os.pathsep+os.environ.get('PYTHONPATH',''),
         OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',NUMBA_CACHE_DIR=str(dst/'numba-cache'))
def call(argv):
    subprocess.run([sys.executable,*argv],cwd=work,env=env,check=True)
if not args.frozen_optics:
    optical=dst/'field-certificate'
    subprocess.run([sys.executable,'-c',
        "from global_cover import Arc,fit_cover;fit_cover(Arc('juice'),6.,2.,0.,0.,0.,audit_tag='nominal')"],
        cwd=optical,env=env,check=True)
    call(['differential.py'])
    call(['contract.py'])
call(['audit_partition.py','optical-cover.json'])
with ThreadPoolExecutor(max_workers=2) as pool:
    list(pool.map(call,[['final_cover.py',str(i),'2'] for i in range(2)]))
call(['summarize.py'])
result=json.loads((work/'RESULTS.json').read_text())
print(json.dumps({k:result[k] for k in ['unresolved','final_classification','worst_contained_face_upper_km']},indent=2))
