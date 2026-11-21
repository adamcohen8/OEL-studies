"""Regenerate fixed-state, modified-state and integration checks with the bundled maps."""
import json,sys,subprocess
from pathlib import Path
D=Path(__file__).resolve().parent
out=Path(sys.argv[1]).resolve()
out.mkdir(parents=True,exist_ok=False)
for lab in ('A','B','C','A+.003'):
    for deg,step,method,suffix in [(0,10,'rkf78','0'),(256,5,'rkf78','256-refined'),
                                   (1200,10,'rkf78','1200'),(1200,5,'dopri5','1200-refined')]:
        args=[sys.executable,str(D/'lunar_replay.py'),lab,str(deg),str(step),method]
        data=json.loads(subprocess.check_output(args))
        (out/f'{lab}-{suffix}.json').write_text(json.dumps(data,indent=2)+'\n')
    print(f'Replayed {lab}',flush=True)
args=[sys.executable,str(D/'check_geometry.py'),'--lunar-results',str(out)]
(out/'geometry.json').write_bytes(subprocess.check_output(args))
