"""Recompute pointing and interval geometry with the frozen optical and replayed lunar rays."""
from pathlib import Path
import sys,json
import numpy as np
D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'square-decision'))
from interval_roll import run
def read(path):return json.loads(path.read_text())
def ang(a,b):return np.rad2deg(np.arctan2(np.linalg.norm(np.cross(a,b)),np.dot(a,b)))
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--lunar-results',type=Path,default=D/'results')
data=parser.parse_args().lunar_results
S=D.parent/'square-decision'
boresight=read(S/'RESULTS.json')['boresight']
ls=read(D.parent/'followup/corrected/missions/juice-score.json')['replay'][0]['los']
q=read(S/'nominal-certificate.json')['2']['accepted'][5]['validation']
corner=np.rad2deg(np.arctan(np.sqrt(2)*np.tan(np.deg2rad(2.2))))
low=ang(ls,q['forecast_center'])-q['forecast_radius_as']/3600
out=dict(corner_deg=corner,witness_distance_lower_deg=low,miss_arcmin=60*(low-corner),
    repointing_deg=ang(ls,boresight),nominal_witness_height_floor_km=q['min_moon_floor_km']-1737.4)
assert low>corner and q['training_max_upper_as']<2
j=read(S/'threshold-certificate.json')['2.15']['accepted']
out['nominal_geometry']=run(j,[1,3,4])
assert out['nominal_geometry']['success']
case={}
for lab in ('A','B','C','A+.003'):
    case[lab]=read(data/(lab+'-1200.json'))
for lab in ('A','A+.003'):
    obs=[{'validation':{'forecast_center':case[k]['forecast_los'],'forecast_radius_as':5.}} for k in [lab,'B','C']]
    out[lab+'_three_ray_geometry_with_5arcsec_padding']=run(obs,[0,1,2])
    assert max(case[k]['max_residual_as'] for k in [lab,'B','C'])<2.148
out['surviving_pair_separation_deg']=ang(case['B']['forecast_los'],case['C']['forecast_los'])
assert out['surviving_pair_separation_deg']<2*corner
assert case['A']['highres']['min_clearance_km']<0
for k in ['C','A+.003']:assert case[k]['highres']['min_clearance_km']>0
assert out['A+.003_three_ray_geometry_with_5arcsec_padding']['success']
print(json.dumps(out,indent=2))
