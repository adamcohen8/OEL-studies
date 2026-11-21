from pathlib import Path
import json
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
MISSIONS = HERE.parent / "followup/corrected/missions"
sys.path.insert(0,str(MISSIONS))
import forecast
from forecast import Arc, o
def load(name):
    return json.loads((MISSIONS / name).read_text())
def save(name, data):
    (HERE / name).write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
