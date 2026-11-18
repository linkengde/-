#!/workspace/venvs/dft-labels/bin/python
"""Static run-0 screening for the corrected v12f potentials on the 347k model."""
from pathlib import Path
import importlib.util

MODEL = Path(__file__).resolve().parent
ROOT = MODEL.parents[1]
BASE_SCRIPT = ROOT / 'models/final_v10/run_large_static_preflight.py'
POD_COEFF = MODEL / 'POD/pod_v12f_coefficients.pod'
POD_DESC = MODEL / 'POD/AgTiSiC_POD_v12f.pod'
SNAP_COEFF = MODEL / 'SNAP/AgTiSiC_SNAP_v12f.snapcoeff'
SNAP_PARAM = MODEL / 'SNAP/AgTiSiC_SNAP_v12f.snapparam'
for path in (POD_COEFF, POD_DESC, SNAP_COEFF, SNAP_PARAM):
    if not path.is_file():
        raise FileNotFoundError(path)
spec = importlib.util.spec_from_file_location('static_helper_v12f', BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
base.MODEL = MODEL
base.OUT = MODEL / 'full_static_run0'
original = base.input_text
def input_text(name, folder):
    text = original(name, folder)
    text = text.replace(str(MODEL / 'POD/AgTiSiC_POD.pod'), str(POD_DESC))
    text = text.replace(str(MODEL / 'POD/pod_coefficients.pod'), str(POD_COEFF))
    text = text.replace('AgTiSiC_SNAP_final_v10.snapcoeff', SNAP_COEFF.name)
    text = text.replace('AgTiSiC_SNAP_final_v10.snapparam', SNAP_PARAM.name)
    return text
base.input_text = input_text
base.main()
