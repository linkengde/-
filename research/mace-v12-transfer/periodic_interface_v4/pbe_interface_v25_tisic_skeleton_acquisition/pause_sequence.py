#!/usr/bin/env python3
"""Request a checkpointed pause at the next 20-SCF-step boundary."""
from pathlib import Path
import json
from datetime import datetime, timezone
ROOT=Path(__file__).resolve().parent
(ROOT/'pause.request').write_text(json.dumps({'requested_utc':datetime.now(timezone.utc).isoformat(),'requested_by':'user'},indent=2)+'\n')
print('Pause requested; the active job will stop after writing a matching checkpoint at the next interval.')
