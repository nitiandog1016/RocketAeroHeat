from pathlib import Path
import sys
root=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(root/'scripts'))
from runSU2RocketPhaseV222 import supervise,BASE,DOC,ROOT
locks=[Path(__file__).resolve(),ROOT/'scripts/finalizeSU2RocketPhaseV222.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',DOC/'AUDIT_CONTRACT.md']
raise SystemExit(0 if supervise('final',[sys.executable,str(ROOT/'scripts/finalizeSU2RocketPhaseV222.py')],30,256,locks,BASE,'offline') else 1)
