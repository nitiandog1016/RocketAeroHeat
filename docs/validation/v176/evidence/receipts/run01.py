from pathlib import Path
import sys,os
sys.path.insert(0,str(Path.cwd()/'scripts'))
from validateIdealAir import run
out=Path('outputs/v176_validation')
run(out,'main01_supervision',[sys.executable,'scripts/runMatlabResearch.py','main01',"addpath(pwd); addpath('tests'); testAerothermalEvidenceLedger('outputs/v176_validation/unit01'); testAerothermalEvidenceWorkbench('outputs/v176_validation/main01',150);",'--folder','outputs/v176_validation','--input','outputs/v104_validation/baseline02/current_case.mat','--input','outputs/v167_validation/main01/result.mat'],300)
