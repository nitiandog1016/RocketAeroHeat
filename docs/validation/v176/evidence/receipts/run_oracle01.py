from pathlib import Path
import sys
sys.path.insert(0,str(Path.cwd()/'scripts'))
from validateIdealAir import run
folder=Path('outputs/v176_validation')
main=folder/'main02'
command=[sys.executable,'tests/verifyAerothermalEvidenceLedger.py','--fixture',str(folder/'unit01/fixture.json'),'--engineering-json',str(main/'engineering_data.json'),'--research-json',str(main/'research_data.json'),'--csv']
command += [str(main/name) for name in ['engineering_pressure.csv','engineering_heatFlux.csv','research_pressure.csv','research_heatFlux.csv','research_wallSolved.csv']]
command += ['--output',str(folder/'oracle01')]
run(folder,'oracle01_supervision',command,90)
