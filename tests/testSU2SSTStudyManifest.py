"""核对未来准备器更新任务控制字段，不改写本轮已运行的原始输入。"""
from pathlib import Path
import sys,json,copy
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from prepareSU2CheckpointContinuation import read,sha
from prepareSU2SSTNumericalStudy import current_manifest

def main():
    folder=ROOT/'outputs/v127_validation/control';source=read(folder/'case_manifest.json');before=copy.deepcopy(source)
    corrected=current_manifest(source,source['sst'],source['initialization'],source['sstContractSHA256'])
    assert source==before and source['flight']['iterations']==50 and source['flight']['outputEvery']==25
    assert corrected['flight']['iterations']==150 and corrected['flight']['outputEvery']==75
    assert corrected['sourceRunControls']==dict(iterations=50,outputEvery=25,namedCheckpoints=True)
    for k in ['Mach','alpha','beta','temperature','pressure','wallTemperature']:
        assert corrected['flight'][k]==before['flight'][k]
    assert corrected['sst']==source['sst'] and corrected['initialization']==source['initialization']
    target=ROOT/'outputs/v127_validation/manifest_control_check.json';assert not target.exists()
    target.write_text(json.dumps(dict(passed=True,physicalInputsPreserved=True,originalNativeInputsModified=False,
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),ROOT/'scripts/prepareSU2SSTNumericalStudy.py',folder/'case_manifest.json']}),indent=2)+'\n',encoding='utf-8')
    print('Manifest control metadata corrected; original inputs unchanged.')

if __name__=='__main__':main()
