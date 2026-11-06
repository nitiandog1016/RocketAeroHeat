"""模块用途。复用已有主界面，只验收新限步结果/方法及旧结果兼容，独立复算导出。"""
from pathlib import Path
import sys
from su2ThermalRuntime import ROOT,read
from validateIdealAir import run


def main():
    base=ROOT/'outputs/v165_validation'
    assert read(base/'review01/report.json')['passed']
    assert not (base/'main01').exists()
    command=[sys.executable,str(ROOT/'scripts/runMatlabResearch.py'),'main01',
        "addpath(pwd); addpath(fullfile(pwd,'tests')); testSU2ConservativeStepWorkbench(fullfile(pwd,'outputs','v165_validation','main01'));",
        '--folder','outputs/v165_validation/matlab']
    for name in ('outputs/v118_validation/actual01/surface.mat',
                 'outputs/v56_validation/state_validity/state_validity_actual_rocket.mat',
                 'outputs/v145_validation/mesh/hybrid_volume.mat',
                 'outputs/v145_validation/mesh/mesh_input.mat',
                 'outputs/v145_validation/mesh/mesh_report.json',
                 'outputs/v145_validation/mesh/mesh_process.json',
                 'outputs/v165_validation/real01/case_manifest.json',
                 'outputs/v165_validation/real01/engine_process.json',
                 'outputs/v165_validation/real01/audit49/report.json'):
        # 不存在的输入在创建MATLAB进程前拒绝，不猜测替代路径。
        assert (ROOT/name).is_file(),name
        command+=['--input',name]
    run(base,'main01_parent',command,300)
    result=read(base/'matlab/main01_process.json')
    assert result['exitCode']==0 and all(result['sourceStableAfter'].values())
    sys.path.insert(0,str(ROOT/'tests'))
    from verifySU2IdealAirFullDomain import main as verify
    verify(base,'main01',49,150,('OTH','共同守恒方向有限半步','按当前k重算','限步不证明全局守恒'))


if __name__=='__main__':main()
