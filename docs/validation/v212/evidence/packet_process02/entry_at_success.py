"""一次有预算的离线接口验证；复用已有进程监督，不构建或运行流动引擎。"""
from pathlib import Path
import sys, json
from runEilmerProcess import run
ROOT=Path(__file__).resolve().parents[1]

if __name__=='__main__':
    locks=[ROOT/p for p in ['scripts/sequentialRegionSnapshot.py','tests/testSequentialRegionSnapshot.py',
        'scripts/phaseAwareViscousAir.py','scripts/phaseAwareDryAir.py','data/thermo/phase_aware_dry_air_v1.json',
        'docs/validation/v212/PLAN.md','scripts/runSequentialSnapshotAuditV212.py']]
    work=ROOT/'outputs/v212_validation/packet_work02';work.mkdir(exist_ok=False)
    result=run([sys.executable,str(ROOT/'tests/testSequentialRegionSnapshot.py'),str(work/'audit')],work,
        ROOT/'outputs/v212_validation/packet_process02',seconds=60,output_mib=16,locks=locks,
        memory_mib=512,group_rss_mib=512)
    print(json.dumps({k:result.get(k) for k in ['passed','nativeExitCode','elapsedSeconds','peakObservedGroupRSSBytes','inputsStable']}))
