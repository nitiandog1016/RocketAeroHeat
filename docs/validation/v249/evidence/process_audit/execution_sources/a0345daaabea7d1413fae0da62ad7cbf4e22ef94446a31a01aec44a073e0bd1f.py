"""模块职责：补充退出后的强制容量核对，避免快速末输出逃过周期采样；原监督回执不改。"""
from pathlib import Path
from runSU2RocketPhaseV222 import read,save,size_bytes


def terminal_capacity(base,process,limit_mib):
    base=Path(base).resolve();process=Path(process).resolve()
    if process.parent!=base:raise ValueError('CAPACITY_PROCESS_OUTSIDE_BASE')
    receipt=read(process/'report.json')
    if receipt['status']!='terminal' or receipt['liveGroupAfterCleanup']:raise ValueError('CAPACITY_REQUIRES_TERMINAL')
    # 直接枚举真实BASE，不复用可能被调用方改写的扫描口径；硬链接仍按各路径逻辑长度计。
    total=sum(p.stat().st_size for p in base.rglob('*') if p.is_file())
    report=dict(passed=total<=limit_mib*1024**2,logicalBytes=total,limitBytes=int(limit_mib*1024**2),
      originalSupervisorPassed=receipt['passed'],originalReportPreserved=True,measurementAfterExit=True)
    save(process/'terminal_capacity.json',report)
    return report


def supervise(guard,base,stage,argv,seconds,limit_mib,locks,cwd,kind='audit'):
    original=guard.supervise(stage,argv,seconds,limit_mib,locks,cwd,kind)
    check=terminal_capacity(base,Path(base)/('process_'+stage),limit_mib)
    return original and check['passed']
