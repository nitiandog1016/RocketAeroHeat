"""显式原生完成审核：失败文本、时间身份、保存场和载荷分别核对。"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
from eilmerExplicitProfile import audit_explicit_config,DEVELOPMENT_PROFILE


def require(ok, message):
    if not ok:
        raise ValueError(message)


def audit_explicit_completion(case, stdout, manifest):
    # 模块一：实际二进制和停止信息；不能由包装器或原生exit0推断成功。
    case=Path(case);text=Path(stdout).read_text(errors='replace');condition=manifest['condition']
    cfg=json.loads((case/'lmrsim/config').read_text());control=json.loads((case/'lmrsim/control').read_text())
    audit_explicit_config(cfg,condition,control)
    require(not re.search(r'OutOfMemoryError|Memory allocation failed|core\.exception\.|integrateInTime failed|An exception was caught|Begin exception message|Halt set in control file',text), 'native-explicit-failure')
    require(re.findall(r'^Number-type:\s*(\S+)\s*$',text,re.M)==['real'], 'wrong-explicit-number-type')
    reasons=re.findall(r'^STOP-REASON:\s*(.+)$',text,re.M)
    steps=re.findall(r'^FINAL-STEP:\s*(\d+)\s*$',text,re.M)
    times=re.findall(r'^FINAL-TIME:\s*(\S+)\s*$',text,re.M)
    require(len(reasons)==len(steps)==len(times)==1,'missing-explicit-stop-markers')
    step=int(steps[0]);pseudo=float(times[0]);d=condition['explicit']
    development=d['profile']==DEVELOPMENT_PROFILE
    if development:
        from eilmerExplicitAdmission import audit_admission
        audit_admission(manifest,case)
        require(0<step<=d['maxSteps'] and
                ((reasons==['maximum-steps'] and step==d['maxSteps']) or
                 (reasons==['maximum-wall-clock'] and step<d['maxSteps'])), 'incomplete-explicit-step-window')
    else:require(reasons==['maximum-steps'] and step==d['maxSteps'],'incomplete-explicit-step-window')
    # 直接lmr-run默认不打印runsim的verbosity提示；用实际时间步记录和保存元数据核对分派。
    observed_steps=[int(v) for v in re.findall(r'^Step=\s*(\d+)\s+t=',text,re.M)]
    require(observed_steps==list(range(1,step+1)),'wrong-explicit-dispatch')
    require(np.isfinite(pseudo) and 0<pseudo<1.0,'invalid-explicit-pseudo-time')
    # 模块二：逐步残差和准确保存步骤；dt仅表示局部伪时间，不授物理时间身份。
    history=np.atleast_2d(np.loadtxt(case/'lmrsim/transient-residuals.dat'))
    require(history.shape==(step,10) and np.isfinite(history).all(),'invalid-explicit-history')
    require(np.array_equal(history[:,0],np.arange(1,step+1)),'explicit-history-step-mismatch')
    require(np.all(np.diff(history[:,1])>0),'nonincreasing-pseudo-time')
    snapshots=(case/'lmrsim/snapshots/snapshot-times-metadata').read_text()
    entries=re.findall(r"^'(\d+)':\s*\n\s+step:\s*(\d+)\s*\n\s+time:\s*(\S+)",snapshots,re.M)
    expected_steps=[0]+([d['checkpointStep']] if step>=d['checkpointStep'] else [])+[step]
    require([int(e[1]) for e in entries]==expected_steps, 'missing-or-duplicate-explicit-snapshot')
    require([int(e[0]) for e in entries]==list(range(len(entries))),'explicit-snapshot-identity-mismatch')
    if step==d['checkpointStep']:
        # 原生手动检查点未置output_just_written；恰好停在检查点时终帧可能再次写入。
        paths=[case/'lmrsim/snapshots'/('%04d'%i)/'fluid-0000' for i in [1,2]]
        require(hashlib.sha256(paths[0].read_bytes()).digest()==hashlib.sha256(paths[1].read_bytes()).digest(),'duplicate-final-snapshot-differs')
    # 源码打印FINAL-TIME为%g，不能要求文本六位与原始double逐字相等。
    saved_time=float(entries[-1][2])
    require(abs(saved_time-pseudo)<=5.1e-6*saved_time,'explicit-time-mismatch')
    from readEilmerNative import read_field
    dt_ranges=[]
    for idx, saved_step, _ in entries:
        field=read_field(case,int(idx),manifest['mesh']['cells'])
        require(all(np.all(field[k]>0) for k in ['rho','p','T','vol']),'invalid-explicit-physical-state')
        if int(saved_step)>0:
            dt=field['dt_local'];require(np.all(dt>0) and np.all(dt<=1e-3),'invalid-explicit-local-dt')
            dt_ranges.append([float(dt.min()),float(dt.max())])
    # 模块三：检查每个保存的热壁载荷，完整原面覆盖和正负热流均保留。
    meta=np.atleast_2d(np.loadtxt(case/'lmrsim/loads/loads-metadata',skiprows=1))
    # 原生initLoadsFiles只建立目录/表头；第一次实际载荷为检查点，不虚构初始壁载荷。
    expected_loads=([d['checkpointStep']] if step>=d['checkpointStep'] else [])+([step] if step!=d['checkpointStep'] else [])
    require(meta.shape==(len(expected_loads),3) and np.array_equal(meta[:,2],expected_loads),'explicit-load-steps-mismatch')
    require(np.array_equal(meta[:,0],np.arange(len(expected_loads))),'duplicate-explicit-load-index')
    for row in meta:
        count=0
        for path in (case/'lmrsim/loads'/('%04d'%int(row[0]))).glob('*-wall.dat'):
            a=np.atleast_1d(np.genfromtxt(path,names=True,deletechars=''))
            require(all(np.isfinite(a[k]).all() for k in a.dtype.names),'nonfinite-explicit-loads')
            count+=len(a)
        require(count==manifest['mesh']['wallFaces'],'incomplete-explicit-wall-coverage')
    require(abs(meta[-1,1]-saved_time)<=5.1e-6*saved_time,'explicit-load-time-mismatch')
    return dict(passed=True,stopReason=reasons[0],finalStep=step,wallFaces=count,
                savedSteps=expected_steps,localDtRanges=dt_ranges,requestedSteps=d['maxSteps'],completedRequestedSteps=step==d['maxSteps'],
                finalPseudoTime=saved_time,timeMeaning='local-pseudo-time-not-flight-time',
                physicsQualified=False,newQualifiedAreaM2=0,
                nativeLogSHA256=hashlib.sha256(Path(stdout).read_bytes()).hexdigest(),
                scope='Explicit interface completion only; not steady convergence, gradient replay or physical qualification')
