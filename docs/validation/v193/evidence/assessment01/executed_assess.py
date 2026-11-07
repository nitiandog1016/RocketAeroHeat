# -*- coding: utf-8 -*-
"""实际Euler网格资源量级与来源归档；估计不自动准入原生运行。"""
from pathlib import Path
import argparse
import json
import shutil
import time
import numpy as np
from readEilmerNative import sha
from tetgenFixedInterface import read_native
from reviewSU2ExplicitQuality import quality

ROOT=Path(__file__).resolve().parents[1]


def assess(base,out):
    start=time.monotonic();base,out=map(Path,(base,out));out.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(__file__,out/'executed_assess.py')
    mesh=json.loads((base/'mesh01/report.json').read_text())
    assert mesh['geometryPassed'] and not mesh['qualityScreenPassed']
    n=mesh['cells']+mesh['ghostCells'];assert n*8*5*8==mesh['stageArrayPayloadBytes']
    calibration=[]
    for name,c in [('coarse',2400+2720),('fine',9600+10240)]:
        path=ROOT/('outputs/v192_validation/%s_run01/native/report.json'%name)
        r=json.loads(path.read_text());assert r['passed']
        calibration.append(dict(profile=name,cellObjectsIncludingGhost=c,peakRSSBytes=r['peakObservedGroupRSSBytes'],
                                receiptSHA256=sha(path)))
    estimates=[v['peakRSSBytes']/v['cellObjectsIncludingGhost']*n for v in calibration]
    Q,T,B,ref,mark=read_native(base/'repair01/native/native_result.dat');q=quality(Q[T])
    # 固定官方源只读留档：reconstructmesh接受显式边约束；现包装没有填写该字段。
    specs=[(ROOT/'tools/runtime/tetgen/1.5.1/source/tetgen.cxx',[(23676,23760)]),
           (ROOT/'third_party/tetgen/fixed_tetra_repair.cpp',[(31,74)]),
           (ROOT/'scripts/optimizeSU2FanExterior.py',[(30,75)]),
           (ROOT/'airThermalProperties.m',[(1,58)]),
           (ROOT/'scripts/eilmerAirThermo.py',[(1,61)])]
    excerpts=[]
    for path,ranges in specs:
        lines=path.read_text(encoding='utf-8').splitlines()
        excerpts.append(dict(path=str(path),sha256=sha(path),ranges=[dict(start=a,end=min(b,len(lines)),
            text='\n'.join('%d: %s'%(i,lines[i-1]) for i in range(a,min(b,len(lines))+1))) for a,b in ranges]))
    (out/'source_excerpts.json').write_text(json.dumps(excerpts,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    result=dict(passed=True,meshCandidateAccepted=False,nativePreparationAuthorized=False,CFDAuthorized=False,
        actualCells=mesh['cells'],actualFaces=mesh['faces'],actualVertices=mesh['vertices'],actualGhosts=mesh['ghostCells'],
        exactStageArrayPayloadBytes=mesh['stageArrayPayloadBytes'],calibration=calibration,
        illustrativeRSSScalingBytes=estimates,fullPeakBytes=None,fitProven=False,
        scalingLimit='RSS per cell+ghost only; faces/stencils/distribution differ; not a bound or admission criterion',
        rejectedRefineMinimumSICN=float(q.min()),rejectedRefineCellsBelow01=int(np.sum(q<.1)),
        boundaryConstraintFinding='current refine wrapper supplies triangles but no explicit edgelist; upstream supports marked edges; not a proven quality remedy',
        sourceSHA256={str(p):sha(p) for p in [Path(__file__),base/'mesh01/report.json',base/'repair01/native/native_result.dat']},
        physicalQualification=False,newSolvedHeatAreaM2=0,elapsedSeconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['passed','illustrativeRSSScalingBytes','rejectedRefineMinimumSICN','rejectedRefineCellsBelow01','elapsedSeconds']}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();assess(a.base,a.out)
