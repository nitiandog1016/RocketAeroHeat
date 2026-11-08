"""模块职责：只读核对准确来源、原资产和并发文件，归档小证据。"""
from pathlib import Path
import ast,hashlib,json,shutil,time

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v219_validation'
DOC=ROOT/'docs/validation/v219';EVIDENCE=DOC/'evidence'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as stream:
        for b in iter(lambda:stream.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf8'))
def save(p,value):p.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf8')

def main():
    start=time.monotonic();out=BASE/'final03';out.mkdir(exist_ok=False);EVIDENCE.mkdir(exist_ok=True)
    sources=[]
    for part in ('source01','source02'):
        report=read(BASE/part/'report.json');assert report['passed']
        for row in report['sources']:
            assert sha(ROOT/row['path'])==row['sha256']
            sources.append(row)
    assert len(sources)==29 and len({r['relative'] for r in sources})==29
    protected=read(ROOT/'outputs/v176_validation/protected_initial.json')
    for name,digest in protected.items():assert sha(ROOT/name)==digest,name
    assets={
      'outputs/v144_validation/actual/hybrid_volume.mat':'22b35b72396c1a7d14195af9bd24d3c6328138484d4009638c1ab48cddec1d6c',
      'outputs/v167_validation/real01/mesh.su2':'00f96f3a0cb87a11c4218796cd1f5f25aebffb453c794e0cf97f7bd4abde3181',
      'data/thermo/phase_aware_dry_air_v1.json':'7fdde65001d0d0aefe6b40136f269c5f512b8a94e327748854987cc831fe3985'}
    for name,digest in assets.items():assert sha(ROOT/name)==digest,name
    scripts=['inspectSU2RouteV219.py','inspectSU2RouteSupplementV219.py','runSU2RouteV219.py',
             'su2EnergyReference.py','verifySU2EnergyReferenceV219.py','finalizeSU2RouteV219.py']
    for name in scripts:ast.parse((ROOT/'scripts'/name).read_text(encoding='utf8'))
    energy=read(BASE/'energy01/report.json');assert energy['passed'] and energy['rejectionCount']==20
    processes=[]
    for stage in ('source','supplement','energy'):
        folder=BASE/('process_'+stage);r=read(folder/'report.json')
        assert r['passed'] and r['inputsStable'] and not r['liveGroupAfterCleanup']
        assert sha(folder/'executed_task.py')==r['inputSHA256'][r['argv'][-1]]
        # executed_entry是被导入的v216监督模块，实际v219调用入口另存。
        assert sha(folder/'executed_entry.py')==r['inputSHA256'][str(ROOT/'scripts/runCompleteRegionDecisionV216.py')]
        caller=BASE/('entry_'+stage)/'executed_caller.py'
        assert sha(caller)==r['inputSHA256'][str(ROOT/'scripts/runSU2RouteV219.py')]
        # 只读归档恢复：已有副本必须逐字相同，再允许复用目录，不删除首失败。
        for target,origin in [(EVIDENCE/folder.name,folder),(EVIDENCE/caller.parent.name,caller.parent)]:
            if target.exists():
                for p in target.rglob('*'):
                    if p.is_file():assert sha(p)==sha(origin/p.relative_to(target))
            else:shutil.copytree(origin,target)
        processes.append(dict(stage=stage,elapsedSeconds=r['elapsedSeconds'],peakRSS=r['peakObservedGroupRSSBytes']))
    for part in ('source01','source02','energy01'):
        (EVIDENCE/part).mkdir(exist_ok=True);shutil.copyfile(BASE/part/'report.json',EVIDENCE/part/'report.json')
    # 首计划精确字节恢复只接受原执行SHA相等，不把后续修订追称原计划。
    initial=(DOC/'PLAN.md').read_text(encoding='utf8').split('## 执行中补充')[0].rstrip()+'\n'
    old_hash=read(BASE/'process_source/report.json')['inputSHA256'][str(DOC/'PLAN.md')]
    possibilities=[initial.encode('utf8'),initial.replace('\n','\r\n').encode('utf8')]
    original=next((b for b in possibilities if hashlib.sha256(b).hexdigest()==old_hash),None)
    assert original is not None,'original PLAN bytes not recovered'
    (EVIDENCE/'initial_PLAN.md').write_bytes(original)
    save(EVIDENCE/'source_bindings.json',dict(sourceCommit=read(BASE/'source01/report.json')['sourceCommit'],sources=sources))
    # 公开只保留必要准确行段；完整源/二进制继续留E盘，固定上游与整文件SHA可复核。
    excerpts={
      'Common/src/CConfig.cpp':[(3790,3814),(9758,9771)],
      'Common/include/option_structure.hpp':[(1330,1370)],
      'SU2_CFD/include/fluid/CSU2IdealAir.hpp':[(1,63)],
      'SU2_CFD/src/solvers/CEulerSolver.cpp':[(2555,2604),(5285,5340),(5430,5490),(5491,5585)],
      'SU2_CFD/src/solvers/CNSSolver.cpp':[(630,790)],
      'SU2_CFD/src/variables/CNSVariable.cpp':[(144,230)],
      'SU2_CFD/src/solvers/CTurbSASolver.cpp':[(264,294),(530,567),(604,657)],
      'SU2_CFD/include/solvers/CScalarSolver.inl':[(370,391)],
      'SU2_CFD/src/numerics/flow/convection/roe.cpp':[(757,834),(985,1020)],
      'SU2_CFD/src/drivers/CDriver.cpp':[(1430,1450),(1665,1687)],
      'SU2_CFD/include/numerics/scalar/scalar_edge_flux.hpp':[(410,468)],
      'SU2_CFD/include/numerics/turbulent/turb_sa_edge_flux.hpp':[(34,58),(113,170)]}
    snippet_rows=[]
    for name,ranges in excerpts.items():
        row=next(x for x in sources if x['relative']==name);lines=(ROOT/row['path']).read_text(encoding='utf8').splitlines()
        for first,last in ranges:
            assert first>=1 and last<=len(lines),(name,len(lines))
            snippet_rows.append(dict(source=name,sourceSHA256=row['sha256'],firstLine=first,lastLine=last,
                                     text='\n'.join(lines[first-1:last])+'\n'))
    save(EVIDENCE/'source_excerpts.json',snippet_rows)
    for name in ('COPYING','AUTHORS.md'):
        row=next(x for x in sources if x['relative']==name);shutil.copyfile(ROOT/row['path'],EVIDENCE/name)
    for stage in ('final','final2'):
        failed=read(BASE/('process_'+stage)/'report.json');assert not failed['passed'] and failed['nativeExitCode']==1
        assert failed['inputsStable'] and not failed['liveGroupAfterCleanup']
        shutil.copytree(BASE/('process_'+stage),EVIDENCE/('process_'+stage))
        shutil.copytree(BASE/('entry_'+stage),EVIDENCE/('entry_'+stage))
    size=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
    assert size<64*1024**2
    report=dict(passed=True,protectedCount=len(protected),assets=assets,sourceCount=len(sources),
        scriptsSyntaxChecked=scripts,priorProcesses=processes,sourceOutputGuardMatchesOriginalPlan=False,
        originalOnlineOutputLimitMiB=128,originalPlannedOutputMiB=64,observedOutputBytes=size,
        originalPlanSHA256=old_hash,originalPlanBytesRecovered=True,designAccepted=True,
        buildAdmission=False,nativeAdmission=False,CFDAdmission=False,physicalQualification=False,
        firstArchiveFailure='CSU2IdealAir excerpt requested line 77 but exact file has 63 lines',
        secondArchiveFailure='Python 3.7 shutil.copytree has no dirs_exist_ok keyword',
        recovery='archive-only; no source/native/energy test repeated',
        newQualifiedAreaM2=0,seconds=time.monotonic()-start)
    save(out/'report.json',report);shutil.copyfile(out/'report.json',EVIDENCE/'final_report.json')
    print(json.dumps({k:report[k] for k in ('passed','protectedCount','sourceCount','observedOutputBytes','seconds')}))

if __name__=='__main__':main()
