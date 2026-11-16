"""模块职责：导出完整原箭体上的前体研究字段、区域外未知值及可追溯方法，不重算流动。"""
from pathlib import Path
import json,re,subprocess,sys
import numpy as np
import runForebodySpatialV243 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC
CLOSED=True


def export():
    out=BASE/'workbench01';out.mkdir(exist_ok=False)
    assert read(BASE/'process_flow_audit2/report.json')['passed'];a=read(BASE/'flow_audit02/report.json')
    model=read(ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json');n=len(model['faces']);assert n==19752
    candidates=sorted((BASE/'flow_audit02').glob('original_faces_*.csv'))
    selected={0,len(candidates)//2,len(candidates)-1};files=[p for i,p in enumerate(candidates) if i in selected]
    frames={int(p.stem.split('_')[-1]):np.loadtxt(str(p),delimiter=',',skiprows=1) for p in files}
    ids=frames[0][:,0].astype(int)-1;mask=np.zeros(n,bool);mask[ids]=True
    region_state=mask.astype(int).tolist();method_index=(~mask).astype(int).tolist()
    absent=['' if v else '本次前体区域之外，未求解' for v in mask]
    method='SU2 8.6 · 同相态Air Euler · Roe/二阶MUSCL · 三维前体平面出口 · 原生全出射特征'
    scope='M5，α20°；细化前体Euler，17736计算子面按实际面积聚合到4434原父面。残差/净收支未合格；粗细分支均无真实边界层外缘或壁热资格。'
    fields=[];ranges={}
    specs=[('pressure','壁面压力','Pa',2),('fluid','Euler壁邻流体温度','K',3),('H','物理总焓','J/kg',4),
      ('S','物理熵','J/(kg·K)',5),('Mach','局部Mach数','1',9),('Hdeviation','总焓相对来流偏差','%',None)]
    for key,label,unit,col in specs:
        values={s:f[:,col] if col is not None else 100*(f[:,4]/a['freestreamHJkg']-1) for s,f in frames.items()}
        lo=min(float(v.min()) for v in values.values());hi=max(float(v.max()) for v in values.values())
        if key=='Hdeviation':lo,hi=-max(abs(lo),abs(hi)),max(abs(lo),abs(hi))
        for step,v in values.items():
            full=[None]*n
            for i,value in zip(ids,v):full[int(i)]=float(value)
            name='%s_%d'%(key,step);ranges[name]=[lo,hi]
            fields.append(dict(key=name,label='%s · 更新%d'%(label,step),unit=unit,values=full,status=region_state,
              methodIndex=method_index,methodNames=[method,'未求解：本次区域之外'],missingReason=absent,
              scope=scope,computedStateLabel='外流研究值（未获外缘资格）',stateNote='区域分支更新%d；不是物理时间；%s。'%(step,'均匀初场' if step==0 else '完整原生推进'),
              successfulUpdates=step,comparisonVariable=key,availability='regionalResearchOnly'))
    for key,label,unit,value,status,note in [('heatFlux','入壁对流热流','W/m²',None,0,'未求解：Euler不提供黏性壁热'),
      ('wallInput','给定壁温','K',300,2,'给定300 K；本Euler计算未施加热壁'),('wallSolved','求解固体壁温','K',None,0,'未求解：没有执行固体导热')]:
        fields.append(dict(key=key,label=label,unit=unit,values=[value]*n,status=[status]*n,methodIndex=[0]*n,
          methodNames=[note],scope=scope,stateNote=note,availability='notProvidedInSnapshot' if status==0 else 'given'))
    latest=max(frames);model.pop('heatingGaps',None)
    model.update(fields=fields,comparisonColorRanges=ranges,defaultField='pressure_%d'%latest,
      missingReason=['Euler不提供黏性壁热；工程边界层入口尚未合格']*n,
      resultCondition=dict(caption='M5 · α20° · 前体区域更新0—%d · 全箭原面保留'%latest),
      caseNote='v243细化区域：完成%d更新；原父面压力面积变化%.6g%%。未取得物理资格；粗细分支差异含不同迭代误差，不是网格无关证明。'%(a['completedUpdates'],100*a['finalPressureAreaChange']),
      researchIdentity=dict(productionQualified=False,physicalQualification=False,branch='phase-air-refined-forebody-euler-v243',
        successfulUpdates=sorted(frames),originalFaceIDs=(ids+1).tolist(),regionAreaM2=a['regionAreaM2'],
        frameSemantics='steady branch updates, not physical time',actualCompletedUpdates=latest,stopReason=a['stopReason'],
        sourceFaceCSV_SHA256={p.name:sha(p) for p in files},auditSHA256=sha(BASE/'flow_audit02/report.json'),
        originalPressureTrendWrong=False,displayUsesCorrectGlobalMapping=True,newQualifiedAreaM2=0))
    payload=json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    target=out/'forebody_research.json';target.write_text(payload+'\n',encoding='utf8')
    (out/'case-data.js').write_text('window.ROCKET_CASE='+payload+';\n',encoding='utf8')
    html=(ROOT/'web/index.html').read_text(encoding='utf8')
    html=re.sub(r'((?:src|href)=")((?:viewer|section|export)\.(?:js|css)(?:\?[^\"]*)?)(")',r'\1../../../web/\2\3',html)
    (out/'index.html').write_text(html,encoding='utf8')
    restored=read(target);g=np.load(str(ROOT/'outputs/v230_validation/audit02/geometry.npz'))
    assert restored['vertices']==g['vertices'].tolist() and restored['faces']==g['faces'].tolist()
    count=0
    for f in restored['fields'][:-3]:
        step=f['successfulUpdates'];col=dict((k,c) for k,_,_,c in specs)[f['comparisonVariable']]
        expected=frames[step][:,col] if col is not None else 100*(frames[step][:,4]/a['freestreamHJkg']-1)
        assert np.array_equal(np.array(f['values'],object)[ids].astype(float),expected)
        assert all(v is None for i,v in enumerate(f['values']) if not mask[i]);count+=len(ids)
    test=subprocess.run(['node',str(ROOT/'tests/testForebodyWorkbenchV243.js')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(test.stdout);assert test.returncode==0,test.stdout.decode('utf8',errors='replace')
    save(out/'report.json',dict(passed=True,fields=len(fields),updates=sorted(frames),exactResearchValues=count,
      regionFaces=len(ids),originalFaces=n,regionAreaM2=float(np.sum(g['area'][ids])),originalAreaM2=float(np.sum(g['area'])),
      unknownOutsideRegionPreserved=True,actualFrontendPureFunctionsPassed=True,browserInteractionTested=False,sharedWebChanged=False,
      physicalQualification=False,payloadSHA256=sha(target)))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v243 spatial export closed; use NEXT_DECISION')
    if len(sys.argv)>1 and sys.argv[1]=='task':export()
    else:
        run.guard.BASE=BASE;run.guard.CASE=BASE;run.guard.DOC=DOC
        size=run.guard.size_bytes;run.guard.size_bytes=lambda p:size(BASE) if Path(p).name.startswith('process_') else 0
        locks=[Path(__file__).resolve(),ROOT/'scripts/runForebodySpatialV243.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
          BASE/'flow_audit02/report.json',BASE/'process_flow_audit2/report.json',ROOT/'tests/testForebodyWorkbenchV243.js',
          ROOT/'outputs/v230_validation/audit02/geometry.npz',ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json']
        locks += list((BASE/'flow_audit02').glob('original_faces_*.csv'))+[ROOT/'web'/p for p in ('index.html','viewer.js','viewer.css','section.js','export.js')]
        raise SystemExit(0 if run.guard.supervise('export',[sys.executable,str(Path(__file__).resolve()),'task'],60,384,locks,BASE,'audit') else 1)
