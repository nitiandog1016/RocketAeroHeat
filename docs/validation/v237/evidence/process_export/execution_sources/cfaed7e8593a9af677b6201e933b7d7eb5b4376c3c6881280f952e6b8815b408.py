"""模块职责：将真实重启动各阶段映射到原箭体工作台，提供固定色标和可核对的研究图。"""
from pathlib import Path
import csv,json,re,subprocess,sys
import numpy as np
import runPhaseBranchV237 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC


def main():
    audit=read(BASE/'real_audit01/report.json');assert audit['passed']
    out=BASE/'workbench01';out.mkdir(exist_ok=False);g=np.load(str(run.GEOMETRY));n=19752
    # 模块一：保持现有工作台字段接口；更新号仅是迭代标签，不注入时间历程。
    model=read(ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json')
    files=sorted((BASE/'real_audit01').glob('original_faces_*.csv'));frames={int(p.stem.split('_')[-1]):np.loadtxt(str(p),delimiter=',',skiprows=1) for p in files}
    assert 0 in frames and len(frames)>1
    fields=[];ranges={};source={p.name:sha(p) for p in files}
    method='SU2 8.6 · 相态标准Air/线性输运 · general Roe/隐式Euler · SA-neg WITHFT2 · 相态可行性步长守卫'
    scope='M5，30 km基准环境，α20°，β0°，给定Tw300 K；v222检查点199的具名重启动分支，CFL历史重置；未收敛研究，空间/外域/热验证未通过。'
    specs=[('pressure','壁面压力','Pa',2,'当前primitive壁节点到原三角面顶点均值'),
           ('heatFlux','入壁对流热流','W/m²',3,'当前primitive的k×(法向邻点温度−给定Tw)/法向投影距离；原三角面顶点均值；正值入壁'),
           ('wallFluid','壁节点流体温度','K',4,'当前流动primitive；原三角面顶点均值'),
           ('wallDeviation','壁流体温度减给定Tw','K',None,'当前壁流体温度减去给定300 K；不是固体壁温误差')]
    for key,label,unit,col,detail in specs:
        vals={step:a[:,col] if col is not None else a[:,4]-a[:,5] for step,a in frames.items()}
        lo=min(float(v.min()) for v in vals.values());hi=max(float(v.max()) for v in vals.values())
        if key=='wallDeviation':lo,hi=-max(abs(lo),abs(hi)),max(abs(lo),abs(hi))
        for step,values in vals.items():
            name='%s_%d'%(key,step);ranges[name]=[lo,hi]
            fields.append(dict(key=name,label='%s · 分支更新%d'%(label,step),unit=unit,values=values.tolist(),status=[1]*n,
                methodIndex=[0]*n,methodNames=[method+'；'+detail],scope=scope,computedStateLabel='未收敛研究数值',
                stateNote='具名重启动分支更新%d；不是物理时间，也不是原轨迹的精确续算。'%step,
                successfulUpdates=step,comparisonVariable=key,availability='available'))
    old={f['key']:f for f in model['fields']}
    for key in ('wallInput','wallSolved'):old[key]['scope']=scope;fields.append(old[key])
    latest=max(frames);model.update(fields=fields,comparisonColorRanges=ranges,defaultField='heatFlux_%d'%latest,
        resultCondition=dict(caption='M5 · α20° · 给定Tw300 K · 具名重启动更新0—%d · 同变量共同色标'%latest),
        caseNote='v237原19,752面相态守卫分支研究。变量下拉框选择阶段；“对照统一色标”比较同一物理量。给定壁温、壁流体温度和未知固体壁温分开。研究数值覆盖100%不代表全面积气动热合格。',
        researchIdentity=dict(productionQualified=False,physicalQualification=False,branch='v222-state199-native-reset-v237-phase-feasibility',
            successfulUpdates=list(frames),exactHistoricalContinuation=False,frameSemantics='steady branch updates, not physical time',
            sourceFaceCSV_SHA256=source,auditSHA256=sha(BASE/'real_audit01/report.json'),stopReason=audit['nativeCompletion']['stopReason'] if audit['nativeCompletion'] else 'native-task-failed',actualCompletedUpdates=audit['successfulUpdates'],latestAvailableSurfaceUpdate=latest))
    payload=json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    target=out/'rocket_restart_comparison.json';target.write_text(payload+'\n',encoding='utf8')
    (out/'case-data.js').write_text('window.ROCKET_CASE='+payload+';\n',encoding='utf8')
    html=(ROOT/'web/index.html').read_text(encoding='utf8')
    html=re.sub(r'((?:src|href)=")((?:viewer|section|export)\.(?:js|css)(?:\?[^\"]*)?)(")',r'\1../../../web/\2\3',html)
    (out/'index.html').write_text(html,encoding='utf8')
    restored=read(target)
    assert restored['faces']==g['originalFaces'].tolist() and restored['vertices']==g['originalVertices'].tolist()
    for f in restored['fields'][:-2]:
        a=frames[f['successfulUpdates']];col={'pressure':2,'heatFlux':3,'wallFluid':4}.get(f['comparisonVariable'])
        expected=a[:,col] if col is not None else a[:,4]-a[:,5];assert np.array_equal(f['values'],expected)
    # 模块三：调用现有前端实际纯函数，核对源值、状态与色标；不冒称浏览器验收。
    result=subprocess.run(['node',str(ROOT/'tests/testRestartBranchV237.js')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout.decode('utf8',errors='replace')
    save(out/'report.json',dict(passed=True,faces=n,fields=len(fields),updates=ordered,originalAreaM2=float(sum(g['faceArea'])),
        payloadSHA256=sha(target),dataRoundTripExact=True,
        scientificBrowserInteractionTested=False,sharedGUIChanged=False,physicalQualification=False))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.run.configure();locks=[Path(__file__).resolve(),ROOT/'scripts/runPhaseBranchV237.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            BASE/'real_audit01/report.json',run.GEOMETRY,ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json',
            ROOT/'tests/testRestartBranchV237.js']+[p for p in (ROOT/'web').iterdir() if p.name in ('index.html','viewer.js','viewer.css','section.js','export.js')]
        locks += list((BASE/'real_audit01').glob('original_faces_*.csv'))+[run.CASE/'stages.csv']
        raise SystemExit(0 if run.run.guard.supervise('export',[sys.executable,str(Path(__file__).resolve()),'task'],60,1024,locks,BASE,'audit') else 1)
