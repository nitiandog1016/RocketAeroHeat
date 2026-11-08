"""模块职责：把已审核原面数据交给现有科研工作台，不改共享GUI、不授物理资格。"""
from pathlib import Path
import json,sys
import numpy as np
from runSU2RocketPhaseV222 import ROOT,BASE,CASE,DOC,read,save,sha,supervise

def main():
    audit=read(BASE/'audit01/report.json');assert audit['passed']
    last=audit['records'][-1];step=last['iteration'];out=BASE/'workbench01';out.mkdir(exist_ok=False)
    source=BASE/'audit01'/('original_faces_%06d.csv'%step);a=np.loadtxt(str(source),delimiter=',',skiprows=1)
    g=np.load(str(BASE/'geometry01/cache.npz'));n=len(a)
    method='SU2 8.6 · general Roe/implicit Euler · SA-neg WITHFT2 · 相态标准Air/线性输运'
    scope='M5，p∞1197.0279007383137 Pa，T∞226.50908361133003 K，攻角20°，给定Tw300 K；60—2000 K及相态联合域。未收敛，空间/外域误差及公开热验证未完成。'
    fields=[]
    for key,label,unit,column,status,algorithm in (
        ('pressure','壁面压力','Pa',2,1,method+'；原生壁节点到原三角面顶点均值'),
        ('heatFlux','入壁对流热流','W/m²',3,1,method+'；原生等温壁法向温度梯度；正值入壁'),
        ('wallFluid','壁节点流体温度','K',4,1,method+'；流动更新后原始量，面顶点均值'),
        ('wallInput','给定壁温','K',5,2,'等温壁输入边界条件'),
        ('wallSolved','独立固体壁温','K',None,0,'未求解：未运行固体导热或独立热平衡')):
        fields.append(dict(key=key,label=label,unit=unit,values=[None]*n if column is None else a[:,column].tolist(),
            status=[status]*n,methodIndex=[0]*n,methodNames=[algorithm],scope=scope,
            computedStateLabel='未收敛研究数值',stateNote='研究场，不计正式合格求解面积。' if status==1 else
                '给定边界，不能当作求解壁温。' if status==2 else '未执行固体壁温求解；不以给定300K补齐。',
            availability='notProvidedInSnapshot' if status==0 else 'available'))
    model=dict(kind='aerothermal-model',schemaVersion=1,vertices=g['originalVertices'].tolist(),faces=g['originalFaces'].tolist(),
        faceArea=a[:,1].tolist(),totalArea=float(sum(a[:,1])),region=['完整原壁面']*n,
        regions=[dict(name='完整原壁面',area=float(sum(a[:,1])),faces=n)],missingReason=['']*n,heatingGaps=[],fields=fields,
        defaultField='heatFlux',coordinateNote='SI · X主轴 · 鼻尖+X · 真实几何比例',
        flowDirection=[-float(np.cos(np.deg2rad(20))),0,float(np.sin(np.deg2rad(20)))],
        resultCondition=dict(caption='M5 · 30 km基准环境 · α20° · β0° · 给定Tw300 K · 更新%d · 未收敛研究场'%(step+1)),
        caseNote='v222真实完整箭体。研究数值面积100%不等于全面积气动热已验证。独立固体壁温未求解；负热流保留。',
        researchIdentity=dict(productionQualified=False,physicalQualification=False,model='rah-su2-phase-air-saneg-eplus300k-v1',
            iteration=step,successfulUpdates=step+1,originalFaceCSV_SHA256=sha(source),auditSHA256=sha(BASE/'audit01/report.json')))
    target=out/'rocket_phase_research.json';target.write_text(json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False)+'\n',encoding='utf8')
    # 导出后逐数回读，确认压力/有符号热流/温度没有插值或补零改变。
    restored=read(target)
    assert restored['faces']==g['originalFaces'].tolist() and restored['vertices']==g['originalVertices'].tolist()
    for field,column in zip(restored['fields'][:4],(2,3,4,5)):assert np.array_equal(field['values'],a[:,column])
    assert all(v is None for v in restored['fields'][-1]['values'])
    save(out/'report.json',dict(passed=True,originalFaces=n,areaM2=float(sum(a[:,1])),exportSHA256=sha(target),
        sourceSHA256=sha(source),physicalQualification=False,solidWallTemperatureSolved=False,
        savedNegativeHeatFaces=int(np.sum(a[:,3]<0)),browserInteractionTested=False,sharedGUIChanged=False))

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        locks=[Path(__file__).resolve(),BASE/'audit01/report.json',BASE/'geometry01/cache.npz',DOC/'PLAN.md']
        raise SystemExit(0 if supervise('export',[sys.executable,str(Path(__file__).resolve()),'task'],30,256,locks,BASE,'audit') else 1)
