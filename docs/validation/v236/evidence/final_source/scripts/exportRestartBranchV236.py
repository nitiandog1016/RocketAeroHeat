"""模块职责：将真实重启动各阶段映射到原箭体工作台，提供固定色标和可核对的研究图。"""
from pathlib import Path
import csv,json,re,subprocess,sys
import numpy as np
import runRestartBranchV236 as run
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC


def main():
    audit=read(BASE/'audit03/report.json');assert audit['passed']
    out=BASE/'workbench01';out.mkdir(exist_ok=False);g=np.load(str(run.GEOMETRY));n=19752
    # 模块一：保持现有工作台字段接口；更新号仅是迭代标签，不注入时间历程。
    model=read(ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json')
    files=sorted((BASE/'audit03').glob('original_faces_*.csv'));frames={int(p.stem.split('_')[-1]):np.loadtxt(str(p),delimiter=',',skiprows=1) for p in files}
    assert 0 in frames and audit['successfulUpdates'] in frames
    fields=[];ranges={};source={p.name:sha(p) for p in files}
    method='SU2 8.6 · 相态标准Air/线性输运 · general Roe/隐式Euler · SA-neg WITHFT2'
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
    latest=audit['successfulUpdates'];model.update(fields=fields,comparisonColorRanges=ranges,defaultField='heatFlux_%d'%latest,
        resultCondition=dict(caption='M5 · α20° · 给定Tw300 K · 具名重启动更新0—%d · 同变量共同色标'%latest),
        caseNote='v236原19,752面真实重启动研究。变量下拉框选择阶段；“对照统一色标”比较同一物理量。给定壁温、壁流体温度和未知固体壁温分开。研究数值覆盖100%不代表全面积气动热合格。',
        researchIdentity=dict(productionQualified=False,physicalQualification=False,branch='v222-state199-native-reset-v236',
            successfulUpdates=list(frames),exactHistoricalContinuation=False,frameSemantics='steady branch updates, not physical time',
            sourceFaceCSV_SHA256=source,auditSHA256=sha(BASE/'audit03/report.json'),stopReason=audit['nativeCompletion']['stopReason']))
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
    # 模块二：科学对照图保持原几何长宽比例，各变量初末共同色标，不修饰或补齐研究结果。
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize,SymLogNorm
    tri=g['originalVertices'][g['originalFaces']];right=np.array([1.,0.,.10]);right/=np.linalg.norm(right)
    up=np.array([-.025,1.,.25]);up-=right*np.dot(up,right);up/=np.linalg.norm(up)
    projected=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(4,2,figsize=(14,9.5),facecolor='white',gridspec_kw=dict(height_ratios=[1,1,1,1.4]))
    for row,(key,col,scale,label,cmap) in enumerate((('pressure',2,1e-3,'Wall pressure [kPa]','viridis'),
          ('heatFlux',3,1e-6,'Heat flux into wall [MW/m²]','coolwarm'),('wallDeviation',None,1,'Wall fluid T minus given Tw [K]','coolwarm'))):
        lo,hi=np.array(ranges[key+'_0'])*scale
        if key=='heatFlux':limit=max(abs(lo),abs(hi));norm=SymLogNorm(linthresh=.01,vmin=-limit,vmax=limit,base=10)
        else:norm=Normalize(lo,hi)
        for column,step in enumerate((0,latest)):
            ax=axs[row,column];a=frames[step];values=(a[:,col] if col is not None else a[:,4]-a[:,5])*scale
            poly=PolyCollection(projected[order],linewidths=0,facecolors=plt.get_cmap(cmap)(norm(values[order])))
            ax.add_collection(poly);ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.margins(.02);ax.set_axis_off()
            ax.set_title('%s | branch update %d'%(label,step),loc='left',fontsize=9)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axs[row,:].tolist(),fraction=.016,pad=.02)
    stages=[]
    with (run.CASE/'stages.csv').open() as f:
        for r in csv.DictReader(f):
            if r['stage']=='after_time':stages.append(r)
    for j,name in enumerate(('mass','x momentum','y momentum','z momentum','physical energy','SA')):
        selected=[r for r in stages if int(r['system'])==(0 if j<5 else 2)]
        axs[3,0].plot([int(r['update']) for r in selected],[np.log10(max(float(r['rms%d'%j]),1e-300)) for r in selected],label=name,lw=1)
    axs[3,0].axhline(-8,c='black',ls='--',lw=.7);axs[3,0].set_xlabel('Branch updates');axs[3,0].set_ylabel('Physical component log10 RMS');axs[3,0].legend(fontsize=7,ncol=3);axs[3,0].grid(alpha=.2)
    ordered=sorted(frames);axs[3,1].plot(ordered,[g['faceArea']@frames[s][:,3]/1e6 for s in ordered],marker='o',ms=3,color='#265b79')
    axs[3,1].set_xlabel('Branch updates');axs[3,1].set_ylabel('Original-surface heat power [MW]');axs[3,1].grid(alpha=.2)
    fig.suptitle('RocketAeroHeat | M5 / alpha 20 deg / given Tw 300 K | original 19,752 faces\nNamed restart branch: current-state wall diagnostics; NOT physically qualified',fontsize=12,y=.98)
    fig.text(.07,.015,'SU2 8.6 | phase-aware Air | SA-neg WITHFT2 | fixed comparison scales\nIteration index is not physical time. Solid wall temperature is unsolved. Signed heat flux retained.',fontsize=9,color='#444444')
    fig.subplots_adjust(top=.89,bottom=.11,hspace=.4,wspace=.25)
    fig.savefig(str(out/'rocket_restart_comparison.png'),dpi=160);plt.close(fig)
    # 模块三：调用现有前端实际纯函数，核对源值、状态与色标；不冒称浏览器验收。
    result=subprocess.run(['node',str(ROOT/'tests/testRestartBranchV236.js')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout.decode('utf8',errors='replace')
    save(out/'report.json',dict(passed=True,faces=n,fields=len(fields),updates=ordered,originalAreaM2=float(sum(g['faceArea'])),
        payloadSHA256=sha(target),imageSHA256=sha(out/'rocket_restart_comparison.png'),dataRoundTripExact=True,
        scientificBrowserInteractionTested=False,sharedGUIChanged=False,physicalQualification=False))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        run.configure_guard();locks=[Path(__file__).resolve(),ROOT/'scripts/runRestartBranchV236.py',ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
            BASE/'audit03/report.json',run.GEOMETRY,ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json',
            ROOT/'tests/testRestartBranchV236.js']+[p for p in (ROOT/'web').iterdir() if p.name in ('index.html','viewer.js','viewer.css','section.js','export.js')]
        locks += list((BASE/'audit03').glob('original_faces_*.csv'))+[run.CASE/'stages.csv']
        raise SystemExit(0 if run.guard.supervise('export',[sys.executable,str(Path(__file__).resolve()),'task'],60,1024,locks,BASE,'audit') else 1)
