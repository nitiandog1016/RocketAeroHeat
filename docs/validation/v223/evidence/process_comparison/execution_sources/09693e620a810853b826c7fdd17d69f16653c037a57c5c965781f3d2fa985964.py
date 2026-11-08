"""模块职责：把四份真实原面研究场组成统一色标比较，显式区分迭代、给定壁温和未知固体温度。"""
from pathlib import Path
import json,sys,subprocess
import numpy as np
import runSU2RocketPhaseV222 as guard
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v223_validation';DOC=ROOT/'docs/validation/v223';OLD=ROOT/'outputs/v222_validation'
STEPS=(50,100,150,199)

def main():
    out=BASE/'workbench01';out.mkdir(exist_ok=False)
    assert read(BASE/'audit02/report.json')['passed']
    model=read(OLD/'workbench01/rocket_phase_research.json')
    g=np.load(str(OLD/'geometry01/cache.npz'));n=len(model['faces'])
    frames=[];sources={};fields=[];ranges={}
    template={f['key']:f for f in model['fields']}
    specs=(('pressure','壁面压力','Pa',2),('heatFlux','入壁对流热流','W/m²',3),
           ('wallDeviation','壁节点流体温度 − 给定壁温','K',None))
    for step in STEPS:
        p=OLD/'audit01'/('original_faces_%06d.csv'%step)
        a=np.loadtxt(str(p),delimiter=',',skiprows=1)
        assert np.array_equal(a[:,0],np.arange(1,n+1)) and np.array_equal(a[:,1],model['faceArea'])
        assert np.all(a[:,5]==300) and np.isfinite(a).all()
        frames.append(a);sources[p.name]=sha(p)
    for key,label,unit,column in specs:
        values=[a[:,column] if column is not None else a[:,4]-a[:,5] for a in frames]
        lo=min(float(v.min()) for v in values);hi=max(float(v.max()) for v in values)
        if key=='wallDeviation':lo,hi=-max(abs(lo),abs(hi)),max(abs(lo),abs(hi))
        assert lo<hi
        for step,v in zip(STEPS,values):
            field=dict(template['wallFluid' if column is None else key]);name='%s_%d'%(key,step+1)
            field.update(key=name,label='%s · 第%d次更新'%(label,step+1),unit=unit,values=v.tolist(),
                successfulUpdates=step+1,nativeIteration=step,comparisonVariable=key,
                stateNote='第%d次成功更新的未收敛研究场；稳态迭代编号不是物理时间。'%(step+1))
            if column is None:
                field['methodNames']=[template['wallFluid']['methodNames'][0]+'；减去给定Tw300 K，非固体壁温误差']
            fields.append(field);ranges[name]=[lo,hi]
    model.update(fields=fields+[template['wallInput'],template['wallSolved']],comparisonColorRanges=ranges,
        defaultField='wallDeviation_200',
        resultCondition=dict(caption='M5 · α20° · 给定Tw300 K · 第51/101/151/200次更新对比 · 同变量共同色标'),
        caseNote='v223真实原面四阶段比较。变量下拉框选择阶段；色标选“对照共同范围”。温差是流体壁节点温度减给定边界温度。研究数值100%不代表物理资格；无物理时间序列。')
    model['researchIdentity']=dict(productionQualified=False,physicalQualification=False,
        model='rah-su2-phase-air-saneg-eplus300k-v1',successfulUpdates=[s+1 for s in STEPS],
        frameSemantics='steady solver updates; not physical time',sourceFaceCSV_SHA256=sources,
        auditSHA256=sha(BASE/'audit02/report.json'),frozenOperatorFlowUpdates=0)
    target=out/'rocket_phase_comparison.json'
    target.write_text(json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False)+'\n',encoding='utf8')
    restored=read(target)
    assert restored['faces']==g['originalFaces'].tolist() and restored['vertices']==g['originalVertices'].tolist()
    for f in restored['fields'][:12]:
        a=frames[STEPS.index(f['nativeIteration'])]
        column={'pressure':2,'heatFlux':3}.get(f['comparisonVariable'])
        expected=a[:,column] if column is not None else a[:,4]-a[:,5]
        assert np.array_equal(expected,f['values'])
    # 模块二：同一真实几何、同一投影和同一色标，比较壁温差变化，不改变数据。
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize
    tri=g['originalVertices'][g['originalFaces']]
    right=np.array([1.,0.,.10]);right/=np.linalg.norm(right)
    up=np.array([-.025,1.,.25]);up-=right*np.dot(up,right);up/=np.linalg.norm(up)
    projected=tri@np.column_stack([right,up]);order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(4,2,figsize=(14,8.5),facecolor='white')
    for row,(step,a) in enumerate(zip(STEPS,frames)):
        for col,(key,cmap,scale) in enumerate((('heatFlux','viridis',1e-6),('wallDeviation','coolwarm',1))):
            ax=axs[row,col];vals=(a[:,3] if col==0 else a[:,4]-a[:,5])*scale
            lo,hi=np.array(ranges['%s_%d'%(key,step+1)])*scale;norm=Normalize(lo,hi)
            poly=PolyCollection(projected[order],linewidths=0,facecolors=plt.get_cmap(cmap)(norm(vals[order])))
            ax.add_collection(poly);ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.margins(.02);ax.set_axis_off()
            ax.set_title('Update %d'%(step+1),loc='left',fontsize=10,pad=0)
    for col,(key,cmap,scale,label) in enumerate((('heatFlux','viridis',1e-6,'Heat flux into wall [MW/m²]'),
                                               ('wallDeviation','coolwarm',1,'Wall fluid T - given Tw [K]'))):
        lo,hi=np.array(ranges[key+'_200'])*scale
        cax=fig.add_axes([.09+col*.48,.11,.36,.018])
        fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(lo,hi),cmap=cmap),cax=cax,orientation='horizontal',label=label)
    fig.suptitle('RocketAeroHeat | original 19,752 faces | M5 / alpha 20 deg / given Tw 300 K\nUnconverged research: fixed ranges across solver updates, not physical time',fontsize=12,y=.97)
    fig.text(.07,.02,'SU2 8.6 | phase-aware standard Air | SA-neg WITHFT2 | geometry aspect ratio preserved\nWall-fluid temperature deviation is not a solved solid temperature. No physical qualification.',fontsize=9,color='#444444')
    fig.subplots_adjust(top=.86,bottom=.19,left=.06,right=.97,hspace=.18,wspace=.12)
    fig.savefig(str(out/'rocket_phase_comparison.png'),dpi=160);plt.close(fig)
    # 模块三：只调用现有前端纯合同，不打开浏览器、不改并发UI源。
    test=ROOT/'tests/testSU2PhaseComparisonV223.js'
    run=subprocess.run(['node',str(test)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=15)
    (out/'frontend_test.txt').write_bytes(run.stdout)
    assert run.returncode==0,run.stdout.decode('utf8',errors='replace')
    save(out/'report.json',dict(passed=True,faces=n,areaM2=float(sum(model['faceArea'])),fields=len(model['fields']),
        updates=[s+1 for s in STEPS],exportSHA256=sha(target),imageSHA256=sha(out/'rocket_phase_comparison.png'),
        sourceFaceCSV_SHA256=sources,dataRoundTripExact=True,sharedFrontendSHA256=sha(ROOT/'web/viewer.js'),
        browserInteractionTested=False,sharedGUIChanged=False,physicalQualification=False))
    print(json.dumps(dict(passed=True,faces=n,fields=len(model['fields']),physicalQualification=False)))

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='task':main()
    else:
        guard.BASE=BASE;guard.CASE=BASE;guard.DOC=DOC
        actual=guard.size_bytes;guard.size_bytes=lambda unused:actual(BASE)
        locks=[Path(__file__).resolve(),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',
               BASE/'audit02/report.json',ROOT/'tests/testSU2PhaseComparisonV223.js',ROOT/'web/viewer.js',
               OLD/'geometry01/cache.npz',OLD/'workbench01/rocket_phase_research.json']
        locks += [OLD/'audit01'/('original_faces_%06d.csv'%s) for s in STEPS]
        raise SystemExit(0 if guard.supervise('comparison',[sys.executable,str(Path(__file__).resolve()),'task'],40,256,locks,BASE,'offline') else 1)
