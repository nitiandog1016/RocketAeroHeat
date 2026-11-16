"""模块职责：真实完整循环更新原面研究场接既有WebGL，统一色标、方法和未知温度身份。"""
from pathlib import Path
import csv,json,re,subprocess,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm,SymLogNorm,Normalize
import runNativeILUV251 as run
from terminalCapacityV249 import supervise
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=run.BASE;DOC=run.DOC;AUDIT=BASE/'audit01';CLOSED=True
GEOMETRY=ROOT/'outputs/v222_validation/geometry01/cache.npz'


def export():
    out=BASE/'workbench02';out.mkdir(exist_ok=False);audit=read(AUDIT/'report.json');assert audit['passed']
    g={k:v for k,v in np.load(str(GEOMETRY)).items()};frames={s:np.loadtxt(str(AUDIT/('original_faces_%06d.csv'%s)),delimiter=',',skiprows=1) for s in audit['completeSavedUpdates']}
    recovery=read(BASE/'recovery01/report.json');assert recovery['passed'] and not recovery['nativeFullFrame']
    recovered=recovery['update'];recovered_path=BASE/'recovery01'/('original_faces_%06d.csv'%recovered)
    assert recovered not in frames and sha(recovered_path)==recovery['faceCSV_SHA256']
    frames[recovered]=np.loadtxt(str(recovered_path),delimiter=',',skiprows=1)
    sources={s:('pre-update-old-six-state-reconstruction' if s==recovered else 'native-full-34-column-frame') for s in frames}
    stages=sorted(frames);last=stages[-1]
    # 已完成次数与已保存帧是两种身份，异常退出时不能用较早保存帧冒充末态。
    stops={'linear-tolerance-rejection':'更新前线性精度拒绝（未应用该方向）','wallclock-reserve-stop':'时间预算结束','significant-deterioration':'窗口恶化停止',
      'two-windows-no-benefit':'两窗无收益停止','maximum-updates':'更新上限',
      'native-convergence-stop-unqualified':'原生停止条件触发，未授物理资格',
      'native-domain-or-runtime-failure':'数值或运行异常停止'}
    stop_label=stops[audit['stopReason']]
    if audit['recoveryRoles'].get('jv'):stop_label='差分试探越域，本分支停止'
    model=read(ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json');old={f['key']:f for f in model['fields']};fields=[];ranges={}
    method='SU2 8.6 / 原生Newton-Krylov，原生FGMRES / 块ILU(0) / 更新前原生容差检查 / 相态Air线性输运 / Roe对流通量/原生黏性RANS / SA-neg WITHFT2'
    scope='原19752面，M5、30 km、α20°、给定Tw300 K；从v248200六量建立的新方法有限发展，未收敛、非热资格。'
    specs=[('pressure','壁面压力','Pa',2,'当前壁节点primitive到原三角面顶点均值'),
      ('heatFlux','入壁对流热流','W/m²',3,'分子k×(当前法向邻点温度−给定Tw)/原法向距离；原面顶点均值，正值入壁'),
      ('wallFluid','壁节点流体温度','K',4,'当前fine primitive；不是固体壁温'),
      ('wallDeviation','壁流体温度减给定Tw','K',None,'壁节点流体温度减300 K；不是已解固体温度或热平衡误差')]
    for key,label,unit,col,detail in specs:
        vals={s:a[:,col] if col is not None else a[:,4]-a[:,5] for s,a in frames.items()};lo=min(float(min(v)) for v in vals.values());hi=max(float(max(v)) for v in vals.values())
        if key=='wallDeviation':lo,hi=-max(abs(lo),abs(hi)),max(abs(lo),abs(hi))
        for s,v in vals.items():
            name=key+'_'+str(s);ranges[name]=[lo,hi]
            fields.append(dict(key=name,label=label+' · NK更新'+str(s),unit=unit,values=v.tolist(),status=[1]*19752,
                methodIndex=[0]*19752,methodNames=[method+'；'+detail],scope=scope,computedStateLabel='未收敛研究数值',
                stateNote=('更新%d旧六量独立反算；不是原生完整末场，热流为法向差分研究诊断；不是物理时间。'%s if s==recovered else '真实NK更新%d；不是物理时间。初始化值未经过壁条件收敛；研究覆盖不等于热资格。'%s),
                sourceKind=sources[s],successfulUpdates=s,comparisonVariable=key,availability='available'))
    for key in ('wallInput','wallSolved'):old[key]['scope']=scope;fields.append(old[key])
    model.update(fields=fields,comparisonColorRanges=ranges,defaultField='heatFlux_'+str(last),
      resultCondition=dict(caption='M5 · α20° · 给定Tw300 K · 已完成%d次，完整帧%s，旧态恢复%d · %s · 共同色标'%(audit['completedUpdates'],audit['completeSavedUpdates'],recovered,stop_label)),
      caseNote='v251原真实箭体原生Newton-Krylov有限发展。'+stop_label+'；实际完成%d次，原生完整保存帧%s，更新%d来自旧六量恢复。'%(audit['completedUpdates'],audit['completeSavedUpdates'],recovered)+'压力/热流均为未收敛研究；给定Tw、壁流体T和未知固体Tw分开。尚无新增合格热面积。',
      researchIdentity=dict(initialization='v248200-fixed-native-NK-ILU0',physicalQualification=False,productionQualified=False,successfulUpdates=stages,
        nativeCompleteSavedUpdates=audit['completeSavedUpdates'],reconstructedUpdates=[recovered],surfaceSources=sources,recoveryAudit=recovery,
        sourceFaceCSV_SHA256={str(s):sha(recovered_path if s==recovered else AUDIT/('original_faces_%06d.csv'%s)) for s in frames},auditSHA256=sha(AUDIT/'report.json'),
        linearAcceptance=audit['nativeLinearHistory'],linearRefusal=audit['linearRefusal'],workingArrayIsFinalResidual=False,terminalResidual=audit['terminal'],nativeDevelopmentAudit=dict(actualFactory=audit['actualFactory'],observedMatrixFreeProducts=audit['observedMatrixFreeProducts'],completedUpdates=audit['completedUpdates'],recoveryRoles=audit['recoveryRoles'],firstUpdateBitEqual=audit['firstUpdateBitEqual'],resourceContractPassed=audit['resourceContractPassed'],stopReason=audit['stopReason'])))
    payload=json.dumps(model,ensure_ascii=False,separators=(',',':'),allow_nan=False);target=out/'rocket_newton_development.json';target.write_text(payload+'\n',encoding='utf8')
    (out/'case-data.js').write_text('window.ROCKET_CASE='+payload+';\n',encoding='utf8')
    html=(ROOT/'web/index.html').read_text(encoding='utf8');html=re.sub(r'((?:src|href)=")((?:viewer|section|export)\.(?:js|css)(?:\?[^\"]*)?)(")',r'\1../../../web/\2\3',html)
    (out/'index.html').write_text(html,encoding='utf8');restored=read(target)
    assert restored['faces']==g['originalFaces'].tolist() and restored['vertices']==g['originalVertices'].tolist()
    for field in restored['fields'][:-2]:
        a=frames[field['successfulUpdates']];col={'pressure':2,'heatFlux':3,'wallFluid':4}.get(field['comparisonVariable']);assert np.array_equal(field['values'],a[:,col] if col is not None else a[:,4]-a[:,5])
    # 模块二：首末完整保存帧按同一色标比较，明确末帧与实际完成更新数的区别。
    tri=g['originalVertices'][g['originalFaces']];right=np.array([1.,0.,.1]);right/=np.linalg.norm(right);up=np.array([-.025,1.,.25]);up-=right*(up@right);up/=np.linalg.norm(up)
    proj=tri@np.column_stack((right,up));order=np.argsort(tri.mean(axis=1)@np.cross(right,up))
    fig,axs=plt.subplots(3,2,figsize=(13,7.5));fig.subplots_adjust(left=.05,right=.86,top=.87,bottom=.11,hspace=.32,wspace=.18)
    compare=[stages[0],stages[-1]]
    for row,(key,col,scale,title,cmap) in enumerate((('pressure',2,1e-3,'Pressure [kPa], log','viridis'),('heatFlux',3,1e-6,'Heat into wall [MW/m2], symlog','coolwarm'),('wallDeviation',None,1,'Wall fluid T - given Tw [K]','coolwarm'))):
        plotted=np.concatenate([(frames[t][:,col] if col is not None else frames[t][:,4]-frames[t][:,5])*scale for t in compare]);lo,hi=float(min(plotted)),float(max(plotted))
        if hi==lo:hi=lo+1
        norm=LogNorm(lo,hi) if key=='pressure' else SymLogNorm(.01,vmin=-max(abs(lo),abs(hi)),vmax=max(abs(lo),abs(hi))) if key=='heatFlux' else Normalize(-max(abs(lo),abs(hi)),max(abs(lo),abs(hi)))
        for colplot,step in enumerate(compare):
            ax=axs[row,colplot];a=frames[step];v=(a[:,col] if col is not None else a[:,4]-a[:,5])*scale
            ax.add_collection(PolyCollection(proj[order],linewidths=0,facecolors=plt.get_cmap(cmap)(norm(v[order]))));ax.autoscale_view();ax.set_aspect('equal',adjustable='datalim');ax.set_axis_off()
            ax.set_title(title+' | update '+str(step)+(' (old-state recovery)' if step==recovered else ''),loc='left',fontsize=10)
        box=axs[row,1].get_position();cax=fig.add_axes([.9,box.y0+.1*box.height,.012,.8*box.height]);fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
    fig.suptitle('RocketAeroHeat v251 | Original 19,752 faces | M5 / alpha 20 deg / given Tw 300 K\nNative Newton-Krylov + ILU(0): bounded development from v248200, NOT a qualified thermal solution',fontsize=12)
    fig.text(.05,.035,'Phase-aware Air / JFNK / ILU0 + Roe / SA-neg WITHFT2 | Update 25: old-six-state recovery, no native final field\nSigned diagnostic heat retained. Given wall T, fluid wall T and unknown solid T remain separate.',fontsize=9)
    fig.savefig(str(out/'native_newton_development.png'),dpi=160);plt.close(fig)
    # 模块三：各式按自身首基准RMS归一化，不将线性残差或Jv扰动混入趋势。
    data=np.loadtxt(str(AUDIT/'baseline_trends.csv'),delimiter=',',skiprows=1,ndmin=2)
    fig,ax=plt.subplots(figsize=(10,4.8));fig.subplots_adjust(left=.1,right=.96,bottom=.18,top=.82)
    for j,label in enumerate(('mass','x momentum','y momentum','z momentum','physical energy','SA (separate after-time stage)')):
        ax.semilogy(data[:,0],np.maximum(1e-300,data[:,j+1]/max(1e-300,data[0,j+1])),label=label,lw=1.3)
    ax.axhline(1,color='gray',lw=.7);ax.set_xlabel('Completed nonlinear updates (not physical time)');ax.set_ylabel('RMS / first baseline RMS')
    ax.grid(True,alpha=.2);ax.legend(fontsize=8,ncol=2)
    fig.suptitle('RocketAeroHeat v251 | Native NK baseline flow + staggered SA\nUnconverged research; derivative probes excluded from residual history',fontsize=11)
    fig.text(.1,.04,'Each equation normalized by its own first value; not the absolute -8 qualification test.',fontsize=9)
    fig.savefig(str(out/'nk_residual_trends.png'),dpi=160);plt.close(fig)
    # 模块四：更新前的精度门单独展示；该量是原生缩放报告，不混同物理残差或实验误差。
    with (BASE/'rocket01/linear_acceptance.csv').open(encoding='utf8') as f:
        gates=[a for a in csv.DictReader(f) if a['system']=='flow']
    x=np.array([int(a['update']) for a in gates]);y=np.array([float(a['scaledReport']) for a in gates])
    accepted=np.array([a['accepted']=='1' for a in gates]);iters=np.array([int(a['iterations']) for a in gates])
    fig,axs=plt.subplots(2,1,figsize=(10,6),sharex=True);fig.subplots_adjust(left=.12,right=.96,bottom=.17,top=.84,hspace=.18)
    axs[0].semilogy(x,np.maximum(y,1e-300),color='#356b8a',lw=1.2,label='Native scaled Krylov estimate')
    axs[0].axhline(1e-5,color='#333333',ls='--',label='Configured acceptance: 1e-5')
    if np.any(~accepted):axs[0].scatter(x[~accepted],y[~accepted],color='#b43131',s=45,label='Refused before flow update',zorder=3)
    axs[0].set_ylabel('Native residual estimate');axs[0].legend(fontsize=9)
    axs[1].plot(x,iters,color='#356b8a',lw=1.2);axs[1].axhline(20,color='#333333',ls='--')
    axs[1].set_xlabel('Attempted nonlinear update (not physical time)');axs[1].set_ylabel('Outer FGMRES iterations')
    for ax in axs:ax.grid(True,alpha=.2)
    fig.suptitle('RocketAeroHeat v251 | Fixed FGMRES + block ILU(0)\nLinear acceptance is separate from nonlinear / thermal qualification',fontsize=12)
    fig.text(.12,.045,'Report = estimated relative residual / native tolerance factor (1 to 2).\nNo independent full-Jacobian accuracy or converged thermal qualification is claimed.',fontsize=9)
    fig.savefig(str(out/'linear_acceptance.png'),dpi=160);plt.close(fig)
    result=subprocess.run(['node',str(ROOT/'tests/testNativeILUV251.js')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
    (out/'frontend_test.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout.decode('utf8',errors='replace')
    save(out/'report.json',dict(passed=True,fields=len(fields),faces=19752,updates=stages,exactResearchValuesChecked=19752*4*len(stages),
      originalAreaM2=float(sum(g['faceArea'])),nativeCompleteSavedUpdates=audit['completeSavedUpdates'],reconstructedUpdates=[recovered],linearAcceptanceRecords=len(gates),physicalQualification=False,solidTwSolved=False,sharedWebGLChanged=False,browserInteractionTested=False))


def main():
    if CLOSED:raise SystemExit('v251 native development export closed; follow NEXT_DECISION')
    g=run.setup_guard()
    locks=[Path(__file__).resolve(),Path(run.__file__),ROOT/'scripts/runSU2RocketPhaseV222.py',DOC/'PLAN.md',GEOMETRY,AUDIT/'report.json',
      ROOT/'outputs/v222_validation/workbench01/rocket_phase_research.json',ROOT/'tests/testNativeILUV251.js',BASE/'rocket01/stages.csv',BASE/'rocket01/linear_acceptance.csv',ROOT/'scripts/terminalCapacityV249.py']+list(AUDIT.glob('original_faces_*.csv'))
    locks += [p for p in (ROOT/'web').iterdir() if p.name in ('index.html','viewer.js','viewer.css','section.js','export.js')]
    locks += [BASE/'recovery01/report.json']+list((BASE/'recovery01').glob('original_faces_*.csv'))
    return supervise(g,BASE,'export2',[sys.executable,str(Path(__file__).resolve()),'task'],45,1536,locks,BASE,'audit')


if __name__=='__main__':
    if CLOSED:raise SystemExit('v251 native development export closed; follow NEXT_DECISION')
    if sys.argv[1]=='task':export()
    else:raise SystemExit(0 if main() else 1)
