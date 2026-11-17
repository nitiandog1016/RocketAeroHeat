"""模块职责：将真正原生返回态映射回全原面，提供压力/热流/节点温度及控制诊断。"""
import json,re,subprocess
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v257_validation';CLOSED=False


def export():
    audit=read(BASE/'audit01/report.json');n=audit['nativeSummary'];assert audit['passed']
    out=BASE/'workbench01';out.mkdir(exist_ok=False);a=np.load(str(BASE/'audit01/surface_values.npz'));fields=[]
    for key,label,unit,values,method in [
        ('kinsolReturnedQ','返回研究态对流热流 · 给定Tw300 K','W/m²',a['terminal_qConvGivenTwWm2'],
         'SU2同状态弱等温热壁：(k+Cp·mu_t/Pr_t)(T邻−Tw给定)/法向距离；原面三顶点均值，正值入壁。KINSOL首轮线性拒绝，返回原未收敛初态；不是新合格热解'),
        ('kinsolReturnedP','返回研究态压力 · 未收敛','Pa',a['terminal_pressurePa'],
         '原完整flow-SA同状态变比热Air物性；原面三顶点压力均值。KINSOL实际接受更新0，返回原未收敛初态'),
        ('kinsolWallFluidT','壁节点流体温度 · 非固体壁温','K',a['terminal_wallFluidTemperatureK'],
         '原完整域弱热壁的节点流体T；给定Tw=300 K另列，求解固体Tw仍未知。原面三顶点均值，未收敛研究态'),
        ('kinsolEnergyResidual','物理能量方程积分残差 · 节点均值','W',a['terminal_physicalEnergyResidual'],
         '同状态成熟空间算子：R_Ephysical=R_Einternal−300000R_rho；原面三顶点均值，不是热通量'),
        ('kinsolDeltaQ','返回态减初态 · 热流差','W/m²',a['terminal_qConvGivenTwWm2']-a['initial_qConvGivenTwWm2'],
         '同一原面、同一给定Tw下返回态与初态弱热壁诊断之差；真实零变化来自未接受更新，不是用零填补未知热流'),
        ('kinsolDeltaP','返回态减初态 · 压力差','Pa',a['terminal_pressurePa']-a['initial_pressurePa'],
         '同一原面的返回研究态压力减初态压力；真实零变化来自更新前线性拒绝，不代表流动已收敛')]:
        fields.append(dict(key=key,label=label,unit=unit,values=values.tolist(),method=method))
    identity=dict(actualObservedStateChanges=audit['actualObservedStateChanges'],nativeNonlinearAttempts=n['nativeNonlinearAttempts'],
        nativeFlag=n['nativeFlag'],nativeFlagName=n['nativeFlagName'],linearGateRefused=n['linearGateRefused'],
        physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,wallTemperatureGivenK=300,solvedSolidWallTemperatureK=None,
        nativeSourceSHA256=sha(ROOT/'native/su2PhaseAwareAir/kinsolCoupledV257.cpp'),sourceStateSHA256=sha(BASE/'rocket01/terminal_state.bin'))
    data=dict(schema='rah-kinsol-coupled-v257',originalFaces=19752,originalAreaM2=1.9646292565715553,fields=fields,identity=identity,
        scope='完整原域KINSOL控制的返回研究态。数据覆盖全部原面，但无新增合格面积、可信外缘或完整CFD入口。',
        caseNote='v257：1次非线性尝试，20次线性迭代，线性精度未达要求，实际接受更新0；返回原未收敛初态。新合格面积0，3904重构警告保留。给定Tw300 K，固体Tw未知。')
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    (out/'coupled.json').write_text(payload+'\n',encoding='utf8');(out/'coupled-data.js').write_text('window.KINSOL_COUPLED='+payload+';\n',encoding='utf8')
    html=(ROOT/'outputs/v252_validation/workbench02/index.html').read_text(encoding='utf8')
    html=html.replace('<p style="margin:10px 20px"><a href="profiles.html">查看沿程热流与实际边界层剖面 →</a></p>',
        '<p style="margin:10px 20px;color:#485868">KINSOL完整域控制 · 线性精度拒绝 · 接受更新0 · 返回原未收敛态。<a href="../../v255_validation/workbench01/index.html">已有工程热候选</a> · <a href="coupled_control.png">控制诊断图</a></p>')
    html=html.replace('<script src="case-data.js"></script>','<script src="../../v252_validation/workbench02/case-data.js"></script><script src="coupled-data.js"></script><script src="../../../web/kinsol-coupled-v257.js"></script>')
    (out/'index.html').write_text(html,encoding='utf8');assert read(out/'coupled.json')==data
    plot(out,a,audit)
    result=subprocess.run(['node',str(ROOT/'tests/testKinsolViewerV257.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout
    save(out/'report.json',dict(passed=True,fields=len(fields),originalFaces=19752,originalAreaM2=1.9646292565715553,
        valuesChecked=len(fields)*19752,roundtripExact=True,actualObservedStateChanges=identity['actualObservedStateChanges'],
        nativeFlag=n['nativeFlag'],physicalQualification=False,newQualifiedAreaM2=0,newHeatSolutions=0,
        browserInteractionTested=False,sharedWebGLChanged=False,oldLargeModelCopied=False,dataSHA256=sha(out/'coupled.json')))


def plot(out,a,audit):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import TwoSlopeNorm
    geom=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'));tri=geom['originalVertices'][geom['originalFaces']];order=np.argsort(tri.mean(axis=1)[:,2])
    fig,(ax,bx)=plt.subplots(2,1,figsize=(11,7),layout='constrained');q=a['terminal_qConvGivenTwWm2'];extent=float(np.max(abs(q)))
    norm=TwoSlopeNorm(vmin=-extent,vcenter=0,vmax=extent)
    ax.add_collection(PolyCollection(tri[:,:,0:2][order],facecolors=plt.get_cmap('coolwarm')(norm(q[order])),edgecolors='none',rasterized=True))
    ax.autoscale_view();ax.set_aspect('equal');ax.set(title='Returned research-state heat flux | 0 accepted updates | Tw given = 300 K',xlabel='X [m]',ylabel='Y [m]')
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='coolwarm'),ax=ax,label='q into wall [W/m²]',shrink=.65)
    text=(BASE/'rocket01/kinsol_info.log').read_text(encoding='utf8');pairs=re.findall(r'cur-iter = (\d+), total-iters = \d+, res-norm = ([0-9.eE+-]+)',text)
    points=np.array(pairs,dtype=float);target=float(audit['linearGates'][0]['target'])
    bx.semilogy(points[:,0],points[:,1],'o-',color='#32597a',markersize=3,label='Native SPGMR scaled linear residual')
    bx.axhline(target,color='#a44d43',linestyle='--',label='Fixed linear target');bx.set(xlabel='Krylov iteration',ylabel='Scaled residual norm',xticks=[1,5,10,15,20])
    bx.grid(alpha=.2);bx.legend(loc='center right');bx.set_title('20-vector control rejected before nonlinear update; no new qualified area')
    fig.savefig(out/'coupled_control.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v257 fixed coupled workbench closed; follow NEXT_DECISION')
    export()
