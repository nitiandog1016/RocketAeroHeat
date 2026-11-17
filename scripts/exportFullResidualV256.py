"""模块职责：完整原面展示真实零时间残差/研究态及所有未知资格，不复制旧大模型。"""
import json,subprocess
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v256_validation';CLOSED=True


def export():
    audit=read(BASE/'audit01/report.json');native=read(BASE/'rocket02/native_completion.json');assert audit['passed']
    out=BASE/'workbench01';out.mkdir(exist_ok=False);a=np.load(str(BASE/'audit01/surface_values.npz'))
    fields=[]
    for key,label,unit,source,method in [
        ('nativePhysicalEnergyResidual','物理能量方程积分残差 · 节点均值','W','physicalEnergyResidual','成熟SU2同一状态完整空间残差；逐节点R_Ephysical=R_Einternal−300000R_rho后取原面三顶点均值'),
        ('nativeMassResidual','质量方程积分残差 · 节点均值','kg/s','massResidual','成熟SU2完整空间质量残差；原面三顶点均值；不是面通量或局部热流'),
        ('nativeSAResidual','SA方程残差 · 节点均值','原生单位','SAResidual','成熟SA-neg WITHFT2同状态空间残差；保留非守恒nu_tilde原生量纲'),
        ('nativeInputPressure','输入研究态压力 · 未收敛','Pa','pressurePa','v248/200输入经本轮物理能量坐标与原相态Air恢复；无非线性或时间推进'),
        ('nativeWallFluidT','壁节点流体温度 · 输入研究态','K','wallFluidTemperatureK','原热壁为弱能量边界；该值是节点流体T，不是给定Tw300 K或求解固体Tw')]:
        fields.append(dict(key=key,label=label,unit=unit,values=a[source].tolist(),method=method))
    data=dict(schema='rah-full-residual-v256',originalFaces=19752,originalAreaM2=1.9646292565715553,fields=fields,
        scope='原完整连通黏性域的零时间接口控制；旧200未收敛态只是输入。没有新气动热解/真实外缘/完整入口资格。',
        caseNote='v256 接口控制：同输入完整残差逐位重复，非法密度试探未污染状态；实际非线性/时间更新0。原19752面保留，所有诊断数值覆盖不计新增合格面积。3904重构警告原样保留；给定Tw300 K，固体Tw未知。',
        identity=dict(physicalQualification=False,newQualifiedAreaM2=0,profileReadyForCFD=False,nonlinearIterations=0,timeUpdates=0,
            wallTemperatureGivenK=300,solvedSolidWallTemperatureK=None,actualProfilesFilled=0,KINSOLExecuted=False,
            nativeSourceSHA256=sha(ROOT/'native/su2PhaseAwareAir/fullDomainResidualV256.cpp'),sourceStateSHA256=audit['inputSHA256']['outputs\\v256_validation\\rocket02\\initial_state.bin']))
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    (out/'residual.json').write_text(payload+'\n',encoding='utf8');(out/'residual-data.js').write_text('window.FULL_RESIDUAL='+payload+';\n',encoding='utf8')
    html=(ROOT/'outputs/v252_validation/workbench02/index.html').read_text(encoding='utf8')
    html=html.replace('<p style="margin:10px 20px"><a href="profiles.html">查看沿程热流与实际边界层剖面 →</a></p>',
        '<p style="margin:10px 20px;color:#485868">完整域零时间接口诊断 · 非线性推进0 · 新合格面积0 · 重构警告3904。<a href="../../v255_validation/workbench01/index.html">查看已有工程热候选</a></p>')
    html=html.replace('<script src="case-data.js"></script>','<script src="../../v252_validation/workbench02/case-data.js"></script><script src="residual-data.js"></script><script src="../../../web/full-residual-v256.js"></script>')
    (out/'index.html').write_text(html,encoding='utf8');assert read(out/'residual.json')==data
    plot(out,a,audit)
    p=subprocess.run(['node',str(ROOT/'tests/testFullResidualViewerV256.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(p.stdout);assert p.returncode==0,p.stdout
    save(out/'report.json',dict(passed=True,fields=len(fields),originalFaces=19752,originalAreaM2=1.9646292565715553,
        diagnosticValues=len(fields)*19752,roundtripExact=True,physicalQualification=False,newQualifiedAreaM2=0,newHeatSolutions=0,
        nonlinearIterations=0,timeUpdates=0,browserInteractionTested=False,sharedWebGLChanged=False,oldLargeModelCopied=False,
        dataSHA256=sha(out/'residual.json'),oldModelSHA256=sha(ROOT/'outputs/v252_validation/workbench02/case-data.js')))


def plot(out,a,audit):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import TwoSlopeNorm,Normalize
    geom=np.load(str(ROOT/'outputs/v253_validation/inlet02/inlet_profiles.npz'));V=geom['originalVertices'];F=geom['originalFaces'];tri=V[F];order=np.argsort(tri.mean(axis=1)[:,2])
    fig,axes=plt.subplots(2,1,figsize=(11,6.5),layout='constrained')
    for ax,key,title,unit,palette in [(axes[0],'physicalEnergyResidual','Physical-energy spatial residual: vertex mean','W','coolwarm'),
        (axes[1],'wallFluidTemperatureK','Wall-node fluid temperature of un-converged input','K','viridis')]:
        values=a[key];extent=float(np.max(abs(values)));norm=TwoSlopeNorm(vmin=-extent,vcenter=0,vmax=extent) if key=='physicalEnergyResidual' else Normalize(values.min(),values.max())
        ax.add_collection(PolyCollection(tri[:,:,0:2][order],facecolors=plt.get_cmap(palette)(norm(values[order])),edgecolors='none',rasterized=True));ax.autoscale_view();ax.set_aspect('equal');ax.set(title=title,xlabel='X [m]',ylabel='Y [m]');ax.grid(alpha=.15)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=palette),ax=ax,label=unit,shrink=.7)
    fig.suptitle('v256 | original 19752 faces | zero-time operator control | NO new qualified heat area')
    fig.savefig(out/'full_residual_diagnostic.png',dpi=160);plt.close(fig)


if __name__=='__main__':
    if CLOSED:raise SystemExit('v256 fixed residual workbench closed; follow NEXT_DECISION')
    export()
