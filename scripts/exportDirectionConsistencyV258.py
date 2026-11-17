"""模块职责：固定方向诊断接回原面工作台，保留差分未通过身份与既有p/q来源。"""
import json,subprocess
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
BASE=ROOT/'outputs/v258_validation';CLOSED=True


def export():
    report=read(BASE/'audit01/report.json');assert report['passed'];out=BASE/'workbench01';out.mkdir(exist_ok=False)
    a=np.load(str(BASE/'audit01/surface_values.npz'));fields=[]
    # 审核载荷的energyTrueDirection键只表示完整F半步差分，不能因该键名宣称精确方向导数。
    for key,label,source,method in [
        ('directionEnergyHalf','物理能量完整F半步差分估计','energyTrueDirection',
         '固定未接受试探差d，(F_E(x0+sigma*d/2)−F_E(x0))/(sigma/2)，再取原面三顶点均值。差分步长一致性未通过，非精确Jacobian方向，非壁面热流'),
        ('directionEnergyApprox','物理能量近似块矩阵方向','energyApproxDirection',
         '同状态原生flow/SA块对角空间矩阵乘C*d，再C逆变换成物理能量坐标，原面三顶点均值。只用于预条件诊断，非完整Jacobian，非壁面热流'),
        ('directionEnergyDifference','近似矩阵减半步差分 · 能量方向差','energyDirectionDifference',
         '两个方向估计的物理能量分量之差，再取原面三顶点均值；因完整F差分未获一致性，不能把此差值当作近似矩阵对精确导数的误差，更不是热流增量')]:
        fields.append(dict(key=key,label=label,unit='W／方向参数',values=a[source].tolist(),method=method))
    identity={k:report[k] for k in ['fixedDirectionHalfStepConsistent','entireJacobianAccuracyProved','nonlinearUpdates','newQualifiedAreaM2','physicalQualification','profileReadyForCFD','sigma']}
    identity.update(KINSOLExecuted=False,terminalIsRestoredInitial=True,sourceDirectionSHA256=sha(BASE/'rocket01/saved_direction.bin'),
        sourceStateSHA256=sha(BASE/'rocket01/initial_state.bin'))
    data=dict(schema='rah-direction-consistency-v258',originalFaces=19752,originalAreaM2=1.9646292565715553,fields=fields,identity=identity,
        scope='原完整域固定失败方向的零更新诊断；半步差分未获一致性，不能用作可信Jacobian或外流/热流解。',
        caseNote='v258：固定方向sigma/half差分不一致，非精确导数；近似块矩阵差异单列。4次F评估后恢复初态，接受更新0、新合格面积0。已有p/q/T仍来自v257未收敛输入，给定Tw与未知固体Tw身份保持。')
    payload=json.dumps(data,ensure_ascii=False,separators=(',',':'),allow_nan=False)
    (out/'direction.json').write_text(payload+'\n',encoding='utf8');(out/'direction-data.js').write_text('window.DIRECTION_CONSISTENCY='+payload+';\n',encoding='utf8')
    html=(ROOT/'outputs/v257_validation/workbench01/index.html').read_text(encoding='utf8')
    html=html.replace('<script src="coupled-data.js"></script>','<script src="../../v257_validation/workbench01/coupled-data.js"></script>')
    html=html.replace('<script src="../../../web/kinsol-coupled-v257.js"></script>',
        '<script src="../../../web/kinsol-coupled-v257.js"></script><script src="direction-data.js"></script><script src="../../../web/direction-consistency-v258.js"></script>')
    html=html.replace('KINSOL完整域控制 · 线性精度拒绝 · 接受更新0 · 返回原未收敛态。','固定方向诊断 · 半步一致性未过 · 接受更新0 · 终态恢复原输入。')
    html=html.replace('href="coupled_control.png">控制诊断图','href="../audit01/direction_consistency.png">方向诊断图')
    (out/'index.html').write_text(html,encoding='utf8');assert read(out/'direction.json')==data
    result=subprocess.run(['node',str(ROOT/'tests/testDirectionViewerV258.js'),str(out)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=20)
    (out/'frontend_test.txt').write_bytes(result.stdout);assert result.returncode==0,result.stdout
    save(out/'report.json',dict(passed=True,newFields=3,originalFaces=19752,originalAreaM2=1.9646292565715553,valuesChecked=3*19752,
        roundtripExact=True,differenceIsNotExactJacobianError=True,existingPressureHeatFieldsUnchanged=True,nonlinearUpdates=0,
        physicalQualification=False,newQualifiedAreaM2=0,browserInteractionTested=False,sharedWebGLChanged=False,oldLargeModelCopied=False,dataSHA256=sha(out/'direction.json')))


if __name__=='__main__':
    if CLOSED:raise SystemExit('v258 fixed direction workbench closed; follow NEXT_DECISION')
    export()
