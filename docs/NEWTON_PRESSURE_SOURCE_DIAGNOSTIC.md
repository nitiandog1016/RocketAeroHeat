# v71：真实箭体Newton拒绝诊断

本阶段已完成软件测试和一次真实失败步的重演，未修复稳态收敛，也未增加物理解面积。测试工况仍为M5、30 km、零攻角完整箭体；三维横流、分离/再附、底热物理验证及全周期目标继续保留。

## 新程序与来源合同

新增完整文件 `regionalFlow/axisymmetric/advanceAxisymmetricSteadyObserved.m`，只用于研究诊断。默认求解器 `advanceAxisymmetricSteady.m` 与冻结v70保持逐字节一致，旧端点的源码绑定、读取审计、导出和demo都保持。这样可以观察旧初态，而不用改写旧链来源或重算一条链来更换SHA。

观察入口与原求解器前五个参数一致；独立第六参数指定 `enabled`（默认false）和 `maxRecords`（默认256）。第一输出是完整原result，第二输出是diagnostics。关闭时不保存候选。开启后每次真实线搜索保存迭代起Q、候选Q、无量纲增量、步/尝试/Newton/回溯序号、alpha、dt、原缺陷、原目标、候选缺陷或错误identifier、错误消息以及线性缺陷。步骤和Newton序号从1开始，attempt/search从0开始。

条目分为accepted、valid-not-decreasing和state-failure。后者的缺陷为NaN，不能当作零或正常残差。容量不足时保留最先发生的条目，同时记录seenRecords、droppedRecords和complete=false；maxRecords=0明确丢弃所有条目。诊断范围仅为线搜索，不包括Jacobian内部差分试探；没有新增替代残差、裁值或异常吞并。

观察入口保留冻结默认求解器的算法正文副本。`tests/testAxisymmetricNewtonObserverSource.py`移除明确列出的观察语句后，对算法正文逐行核对，任何其他变化均失败。未来默认求解器变化时必须重新审核该副本，不能自动宣称相容。运行测试另外比较完整result和调用计数，源结构相同不能替代实际结果测试。

```matlab
[result, diagnostics] = advanceAxisymmetricSteadyObserved( ...
    Q, mesh, gas, boundary, solverOptions, ...
    struct('enabled', true, 'maxRecords', 512));
```

## 已实际核对的结果

九格非均匀闭域四量BE测试中，默认求解器与观察器的关闭、开启、容量1和容量0五种结果完整相同，包括调用计数。原闭域质量约束和墙热保持；六项诊断/输入拒绝通过。

真实重演使用v70第七步的内部接受态作初猜，逐字节绑定原工况、汇总MAT和第八步MAT；不是给定边界，也不是已验证物理解。独立观察结果与冻结第八步的完整result逐值相同：2111次残差、一次伪步拒绝，随后接受dt=5.284620459572102e−9 s。旧链完全没有改写，未将重演计为第十一接受步或新生产解。

56个候选全部重新调用完整RANS核对：9个接受线搜索、47个状态失败、0个有效但不下降。47个错误identifier均为 `AxisFV:PressureSource`。最后失败的第6次Newton迭代中，alpha从1至2^−14的15个候选均失败；其原BE缺陷0.35415599024678257，目标1.0323953677463062e−9。候选实际到达源压力护栏，不能将这一次拒绝解释为单元温度超过2000 K或仅为Armijo下降不足。

## 独立定位：底面第108格的源压力重构

原空间算子在平面质心计算轴对称几何源，压力由r加权体质心的单元压力加WLS梯度线性外推：

`pA = pCell + grad(p) · (planarCentroid - volumeCentroid)`。

独立Python直接从实际面位移重新组装镜像WLS，不读取MATLAB保存的梯度矩阵；单独实现EOS及开放HLLE边界质量方向。对两个已保存完整终态，其pA与原info最大绝对差4.37e−11 Pa。

56个候选的单元守恒状态均有限、密度/压力为正、SA变量非负，单元温度均处在原40～2000 K数值温域内。47个失败候选均出现非正重构源压力。最后15个候选中，最小值均位于真实底面第108格、第一径向层，体质心约(x,r)=(-1.00035522e−5, 0.07511633) m。该格平面质心与体质心的位移约(2.22207e−10, -1.11326e−4) m，首层偏斜比约0.52149。

最后一个候选alpha=2^−14时：单元压力仍为161.3665026 Pa，外推源压力为−0.00068284 Pa；所有单元温度范围约226.509～1765.600 K。该候选包含一个非正pA控制体。负源压力足以触发实际护栏；Python没有重新实现全部RANS通量，不能据此声称已经独立验证所有通量有限性或全部失败机制。

这是对第八步一次拒绝的直接定位，不自动解释v70其余四次拒绝，不证明所有不稳态都来自该处。未改变压力源、线搜索次数、差分步长、温域或任何接受门槛。

## 保留的工作台及验证边界

原十步链已重新逐端点核对；`demoRocketAxisymmetricStudy`在当前源码下重新生成五步骤研究工作台成功。Node源面/气体温度代理/周向曲线/CSV及同宿主事件测试通过。默认main流程不变；没有新增浏览器视觉或WebGL运行验收。

运行：

```matlab
addpath('tests', 'regionalFlow/axisymmetric');
testAxisymmetricNewtonDiagnostics;
testRocketCollarNewtonDiagnostics;
```

工程包内真实前缀在 `validation/v70/normal_collar`，第八端点在 `validation/v70/normal_collar_continued/steps/step_000008.mat`。包内重演应将这两个位置显式传给第二个测试，输出到新的目录。已有结果不可覆盖。原型原位加记录的源码和实际测试记录单独保留在source_epochs；最终交付使用独立观察入口。首次拆分脚本因本机旧Python不支持write_text的newline参数而停止，随后测试因观察器文件尚未生成而exit1；修正脚本后完整复测exit0，失败记录没有删除。

首次包内隔离的MATLAB、两个Node和Python源码审计均实际exit0；最终比较器对MATLAB空字符数组错误调用数值NaN比较，触发Python类型错误，隔离整体仍判失败。原ZIP、发布/失败报告和比较器原源码保存在v71_attempt1_empty_character_comparison及validation/v71/attempts/isolation_empty_character，未将其标为交付成功。比较器现按字符/数值类型分开处理；这不改变候选、求解器、容差或任何物理结果。最终交付以修正后重新归档及隔离报告的实际成功为准。

## 下一步

先比较有数学一致性和适用依据的压力源积分方案，检查实际轴/底面/压力通量与几何源的相容性，包括常压平衡、仿射压力制造场、网格收敛及五量积分账。若引入有界重构，必须明确它改变的离散合同和精度阶，保留原默认路径，不能对负pA事后裁零。用新源、新研究链作真实同初态对照，并独立检查接受域、步长和稳态残差，不能用延长回溯、热流变色或守恒闭合替代物理解。

稳定零攻角完整外域及网格/远场/公开验证后，继续真正20°横流、收缩—后体、分离再附、底面和突起；累计热量仍后置。生产合格热流面积仍66.5626327424%，新增合格面积0；小版本不推送GitHub。
