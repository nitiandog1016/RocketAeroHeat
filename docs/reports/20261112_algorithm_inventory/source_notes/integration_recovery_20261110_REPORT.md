# 三候选独立恢复轮：可交主线审查

## 决定

**集成验收通过，建议“可交主线审查”。** 仅指当前依赖/制造夹具及原判据下的软件组合；不是生产应用或自动合并授权，不新增真实箭体/壁热流/CFD入口物理资格。

新轮只启动1次MATLAB，exit0、timeout=false，**27/27项原生记录通过**；原生入口5.793102300秒、父监督累计18.105345100/90秒，全部数值仅本次MATLAB，未用第二次机会。旧轮唯一8.672643200秒/exit1保持旧失败身份，不改旧监督/准入/回执/报告。父前后有效进程检查为空，不复用或干扰其他会话。最早一次只读进程探针因外层PowerShell插值报.Name错误（exit0不当作成功），已修正后空查询，无MATLAB启动；详见[准备观察](<harness_tasks/integration_recovery_20261110/preparation_observations.json>)。

## 完整交付

- [完整修复验证入口](<harness_tasks/integration_recovery_20261110/integration_validate.m>)、[入口相对旧修复稿的差分](<harness_tasks/integration_recovery_20261110/entry_repair.patch>)、[事前判据/资源限制](<harness_tasks/integration_recovery_20261110/PLAN.md>)、[新轮独立监督器](<harness_tasks/integration_recovery_20261110/run_recovery.py>)。
- 三完整未改候选：[calcBoundaryLayerState.m](<harness_tasks/integration_recovery_20261110/candidates/calcBoundaryLayerState.m>)、[calcTransition.m](<harness_tasks/integration_recovery_20261110/candidates/calcTransition.m>)、[marchHeadBoundaryLayerSegment.m](<harness_tasks/integration_recovery_20261110/candidates/marchHeadBoundaryLayerSegment.m>)；[原三文件组合patch](<harness_tasks/integration_recovery_20261110/combined.patch>)。patch仅git apply --check再次exit0，从未应用。
- [逐项动态结果](<harness_tasks/integration_recovery_20261110/itemized_results.csv>)、[状态汇总](<harness_tasks/integration_recovery_20261110/state_summary.csv>)、[1014行逐逻辑面q/Cf/theta/H及finite/unknown/来源](<harness_tasks/integration_recovery_20261110/per_face_results.csv>)、[A/B差异](<harness_tasks/integration_recovery_20261110/AB_comparison.csv>)、[接缝和单/双输出诊断](<harness_tasks/integration_recovery_20261110/seam_summary.json>)。
- [34项which/函数句柄/SHA](<harness_tasks/integration_recovery_20261110/function_resolution.csv>)、[热主链限定scope的profile计数](<harness_tasks/integration_recovery_20261110/main_chain_profile.csv>)、[实际解析与执行依赖](<harness_tasks/integration_recovery_20261110/dependencies.json>)；快照含预存根目录/测试源，**存在于快照不等于执行**，实际执行另由profile列出。
- [新准入](<harness_tasks/integration_recovery_20261110/evidence/matlab01/admission.json>)、[新回执](<harness_tasks/integration_recovery_20261110/evidence/matlab01/receipt.json>)、[实际执行前源/输入快照索引](<harness_tasks/integration_recovery_20261110/evidence/matlab01/source_input_manifest.json>)、[实际生成输入MAT](<harness_tasks/integration_recovery_20261110/evidence/matlab01/manufactured_inputs.mat>)、[完整原生结果](<harness_tasks/integration_recovery_20261110/evidence/matlab01/native_results.json>)、[原stdout字节](<harness_tasks/integration_recovery_20261110/evidence/matlab01/stdout.bin>)、[严格等值UTF8查看副本](<harness_tasks/integration_recovery_20261110/evidence/matlab01/stdout_utf8.txt>)、[原日志](<harness_tasks/integration_recovery_20261110/evidence/matlab01/matlab.log>)。各阶段5份完整MAT状态位于同轮目录。

## 只修入口，当前主链不回退

restoredefaultpath在**非嵌套**configurePrivatePath局部函数内执行；嵌套note不再令该脚本进入静态工作区。新原生记录实际which=[restoredefaultpath.m](<D:/matlab/toolbox/local/restoredefaultpath.m>)，并记录其SHA和MATLAB版本9.14.0.2206163 (R2023a)。额外只做输出目录校验、阶段/部分结果证据及原判据失效即停。

fixtures、logicalMesh、validateState、snapshot、arrayDiff整个科学区块SHA为`60f203cd32b55c956166e52087037d7d0ebeb0bacf980ad60439244fb5742cf2`，与旧提案逐字节相同；热合同断言子集和接缝输入SHA也保持。原4个小测试在原路径完整执行、断言未改。含真实STL/重网格/UI的完整热流测试从未运行。

当前[calcSurfaceHeating](<calcSurfaceHeating.m>) SHA=`2520252095d438fda56323b854ade02641aae3c58336893731a2f4b5e12a2f3d`，含v224显式相态提前分派，**未降回旧版本**。本轮仍只选原boundaryLayerIntegral/thermallyPerfect/referenceEnthalpy/head，不执行或修改v224相态外流生产器。

| 未改候选 | SHA256 |
|---|---|
| [calcBoundaryLayerState.m](<harness_tasks/integration_recovery_20261110/candidates/calcBoundaryLayerState.m>) | `052a8edc9f6e6b5813a9c8f4bcee2f8be210040224b0e46ebf2a22d429e8a0f1` |
| [calcTransition.m](<harness_tasks/integration_recovery_20261110/candidates/calcTransition.m>) | `016214e2042b033356777093431bff6624e792dac2c8aa4b9f80d2e2bb3347c4` |
| [marchHeadBoundaryLayerSegment.m](<harness_tasks/integration_recovery_20261110/candidates/marchHeadBoundaryLayerSegment.m>) | `ce71bf8c6c3b79c75041d4bfc6e59ce883051d9ad2e620aff4ea377ae62d58cf` |

## 热流主链真实执行，不以组件直调替代

production在新任务中性cwd，candidate在candidates cwd；当前目录优先覆盖只限3候选。阶段切换clear三候选、原热主入口及两个调用方，逐一检查which与函数句柄file。uniform的profile严格只围住一次calcSurfaceHeating，早于任何组件直调/原测试：

| 函数 | 生产A调用数 | 候选B调用数 |
|---|---:|---:|
| calcSurfaceHeating | 1 | 1 |
| calcBoundaryLayerState | 1 | 1 |
| calcTransition | 3 | 3 |
| calcTurbulentBoundaryLayerState | 1 | 1 |
| calcFixedShapeTurbulentBoundaryLayerState | 1 | 1 |
| marchHeadBoundaryLayerSegment | 93 | 93 |

三候选实际FileName和SHA均对应候选目录；热流主入口及其他依赖仍来自当前生产路径。三函数直接调用仅作后续守卫/锁定接缝补充，不冒充主链集成证明。完整实际来源、标签/分类/unknown掩码符合原A/B判据。

## q/Cf/厚度与unknown身份

下表为候选B；原A相应掩码/来源相同。都是制造逻辑面，F记录重复，不是实际STL/真实拓扑或原箭体面积。选用Cf/theta仅按实际层/湍流输出区填写，非输出区保持unknown。独立H仅对湍流状态有限，层流没有被伪造H。

| 夹具 | q有限/总数 | 选用Cf有限 | 选用theta有限 | Head H有限 | 原sourceID计数 |
|---|---:|---:|---:|---:|---|
| uniform | 100/101 | 100/101 | 100/101 | 93/101 | {'0': 8, '7': 1, '2': 92} |
| wall_gap | 99/101 | 99/101 | 99/101 | 92/101 | {'0': 9, '7': 1, '2': 90, '3': 1} |
| external_gap | 18/101 | 18/101 | 18/101 | 11/101 | {'0': 90, '7': 1, '2': 10} |
| forced_turbulent | 100/101 | 100/101 | 100/101 | 100/101 | {'0': 1, '4': 2, '7': 1, '2': 97} |
| deceleration | 39/103 | 39/103 | 39/103 | 39/103 | {'4': 4, '7': 3, '2': 32, '5': 64} |

[逐面表](<harness_tasks/integration_recovery_20261110/per_face_results.csv>)单位q=W/m²、Cf无量纲、theta=m、H无量纲，unknown为CSV空值/JSON null并另给原isfinite/isnan布尔；不是0。完整MAT保留其余热/边界层字段。

| A/B夹具 | 最大q差 W/m² | 最大选用Cf差 | 最大选用theta差 m | 最大H差 |
|---|---:|---:|---:|---:|
| uniform | 0 | 0 | 0 | 0 |
| wall_gap | 0 | 0 | 0 | 0 |
| external_gap | 0 | 0 | 0 | 0 |
| forced_turbulent | 0 | 0 | 0 | 0 |
| deceleration | 0.00560724762 | 3.849984e-08 | 1.53346311e-09 | 1.56843404e-05 |

前4组的有限值满足原atol=1e-12、rtol=1e-10及所有来源/掩码相同。deceleration仅按原计划比较语义和来源；跨接缝有限值有预期修正，不套前4组同值尺度，原结果差完整保留，未调判据。

## 下游传播、局部启动和回退是不同身份

- wall_gap face20 q=unknown、Head Cf/theta未知、保留湍流carrier；face21 source=3（局部平板启动），q=1079.40007498 W/m²恢复有限。原层流Thwaites厚度20后仍unknown，不虚称恢复连续层流积分。
- external_gap face20直接外缘失效，face21虽当地Re<5e5仍下游阻断；q20及全部后继unknown，不按当地低Re重启，nose不填补。
- forced_turbulent观察source4 LT域外平板回退，并与source3无连续前驱时局部启动、source2 Head接续和source7域边混合分别保留标签；未改任何带宽/模型。
- deceleration从face39疑似Head分离，直到face102全部阻断，经过速度恢复和face70壁温空白也不再附；face70/71 source均5。独立分支103的q=967.776130323 W/m²仍有限，数值失效source6未触发，不冒充覆盖。
- 主链非法等成本前驱A/B均由更早RocketAeroHeat:Continuum:Path拒绝。组件补充5种畸形在候选分别得到原登记领域错误ID，原生产接受/异常行为如实保留在原生结果，不把更早主链拒绝算成候选新增功劳。

## 接缝误差和诊断开销

读取旧已独立求解参考H=1.600134154495439，没有新参考ODE。原|H误差|=1.91240657024e-05；候选=6.90558721317e-14，改善276936125.952倍，满足<=2e-6与>=10倍；status0/fraction1、原第一输出字段、单/双输出第一结果完全相同。候选双输出实际定位事件数=1、root trials=29，投影未用。

单输出profile仍执行trace追加(原候选74行)和RK part构造(138行)。单输出0.026489300秒、双输出0.025245600秒均单次带profile/JIT观测，**不是性能基准或惰性诊断实现**。本轮不机会性改nargout/候选算法。

## 主线审查建议、限制与未覆盖

建议审查这三份精确SHA及组合patch，对照当前生产目标基线重新检查上下文，审查来源标签/错误ID与上游调用者兼容性；如后续批准才另行应用，不能据本报告自动覆盖生产。保留当前热主入口，不将本legacy模式验收外推到v224新外流。

原限制保留：准二维工程积分/参考nu不等于完整可压缩边界层；无再附/真实三维分离/横流/激波干扰资格。H=1.6仅验证指定简单横穿，不保证任意ds/切触/步内多次过零。全主链source6数值失效、上Re守卫/上混合带、无路径/fixedShape/其他鼻模及v224分派没有新全面覆盖；原列明未达分支保留在原生结果，不为补覆盖增跑或换条件。

## 边界与归档

773份选定当前生产/原测试/规则/旧任务及旧失败证据SHA保持，Git HEAD和目录外状态保持；三候选/组合patch及执行源输入均不变。仅本新目录写入，旧轮没有第二次启动/计数重置或结果翻案。[最终只读审计](<harness_tasks/integration_recovery_20261110/archive_audit.json>)、[CSV逐值回读](<harness_tasks/integration_recovery_20261110/csv_readback.json>)、[交付SHA清单](<harness_tasks/integration_recovery_20261110/delivery_manifest.json>)可核；选择性SHA不等于全系统IO监控。

物理资格、壁热流资格、solverInletReady恒false，新增真实面积0。无生产patch/v224/UI/Git写、CFD/真实网格/代理/新参考。标准库仅导出已存数据，元数据归档另计，保守父+归档仍须<=120秒；CPU/RSS峰/费用未知，MATLAB整个安装未冻结，不据小控制赋予性能或资源资格。

## 仅元数据校正

审阅时发现MATLAB单个事件struct被JSON编码为对象：初稿把其8个字段误计为8个事件，现按原生记录纠正为1次低支→高支事件（29次root trial）；原生输出/候选/判据完全未改，不重跑MATLAB。CSV回读索引补入本已逐值回读的A/B表，现6表共1102行。初稿元数据副本保留，见[校正记录](<harness_tasks/integration_recovery_20261110/metadata_corrections.json>)。另附[相对原失败入口的完整差分](<harness_tasks/integration_recovery_20261110/entry_from_failed.patch>)，直接显示非嵌套路径恢复修复。
