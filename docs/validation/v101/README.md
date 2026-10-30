# v101验证与交付索引

本阶段为已测试的局部移动激波/ALE内核及科研查看功能，尚未获得真实箭体共同激波、外缘、热流或壁温资格。正式热流覆盖保持66.5626327424%，新增0；总目标继续有效。交付以包含本记录的GitHub源码提交为准，不生成工程ZIP或重复隔离工程。

方法和使用范围见[SHOCK_FITTED_ALE_V101.md](../../SHOCK_FITTED_ALE_V101.md)。全部新增/修改MATLAB文件为项目中的完整源文件，使用中文模块注释。

## 已执行的检查

| 证据 | 实际内容 |
| --- | --- |
| `provenance.json` | 8项MATLAB任务全部终态；首轮预算拒绝测试exit1，其余exit0，所有实际源码快照在各次运行期间稳定；871项准确MATLAB/实际使用UI版本已绑定 |
| `shock_fitted_kernel_summary.json` | 六个移动RH、两物性活塞相容、常gamma三维8步、483原鼻部面/2898单元、严格拒绝及工作数组预算 |
| `independent_audit.json` | 独立Brent/NASA7、实验室xi=D Riemann、Gauss时空几何及全域五量账；没有把近零量的绝对舍入误当成物理误差 |
| `ui_export_independent.json` | 13,560个导出字段逐值核对，原V/F与面号绑定准确；真实箭体未知值全部保持空值 |
| `ale_ui_node_view05.json` | 10个流程事件、43,728面元探针与CSV行、14次CSV生成、10类无关/过时消息拒绝、8个几何子集/重编号核验 |
| `legacy_node.json` | 新通用几何显示下旧研究20流程事件、1092站点与CSV行、52次列导出、24类消息拒绝通过 |
| `browser/actual_browser.json` | 最终view05原生浏览器五步骤、初末帧/变量/范围切换和三维实际点选；活塞p=6788.1126Pa、流体T=428.0112K；真实鼻部面403对应原面19623，场量未知 |
| `browser/legacy_default.txt`、`.png` | 新renderer的旧无字段几何/法向列显示与η=0.6站点切换实际验收 |
| `audit_dependency_binding.json` | 独立热力审核复用Python依赖与实际源码时期/Git版本一致 |
| `research_sources.json` | 公开一手来源与实际阅读范围；本阶段无实验验证 |

原生浏览器曾触发CSV下载按钮，但Playwright下载捕获接口在该绑定中不可用，**未确认实际落盘文件**；43,728行CSV内容校验来自Node执行实际导出逻辑。浏览器没有宣称控制台完全干净。遇到浏览器接口问题后重置会话并使用原生控件恢复验收，没有将其改记为物理计算失败。

初始色标遮挡模型和未显示外壳顶点影响子集范围的问题已修正。最终截图为`browser/piston_final_layout.png`、`browser/real_nose.png`及`browser/full_rocket.png`；早期截图保留其原始阶段身份。工作台是保存研究查看流程，不是新增实时全面积主程序。

## 失败原文与源版本

- `raw/first_budget_refusal/`：首轮`nargin<9`覆盖显式预算，最后拒绝测试失败；改为正确八输入判断后整体通过，未改物理门槛。
- `raw/independent_swept_roundoff/`：恢复端点的近零扫掠固定绝对容差失败；改用实际面积/坐标尺度的舍入传播界。
- `raw/independent_zero_flux/`：近零能量通量逐点归一失败；改为独立逐守恒量通量范数，并保留逐分量绝对差，门槛2e-12。
- `raw/processes/`、`raw/logs/`：八项实际进程记录与输出原字节。
- `raw/UI_versions/`：真实用过的UI模板/脚本/样式版本。没有声称所有无关全工程快照字节均另行归档。

## 本地结果与复现

完整证据仅约1.84MB，保留在`outputs/v101_validation/kernel02/shock_fitted_kernel_evidence.mat`，没有上传被忽略的`outputs`。必要输入`tests/data/shock_layer_transport/parabolic_shock_layer_evidence.mat`已在Git，新增小型数值参考位于`tests/data/shock_fitted_kernel/reference.json`，明确标为数值回归数据。

MATLAB入口：`testShockFittedALEKernel(newFolder, sourceFile)`；独立审核入口：`tests/auditShockFittedALEKernel.py`；界面映射审核：`tests/auditShockFittedExport.py`；实际UI逻辑测试：`tests/testShockFittedStudy.js`。重算/重导出必须使用新目录和新任务名，绑定准确来源，不覆盖现有失败或证据。当前无需等待或重跑这八项任务。

最终查看地址为`http://127.0.0.1:8781/view05/workbench/index.html`。下一步为小规模真实曲面壁面—共同激波的下游特征相容、共享节点运动与非均匀稳定性，而非继续提高旧Cauchy阶数。真实20°鼻部的物理求解及各工程区域缺口仍需完成；原10µm体形门槛不变，NASA7 ALE未实现的部分不能用常gamma顶替。

## 阶段输出回收

最终视图与来源收口后，逐字节核对并回收view01—view04的36个中间工作台文件，文件长度约9.54MB，删除窗口磁盘空闲实增9,601,024字节（约9.60MB）。每个文件均有SHA完全一致的保留副本；独有小型UI字节单独存入`storage/retained_UI_bytes/`，所有数值字段已在view05完整保留。原八项进程记录、失败原文、1.84MB科学MAT与最终UI未删，独立后检查通过。旧view01—view04下载/浏览地址需要按`storage/plan.json`恢复资产，旧目录带`RETIRED_WORKBENCH.json`，不能再宣称旧URL可直接打开。本轮没有累积完整旧矩阵或新工程包；这只是约9.5MB的中间导出清理，不是新增GB级结果回收。
