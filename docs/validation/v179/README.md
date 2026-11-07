# v179验证索引

[结果与限制](../../EILMER_RUN14_CONTROL_V179.md)、[原任务计划](PLAN.md)、[唯一控制预算](CONTROL_RUN_PLAN.md)、[下一步](NEXT_DECISION.md)。

- `evidence/control01_audit02.json`：实际气体模型、初态、物性及完整域几何审核。首次失败见audit01及准确旧源。
- `evidence/control01_run/report.json`：唯一原生任务、输入锁、最终输出SHA；`process00/stdout.txt`保留原生时间停止与831步。
- `evidence/comparison01/report.json`、`wall.csv`、`stations.csv`：原面壁量、逐测站有符号误差、来源与拒绝。全部研究身份。
- `evidence/stability02.json`：固定50步和最终31步窗口，区分前缘全局峰与扩张段峰。首次表头失败及准确旧源单列保留。
- `evidence/measurement_test02.json`：插值/未知/拒绝、切向应力独立旋转真值、全周面积与功率。
- `evidence/decimal01.json`：独立Decimal核对320原生面、76测站及全周热功率，未调用主回读器。
- `evidence/comparison01/comparison.png`及PDF：已目视核对的压力、全壁热流、扩张段热流、残差、剪切和y+科研图。
- `evidence_manifest.json`：逐字SHA、小证据与本地大场位置、29项并发保护核验。原参考表引用v178，不重复拷贝。

原论文和8个大二进制场文件留在`E:\RocketAeroHeat\outputs\v179_validation`，GitHub未备份，独有删除0。其余旧结果及原20°面积台账不变，没有本轮MATLAB/主UI/引擎构建或新增物理资格。
