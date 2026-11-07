# v178证据索引

事前范围见[PLAN.md](PLAN.md)，科学结论和接口说明见[交付文档](../../EILMER_NATIVE_INTERFACE_V178.md)，下次入口见[NEXT_DECISION.md](NEXT_DECISION.md)。本轮没有整箭非线性任务或新增验证面积。

- `evidence/runs/native01_prepare`：原CFL层级尚未正确的准备exit0，不等于配置合同通过。
- `evidence/runs/native02_geometry.json`及`native02`：初场dt未定义先触发拒绝；后续独立核实棱柱−0.625m³，不能继续求解。
- `evidence/runs/native03*`：楔体原序/正确Global CFL，原生准备和唯一最多4步小域运行；末残差未收敛。
- `evidence/runs/interface01`：真实小域回读、独立线性梯度正负零热流及13拒绝；首次只读输出封存时间点单列。
- `evidence/runs/fixture03_process`：旧true断言失败及准确源；`fixture04`是纠正后的独立几何/12拒绝。
- `evidence/runs/real01_prepare`、`real01_geometry.json`、`real01_wall.json`：真实完整域原生准备和原面审核，无真实新热流。
- `evidence/runs/supervision01`：超时与外部取消的真实子进程组结束检查，不是CFD。
- `evidence/runs/matlab01`与`eilmer_entry01_process.json`：MATLAB实际回读、两拒绝和源码锁定；父监督62.5秒。
- `evidence/reference`：复用已过的Run14参考整理及独立换算，本轮未重复执行。
- `evidence/source_versions`：本轮最终源码准确字节；初始/失败版本和实际runner也在各原回执目录保留。

[evidence_manifest.json](evidence_manifest.json)逐项列出小证据SHA和本地大场。准备/运行回执内路径是实际运行时的WSL路径，归档不改成假路径；跨机器复现须新建输入和回执，不能移动归档后称旧运行身份仍成立。源码输入与原生配置的身份绑定不授予物理资格。

大原生网格和初场只在`E:\RocketAeroHeat\outputs\v178_validation\real01`，没有备份GitHub。两份95,253,209字节网格已在SHA逐字一致后硬链接去重，两个路径和SHA保持，独有删除0。29项并发保护文件逐项未变。原v167大网格只读引用，不复制工程包。
