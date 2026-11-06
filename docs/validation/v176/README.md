# v176 验证索引

[功能与边界](../../AEROTHERMAL_EVIDENCE_LEDGER_V176.md) · [事前计划](PLAN.md) · [断言恢复](RECOVERY_PLAN.md) · [实际探针修正](PROBE_REVIEW_PLAN.md) · [原生窗口观察](NATIVE_UI_OBSERVATION.md) · [下一决定](NEXT_DECISION.md)

| 检查 | 实际结果 | 时间 |
|---|---|---:|
| main01（含unit01） | unit通过；UI第43行range数组方向断言失败，完整失败保留 | 57.969s |
| main02（只恢复UI） | 工程/研究原面、数值、区域、固定范围、92段、旧帧拒绝、5CSV通过 | 179.109s |
| probe03（只检查实际发现的可读性改动） | 默认折叠后的真实WebGL、面9327/92段/渲染色标通过 | 109.921s |
| Node初版/折叠修正版 | 118/128项通过，最小DOM不冒称WebGL | 0.172/0.157s |
| Python独立Oracle | 2,147,005项通过，5CSV各19,752行、输入稳定 | 3.344s监督/3.135s程序 |

MATLAB三次均已终态，合计346.999s低于480s；main02/probe03分别含70.732/72.679s窗口观察等待，未从监督总时间扣除。每次源与两份原面MAT输入稳定，未重跑CFD、引擎构建、大体场或旧五步数值链。原生鼠标通过的是折叠修正前的区域→原面/原面→探针与92段曲线；修正后用户Esc停止，未补截图或鼠标，不冒称最终目视验收已完成。

合成6面真值面积1/2/4/8/16/32，p通过7、q通过3、交集3；18拒绝案例、19拒绝记录和3个标量给定Tw兼容案例通过。该合成资格不用于真实模型。真实研究适配只复用既有审核来源，不重新认证v167；两份真实结果的新增物理资格0。

准确执行源码按版本快照、失败日志、unit证明、原面集成与Oracle报告保存在`evidence/`；[清单](evidence_manifest.json)记录文件SHA及局部大输出清单。旧并发查看器源码仅记录SHA，不混入此次提交；全部29个保护文件核对保持。五份CSV/两份原面JSON和完整历史MAT只在本地，GitHub没有备份。

存储目标80MiB已超出：原CSV逐行完整方法范围使工程q导出单项38,619,325字节，全部本轮输出约99MiB（最终精确值见[交付核对](delivery_review.json)）。不删除独有结果、不生成工程包；probe03没有重复大导出。公开归档目标5MiB内，只有必要小证据。开发协调墙钟未计入程序执行时间。

独立复核入口：`tests/testAerothermalEvidenceLedger.m`可在新目录重建合成反例；`tests/testAerothermalEvidenceView.js`运行Node；真实集成需两份本地原面MAT。`tests/verifyAerothermalEvidenceLedger.py`接受unit夹具、工程/研究JSON及五份CSV，在新目录写独立报告。缺本地结果时不能把当前仓库说成可立即重放全部历史大计算。
