# v192证据索引

本轮完成Euler小控制的同场发展、独立细网格计算和空间偏差对照。粗细均达到事前数值收敛指标；膨胀温度误差未随加密降低，不能授予工程外缘或热流资格。结论见[实现报告](../../HYBRID_EULER_REFINEMENT_V192.md)，后续直接按[NEXT_DECISION](NEXT_DECISION.md)进入真实20°全局网格准备与外缘资格接口，不继续第三层统一加密或追加步数。

- 事前预算：[原发展计划](PLAN.md)、[唯一细网格计划](REFINE_PLAN.md)。原计划的保存误解保留，未按实际结果追改。
- [准确原生源码依据](evidence/native_source_evidence.json)、[粗任务复制与来源](evidence/prepare01/report.json)、[6项拒绝及旧v191读取](evidence/preflight01/report.json)。
- 粗任务：[准入](evidence/coarse_admission.json)、[原生监督](evidence/coarse_run01/native/report.json)、[原生输出](evidence/coarse_run01/native/stdout.txt)。从原400步发展到2000步，原历史保留；原生39.027秒，峰RSS67.37MiB。
- 粗保存失败：[首次审核](evidence/coarse_audit01/report.json)、[失败输出](evidence/coarse_audit01.stdout.txt)、[只读恢复](evidence/coarse_audit02/report.json)。实际仅0/200/400/2000步场，缺800/1200/1600；`scheduledSamplingPassed=false`与`planAcceptancePassed=false`保持。实际长窗及末50步数值指标通过，没有补造帧或重跑CFD。
- 细任务：[准备](evidence/fine_prepare01/report.json)、[独立体积/面积/初态](evidence/fine_geometry01/report.json)、[准入](evidence/fine_admission.json)、[原生监督](evidence/fine_run01/native/report.json)、[原生输出](evidence/fine_run01/native/stdout.txt)。9600单元，2400更新，原生275.282秒，峰RSS228.35MiB；实际9帧满足时间间隔合同。
- 细审核失败：[JSON bool_序列化失败输出](evidence/fine_audit01.stdout.txt)、[准确失败源码](evidence/fine_audit01/auditHybridEulerEvolution.py)、[只读恢复](evidence/fine_audit02/report.json)。首次目录没有完整report.json，不能当作已成功报告；修正序列化未重跑CFD。
- 空间对照：[报告](evidence/comparison01/report.json)、[同60区间面积积分](evidence/comparison01/common_axial_strips.csv)、[科研图](evidence/comparison01/refinement_comparison.png)。压缩压力误差+3.7569%→−0.5353%，压缩温度+2.3652%→+1.1199%；膨胀温度+4.5931%→+4.7036%，预定全部空间收益未过。两层网格不证明网格无关。
- 一次只读定位：[法向剖面报告](evidence/profile01/report.json)、[完整剖面CSV](evidence/profile01/normal_profiles.csv)。固定近壁带和熵超量均保留，不外移采样挑选通过值，不由壁邻法向速度推断壁面漏流。
- [原始字节清单](evidence_manifest.json)：150份小证据共5,456,557字节，归档时逐项SHA核对，29保护项保持。本轮本地输出22,581,900字节；无独有删除，无工程包。

12个父监督累计342.781秒，开发、读源和交付协调另计。未运行MATLAB、主UI、computer-use或新整箭CFD；原19752面、面积分母、历史黏性19层与尾缘资产保持。新增箭体热流计算面积、工程外缘资格及物理资格均为0。

清单`localOnly`列出17个仅在本地的原生rawbinary路径、大小和SHA，准确运行库及历史大场也未备份到GitHub。粗网格mesh/surface两份逐字输入复用`reusedPublic`指向的v191公开原件。公开小证据不能替代本地独有体场，不能据此声称GitHub已有完整运行库或所有历史计算结果。历史准入与尝试标记保留，不直接重跑已用任务。
