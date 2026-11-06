# v172 证据索引

本轮具名同场冷重启、自适应反馈和保存帧接口已测试；唯一实际候选物性反解异常，不采纳/不续算/不改默认，物理与实验资格均为0。

- [事前计划](PLAN.md)：准确绑定运行预算及发展筛选，保留原文。
- [前置合同评审](CONTINUATION_REVIEW.md)：属于实现前建议，实际配置/源步/次数以PLAN和真实回执为准。
- [准确自适应源评审](ADAPTIVE_REVIEW.md)、[原列分辨诊断](COLUMN_REVIEW.md)：分别为源语义与发展场诊断，不是网格资格。
- [实际异常复审](FAILURE_REVIEW.md)：残缺栈、真实反馈和未定位的物性调用者。
- [本轮完整说明](../../NACA_ADAPTIVE_DEVELOPMENT_V172.md)、[下一决定](NEXT_DECISION.md)、[存储回执](storage_receipt.json)。

小证据由`evidence_manifest.json`列出原路径、归档路径与SHA256，重复内容引用已有v171归档，不复制全部旧证据。`finalSources`标明最终实现；失败时源码快照独立保留。运行时仅记准确清单/二进制身份，GitHub不包含大场/完整运行包。

主要本地报告位于`outputs/v172_validation/`：`saved_frame_regression01`、`adaptive_contract01`、`adaptive_review01`、`column_fixture01`、`column_real01`、`small01/small02`、`audit_small01`、`verify_small01/02`、`real01`、`audit_real01`、`development_real01`、`verify_real01`、`failure_review01`、`plot_history01`。父级保存准备/原生/审核监督回执及失败日志。

关键身份：原生exit134；normalRunCompleted/interfacePassed=false；仅完整local0/series1502的savedFrameInterfacePassed=true；最后history45/series1547无完整场。验证器passed=true是复核执行通过，不是求解成功。终态壁流体温度及对给定Tw偏差、y+及两完整帧p/q变化待验，首帧值不能代替。

科学图`plot_history01/adaptive_history.png`已生成并实际查看，图示真实CFL/残差/有符号功率，注明冷重启与异常身份。无新MATLAB/主UI验收；通用空间图脚本来源标注只属于实现修改，本轮未重新渲染该脚本。
