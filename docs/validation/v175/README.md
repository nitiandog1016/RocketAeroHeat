# v175有限路线复审证据

本轮按[事前计划](PLAN.md)完成三项有限检查，没有新CFD、引擎构建、大场遍历、MATLAB或UI运行。结论为停止当前NACA发展支线，下一项返回原箭体逐变量区域证据/面积台账。[综合复审](../../reviews/TECHNICAL_ROUTE_REVIEW_V175.md) [下一决策](NEXT_DECISION.md)

| 分项 | 真实程序监督 | 结果 | 范围 |
|---|---:|---|---|
| 冷壁控制 | 0.156秒失败 + 0.234秒通过 | 32输入稳定 | v105/v108准确小回执/源码，不重算平板或相似解 |
| 原生热壁/外边界 | 0.219秒通过 | 18项检查，13输入稳定 | 准确8源与配置；不是新CFD验证 |
| 独立热合同 | 0.156秒失败 + 0.375秒通过 | 9输入稳定，独立代数/矩阵通过 | 正性、温区、凸守恒映射、弱壁身份及两边界重构差异 |

合计1.140秒，各分项在40秒内；阅读/开发/协调墙钟不包含在程序时间中。控制首次因监督目录缺失在子进程创建前拒绝，另列launch_refusal。两次实际失败分别为Python3.7编码别名和缺少math.isqrt，准确源及日志保留；不涉及数值门槛放宽。

- [控制范围](CONVERGED_CONTROL_REVIEW.md)、[原生路线](NATIVE_ROUTE_REVIEW.md)、[独立合同](THERMAL_CONTRACT_REVIEW.md)、[覆盖代码](COVERAGE_CODE_REVIEW.md)给出发现和限制。
- `evidence/`逐字保存本轮32项小运行文件及5份被审MATLAB源，包含失败、准确执行脚本、源码与最终回执；来源/长度/SHA见[evidence_manifest.json](evidence_manifest.json)。没有大场或引擎二进制，没有工程包。
- 输入回执中历史`outputs/...`路径表示原真实位置。未把本地大场说成已经上传；本轮只读的小输入/源能从回执和既有归档追溯。新模型的求解复现仍需对应运行时和本地大输入。
- 独立代数检查可使用`tests/verifyThermalStateContractsV175.py --output <全新目录> --exact-euler-source docs/validation/v175/evidence/native_review01/executed_sources/SU2_CFD/src/solvers/CEulerSolver.cpp --exact-euler-sha256 199e03e4dc258fdd040e3f74510ffd79aff13114caf96c7523c41ab5eedc3f67`，配合项目Python/Numpy；此处仅给重放入口，本轮未再次执行已过检查。

旧目标181550字节引用既有v174快照，SHA不变，没有新复制。当前目标文档与历史分离；没有改变生产数值、原失败、原−8/热流/y+或v174温域，也没有删除独有结果/共享历史。原面台账尚未实现，不计新增计算或验证面积。
