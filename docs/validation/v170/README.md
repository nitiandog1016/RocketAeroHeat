# v170验证证据导航

- [事前计划及一次准备纠错](PLAN.md)
- [20重叠段原文/中心身份复核](SEGMENT_REVIEW.md)
- [外缘定义、合成验证及拒绝范围](EDGE_REVIEW.md)
- [实际网格和准确运行包空间Tw路线审查](MESH_REVIEW.md)
- [本轮实现与限制](../../NACA_SEGMENTS_MESH_V170.md)
- [下一执行起点](NEXT_DECISION.md)
- [准确文件清单](evidence_manifest.json)

小证据按原始字节保存63份、247538字节。清单复用相同SHA的存档引用，不复制同一末版本；失败与当时执行源均保留。原图和新段界独立放在 `tests/data/naca_tn3344/segments_v170`，既有原PDF不重复提交。

最新通过：`evidence/mesh02/report.json`、`evidence/mesh_verify02/report.json`、`evidence/segments02/report.json`、`evidence/edge02/report.json`（相同文件若已被另一条引用收录，以清单archive字段为准）。第一次独立网格验证也通过，第二次是针对审查指出的全参考棱柱/精确连接/法向首邻三项补充；没有新生成网格。

实际 `mesh.su2` 与 `mesh_arrays.npz` 仅本地保留，清单登记准确SHA和字节数，不把GitHub源码同步说成大网格已备份。无独有结果删除、无新工程包。CFD、引擎构建、真实箭体、MATLAB、UI本轮均0；总目标继续，物理/实验资格仍0。
