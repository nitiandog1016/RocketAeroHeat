# v193证据索引

结论见[实现报告](../../HYBRID_ROCKET_EULER_MESH_V193.md)，下一项见[NEXT_DECISION](NEXT_DECISION.md)。本轮生成真实完整Euler网格和采样接口，网格候选未获接受，没有新流场、热面积或物理资格。

- 事前计划：[原PLC准备](PLAN.md)、[唯一refine](REPAIR_PLAN.md)；[外缘与物性合同](EDGE_CONTRACT.md)。
- 原网格：[准备报告](evidence/mesh01/report.json)、[原生PLC审核](evidence/mesh01/plc/interface.json)、[原生进程](evidence/mesh01/plc/process/process.json)、[父监督](evidence/mesh01_process.json)。原生exit0、候选因6个低质单元拒绝，二者不能混称成功。
- 独立核对：[面积/体积/质量/约束](evidence/geometry01/report.json)、[低质单元坐标](evidence/geometry01/bad_cells.csv)、[原面完整采样CSV](evidence/mesh01/wall_sampling_geometry.csv)。
- 修复失败：[拒绝报告](evidence/repair01/report.json)、[原生进程](evidence/repair01/native/process/process.json)、[父监督](evidence/repair01_process.json)。
- 只读定位：[8个原壁连接改变](evidence/repair_review01/report.json)、[新增边界连接](evidence/repair_review01/changed_boundary_triangles.csv)。未经接受的修复场最低质量与10低质数量见[资源/来源评估](evidence/assessment01/report.json)。
- 采样接口：[首次语法失败](evidence/sampling01.stdout.txt)、[准确失败源](evidence/sampling01_failed_source/hybridEdgeSampling.py)、[第二次schema失败](evidence/sampling02.stdout.txt)、[最终检查](evidence/sampling03/report.json)、[480控制面候选](evidence/sampling03/control_candidate.json)。480面为原生已存控制场；19752面只是原箭体源网格几何，不是真实箭体新流场。
- [准确官方/本地源片段](evidence/assessment01/source_excerpts.json)、[正投影定位图](evidence/figure01/rocket_mesh_screen.png)。首次三维图也保留，不作为热云图。
- [证据清单](evidence_manifest.json)：52份小证据共5,245,111字节；8个大网格/原生数组路径仅本地，总本轮输出25,922,114字节。29保护项保持，独有删除0，无工程包。

8父监督累计19.515秒，原生TetGen1.031/.813秒含于相应父监督；补充只读质量与正投影绘图的工具墙钟约1.49/1.69秒另列，开发/读源/交付不计入程序预算。运行库、历史大场与清单`localOnly`文件未上传GitHub，GitHub小证据不等于完整历史备份。不得覆盖历史目录或重跑本轮已关闭候选。
