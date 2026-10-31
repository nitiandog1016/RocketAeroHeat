# v112 证据索引

结论及下一起点见[专题文档](../../SU2_MULTIGRID_CHECKPOINT_UI_V112.md)。新增合格面积0，总目标未完成。

- `PLAN.md`：原数值／UI计划。UI第一次真实任务发现回执竞争后，仅为修复该缺陷新增一次同类5秒任务；多重网格只跑一次，没有延长。
- `engine/mg01/`：一层粗网格实际600循环、176.391秒exit0，数值筛选失败。
- `engine/ui_failed/`：Windows原子替换异常的原文及未闭合回执。实际子进程不再存在，但原终态／退出码缺失，不能补写成成功或正常超时。
- `engine/ui_timeout/`：修复后实际5.078秒超时；原100步检查点及末监测103分列。与M5候选无关。
- `failed_supervisor_observation.json`：旧监督器／子进程不存在的独立后验，原回执未改写。
- `checks/`、`processes/`：三项MATLAB实际任务，1失败／2成功；两次Windows文件锁测试的原文及解码修复；旧成功快照复开回归。
- `independent_audit.json`：原面映射、同一初态、残差／时间筛选、双面积、全部UI面场与48段截面CSV独立数值核对。
- 三张`main_*.png`：实际M4超时研究工作台，含原生窗口切页操作；不是M5物理验证截图。
- `source_versions/`：实际使用的相关代码版本，包括失败时的监督器与UI测试；源码锁稳定性按各任务分别核对。
- `storage_plan.json`、`storage_receipt.json`、`storage_audit.json`：11组逐字重复合并，全部路径／SHA／文件身份通过，独有状态保留。
- `delivery_manifest.json`：小型证据和准确来源绑定。体网格、大重启动仍只在本地，GitHub不构成这些文件的备份；复算须按保留来源写入新目录。

具名检查点及超时研究是显式操作，不能自动宣称结果收敛或区域获得物理资格。当前无通用任意STL／突起体网格或真实19,752面NS完整解。
