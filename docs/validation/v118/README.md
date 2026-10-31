# v118证据：原19,752面/精确30km来流

本目录随源码提交，状态为**已测试的完整域启动研究，物理未合格**。概述和下一起点见[技术说明](../../ACTUAL_ROCKET_FULL_DOMAIN_V118.md)，运行前计划见[PLAN.md](PLAN.md)，预算失误和一次重试见[BUDGET_REVISION.md](BUDGET_REVISION.md)。

- `inputs/`：原表面与精确fs、18层Gmsh设置和准备报告；不是外部CFD数据。
- `mesh/`：唯一Gmsh实际终态、日志和几何报告。完整混合体网格仅本地。
- `cases/engine/`：首次取消任务的配置、来源、监测和取消记录，没有壁面/体场。
- `cases/budget_run/`：300步实际任务及150/299表面场；每个重启动体场SHA登记在终态，本目录不复制约41MB的体场。
- `processes/`：3项MATLAB真实终态/全源码锁/日志；`source_versions/`保存实际参与模块的准确原文。
- `checks/`：准备、原面/CSV/EOS独立审核、双检查点演变、控件与生产WebGL组件结果；来源/存储/未改文件核对。
- `delivery_manifest.json`：实际参与源码、任务和逐文件SHA。历史初场、体网格、运行时及大型输出并未因提交GitHub而备份。

本地现成结果可在项目根目录用MATLAB复开，不运行求解器：

```matlab
d = load('outputs/v118_validation/review01/result.mat','r');
fig = uifigure('Name','RocketAeroHeat · v118研究结果','Position',[50 50 1450 850]);
panel = uipanel(fig,'Units','normalized','Position',[0 0 1 1]);
ui = createSU2ResearchView(panel,d.r);
```

若只有GitHub检出，可将第一行改为读取本目录`checks/review01/result.mat`。此MAT保留原面值及历史审核身份，可以只读显示；原本地路径没有随仓库迁移，不能从缺失的大体场再做原生EOS审核。

重新构造时，读取`inputs/mesh_input.mat`中的V/F和optionsJSON，以新输出目录调用`buildSU2HybridVolume`；读取`inputs/surface.mat`中的fs/flight，设置300步/150间隔后调用`writeSU2RocketCase`，使用`preparation.json`的显式初猜厚度调用`writeSU2HybridInitialState`。需安装项目固定的Gmsh/SU2运行时；新任务重新记录来源与预算，不能覆盖历史目录、直接套用旧SHA或声称新网格逐字一致。不得为了文档复现重复计算本轮已经通过的内容。

所有压力/热流仍为未收敛研究值，给定Tw300K不算求解壁温；新增合格面积0。没有新的主五步骤全链验收或原生截图。三个截面各92段，组件WebGL和数值导出已经实际核对。
