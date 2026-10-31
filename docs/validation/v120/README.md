# v120固定外域近壁对照证据

运行前门槛见[PLAN.md](PLAN.md)，结论见[技术说明](../../SU2_ACTUAL_NEAR_WALL_V120.md)。原19,752面和完整外域保留；两网格结果为未收敛研究值，不能用于宣称全面积物理资格。

- `mesh/`：2µm首位移、32层重分布的真实进程回执与小型输入；外棱柱面、外四面体和远场复用，没有运行Gmsh。
- `cases/`：原网格对照、新网格候选的配置、具名100/199表面、完整历史与回执，以及一次CoolProp运行包能力探测。
- `processes/`：MATLAB网格准备与结果组件的实际进程/源码锁定记录。
- `checks/`：独立初态迁移、原面场/CSV/EOS/热梯度/首层尺度、预设筛选及本地存储检查。
- `source_versions/`：实际使用的源码原字节；`delivery_manifest.json`登记SHA和证据身份。

本地只读复开候选研究场：

```matlab
d = load('outputs/v120_validation/review01/candidate_result.mat','r');
fig = uifigure('Name','RocketAeroHeat · v120近壁研究','Position',[50 50 1450 850]);
panel = uipanel(fig,'Units','normalized','Position',[0 0 1 1]);
ui = createSU2ResearchView(panel,d.r);
```

只持有GitHub证据时，可用本目录`checks/review01/candidate_result.mat`替换第一行路径。该文件可显示完整原面研究结果，不能代替未上传的体网格/守恒体场进行EOS重审或续算。控制结果可将`candidate`换为`control`；二者均为给定Tw300K，求解壁温保持未知。

同网格控制是累计799步，新网格是“源累计599 → 新网格局部199”，不是同网格累计799。迁移仅为初猜，没有全域积分守恒承诺。组件WebGL/探针/截面/导出测试不等于重新验收主五步骤或原生截图；本版没有增加主UI跨网格初态选择器，不能使用同网格续算入口绕过其拒绝保护。

大型体网格、初态、检查点和运行时留本地。仅合并本轮完全相同的最终/具名副本并保留全部路径；不删除独有状态、不处理v117被拒清理目录、不修改OpenCode或共享历史。
