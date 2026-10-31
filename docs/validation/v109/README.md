# v109证据与复现

对应[版本说明](../../SU2_ROCKET_STARTUP_SECTIONS_V109.md)。`PLAN.md`先于计算确定了一个假设、原门槛和预算。`delivery_manifest.json`绑定源码、全部终态、失败及证据SHA。

- `engine/`：唯一三维长任务的原输入、初场说明/独立审核、终态、stdout、原壁面VTK与历史选行。`history_selected_rows.csv`是原文件**选行**，不是完整历史；全历史SHA和行数在manifest。stdout在超时前有未刷出的缓冲，不能由其末行猜测检查点迭代。
- `independent_audit.json`：启动改善、全部残差、检查点3750/终态3921分别标识，EOS和原面映射核对；15组实际CSV共600行逐值/状态审核，最终滑块导出也核对。
- `checks/checkpoint01/`：三维超时研究检查点的2448原面MAT，未授予物理资格。它没有伪造主工作台已确认工况输入，不能当主程序的同工况重开包。
- `checks/sections02/`：主工作台复用v107真实快照的15组截面、实际CSV/JSON、独立Node交线审核。截图和这些测试用M4、α20°、β5°旧研究场；**不是本次M5、α20°、β0°三维数值结果**。
- `checks/final_ui01/`和`main_section_final.png`：最终主窗口说明与原生滑块操作，X从0.7m移到约1.0021065m，交线和曲线同步；导出的当前截面与网页回执一致。`section_before_label.png`为改正旧侧栏说明之前的实际点击证据，不冒充最终版本。
- `processes/`：3项MATLAB终态，2成功/1失败。所有源码/显式输入保持；失败经过定位后只修复和重开UI，没有重算流场。
- `source_versions/`、`failures.json`：准确阶段版本与失败说明。两次Node审核错误只修审核器，保留版本，没有改物理字段或生产截面算法。

在项目根目录用MATLAB重开已交付的小型真实快照，不需运行SU2：

```matlab
addpath('tests');
testSU2Sections('outputs/v109_section_reproduction',false);
```

随后运行`node tests/testSU2SectionView.js outputs/v109_section_reproduction`。此入口验证主五步骤重开、实际WebGL、15组截面、点击原面中心与导出；它不验证数值物理收敛。用户仍从`main`同一工作台的研究结果页使用新截面功能，无另开的永久专题页。

本次三维研究的源几何为`outputs/v108_validation/normal01/normal_volume.mat`，可按v108源码/输入重建。`scripts/prepareSU2RocketInitialState.py <源目录> <新outputs目录>`生成平滑初值；`tests/auditSU2RocketInitialState.py <新目录>`做独立状态检查，均不启动求解。需要不同研究时先写计划，不能重新执行相同720秒任务当新进展。

`outputs/v109_validation/rocket_seed01/initial_flow.csv`、`restart_flow.csv`和完整历史仍只在本地；保留最后完整检查点和源几何依赖，不因GitHub同步删除。原面小MAT、源码和小型证据公开，但这不等于大三维状态已备份。本轮未生成工程ZIP或工程副本；未删除独有状态或共享历史对象。
