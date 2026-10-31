# v114证据目录

当前记录经典锥柱扩张体、原生绝热热边界和统一科研工作台的定向实现与测试，不是完整气动热资格。

- `PLAN.md`：原运行前计划，已绑定单元任务，不回写历史。
- `PLAN_UPDATE.md`：实际引擎启动前明确采用具名检查点和3000迭代上限，因为绝热热负荷恒为零，不能用它唯一识别保存迭代。
- `SU2_THERMAL_SOURCE.json`：固定8.6.0官方源码的热边界/输出依据。零热流是边界输入；`MARKER_PYTHON_CUSTOM`本身不会读温度表，当前运行包不含pysu2。
- `source_pressure.png`及`pressure_digitization_overlay.png`：原图8和8个手工拾取点叠图。原图1给出名义测站位置，数据见`tests/data/schaefer_ferguson_1961/pressure_fig8.json`；没有把低Re低压测点、整理线或理论线计作独立精度依据。
- `source_versions/`：实际任务使用的源版本；不覆盖旧版本。
- `processes/`：4项MATLAB任务终态exit0及对应源锁、日志。`mesh/`和`engine/`保存唯一实际网格/引擎终态、小型输入与历史；SU2未收敛。
- `checks/independent_audit.json`：3648原面、32段CSV、具名2999检查点及8个压力点的分别审核。1个名义测站未通过位置匹配，所有点均未获实验资格。
- `audit_initial_failure.json`与准确源码：首次独立审核因鼻锥点差1.936px超过1.5px失败；保留原点/阈值，最终报告显式不通过，未改成实验通过。
- `checks/area01/report.json`：经典、真实箭体及整流罩的历史清单总面积勘误；新求解器配置/网格SHA与原输入相同。旧清单原字节保留，不能继续使用其中错误wallArea。
- `checks/review01/`：复用相同场验证最终结果组件、绝热身份及旧等温300K兼容。`final_pressure_component.jpg`对应最终组件；`main_adiabatic.jpg`是修正NaN文案前的主窗口，不能冒称最终整窗。
- `checks/storage_receipt.json`：4组逐字重复文件改成硬链接，原路径/SHA保留，执行窗口磁盘空闲实增18,505,728字节。没有删除1000/2000/2999独有状态；共享文件不能原位改写。
- `delivery_manifest.json`：交付源版本、任务、证据SHA、文档及本轮独立文件长度。该长度不是全outputs物理磁盘占用。大型体网格/重启动仍仅本地，无工程包；源PDF可从固定NASA链接按SHA重新取得。

空间变化Tw目前仅完成输入合同及拒绝保护，没有实际进入引擎。热流图11的成组Tw/T0、Taw/T0、St数据仍待数字化与条件匹配。压力点数字化本身也不是实验验证。

小型MAT/CSV可用于核对与复开，完整数值复现仍需按计划重新生成体网格并使用固定SU2/Gmsh运行时。当前测试脚本中使用outputs的只读审计要求保留的本地依赖，不能仅凭GitHub源码宣称大型历史场已备份。
