# v107验证资料

- `delivery_manifest.json`：4项MATLAB终态、6项SU2终态、源码版本、证据SHA、13.84MB本地输出预算及已交付v106共享网格引用。
- `independent_audit.json`：SciPy重新读取原面MAT，逐列比较三份CSV，独立算面积和两个有符号热功率。未收敛、物理资格false。
- `processes/`：实际MATLAB回执与原日志，包含exportapp失败。
- `engines/`：实际运行/取消/关窗/输入改变/配置错误/超时回执和输入。`case03`唯一输入不稳定项是故意改动的测试文件，按失败记录保留。
- `checks/`：主工作台、最终重开、CSV、实际截图观察及负向合同检查。
- `main_workbench_final.png`：最终实际MATLAB主窗口；`main_workbench_before_layout.png`为修正工具栏/页脚之前的观察，不能冒充最终版本。
- `replay_seed/`：小型本程序原面研究快照、输入与实际回执，供无SU2运行的UI重开；不是外部CFD数据，绝不赋予收敛或物理资格。
- `source_versions/`：准确源码原字节，包括失败和未最终交付的阶段版本；使用manifest对应的版本，不把当前文件套到历史运行。

在项目根目录用MATLAB复现最终研究UI（输出目录须尚未存在）：

```matlab
addpath('tests');
testSU2ResearchReplay('outputs/v107_replay_reproduction', ...
    'docs/validation/v107/replay_seed/mapped_research.mat',false);
```

在已安装固定版SU2及配置Python后，`testEmbeddedSU2Workbench`可验证实际完整流程，包含一次短迭代和两次取消；这是接口测试，不是收敛物理基准。需要研究实际箭体收敛时遵守目标中的预算和最新起点，不重复该短迭代作为物理进展。

两次汇总工具的路径错误只修复报告，不重新执行引擎；失败信息见 `python_initial_report_failure.json` 和对应源版本。旧工程状态回归本轮实际使用336面经典模型，不重复19752面长链测试。
