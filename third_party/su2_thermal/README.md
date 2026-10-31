# SU2 8.6.0 原生轴向壁温适配

本目录是 RocketAeroHeat 的受控输入适配，采用 LGPL-2.1-or-later，许可证见 `../su2/COPYING`。上游为 SU2 提交 `7048fd35f56d4bf52e586429d6746a36a10255f8`。这不是未经修改的官方二进制。

`RAHAxialWallTemperature.hpp` 读取严格递增的 SI X—温度表；`SetCustomBoundary.cpp.inc` 在上游已有 CustomBoundaryTemperature 数组初始化时写入每个壁节点温度。仅支持单区、固定网格、串行 MPI/OpenMP、SI、MGLEVEL=0、一个名为 wall 的自定义等温壁。未指定环境变量时保留上游标量行为。空间输入运行由监督器显式设置 `RAH_WALL_TEMPERATURE_FILE`，任意漏项都不得静默回退。

构建脚本另外给 PreprocessPythonInterface 加遗漏标记保护，并修复三个上游 FSI/混合面文件的相对头文件路径。没有修改 NS 方程、壁面热流公式、线性求解器或收敛门槛。准确补丁、工具来源和失败记录在 `docs/validation/v115/`。

复建步骤：

1. `python scripts/prepareSU2ThermalBuild.py` 下载固定源码和编译工具；各归档身份记录在安装目录，交付证据另存下载清单。
2. 用支持 7z 的工具将 `tools/runtime/su2-thermal-build/compiler.7z.exe` 解压至 `tools/runtime/su2-thermal-build/toolchain/`，得到 `w64devkit/bin/g++.exe`；不要修改系统 PATH。本机实测用已安装 Bandizip 的 bz 解压。
3. `python scripts/buildSU2ThermalRuntime.py --attempt NN`；NN 选尚无记录的数字。本机 Ninja 路径在脚本中明确记录，迁移机器时应显式配置并重新记录来源。每次两线程、最多720秒，保留失败终态与已编译对象。
4. 在独立新目录安装已通过的构建：`python scripts/installSU2ThermalRuntime.py --attempt NN`。脚本拒绝覆写既有运行时；若本机已安装，直接复用其完整 SHA 清单，不为了复建删除当前依赖。
5. 运行解析器、常量极限、空间热梯度、拒绝和 MATLAB/UI 测试。仅编译成功不能宣称数值或物理通过。

MATLAB 当前的变壁温模式选择 `8.6.0-rah-axial-wall-v1`；等温、绝热默认仍用官方 `8.6.0`。运行时仅存本机，不随 GitHub 上传；源码与构建配方公开。没有 MPI、动态网格、多区、多重网格或非定常输入适配资格。
