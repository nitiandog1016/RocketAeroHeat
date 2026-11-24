# 在新安装中计算 Short Test 65 双支管银球

这条入口支持 Short 第 65 次风洞试验的名义银球、两支管及湿端盖，Re = 882，14D 远场。两种时间步到达扰动后的同一时刻 2.565681 ms。壁温取自试验，压力、给定壁温下的入壁热流及首层流体温度由程序计算。[试件与测量定义](../v321/SOURCE.md)

## 在 MATLAB 工作台新算

1. 在可写的 E 盘目录检出项目，准备 64 位 Windows、带桌面界面的 MATLAB 和能从 MATLAB 启动的 `python`。本机使用 Python 3.7.7、NumPy 1.21.6、SciPy 1.7.3 和 psutil 5.9.8。首次安装还会用到 Windows 的 `tar.exe`。为项目预留 12 GiB 磁盘空间，供工具链、网格、稳态与物理时间场使用。
2. 从[试件包](../../../cases/short_sphere/supported65.zip)解出 `supported_m.stl`。在项目根目录运行 `rocketAeroHeatWorkbench`，加载 STL，长度比例设为 **1**，确认朝向。网格页选择 **双支管银球 · Short Test 65** 及 **14D**，生成网格并确认。完成原面对应后，边界页显示已登记的风洞来流和七站壁温。
3. 在边界页选择 **两个半步** 或 **一个整步**。在计算页按 **从扰动初场重新计算**。首次启动会生成体网、获取工具链、编译程序、从冷态推进稳态基流并核对其热量收支。随后生成扰动，并分别计算一个整步、两个半步。所选步长决定结果页显示哪个终点。计算目录写在侧栏日志里。
4. 查看结果后可关闭窗口。再次加载同一试件和工况，在计算页按 **重开已完成的本地计算**，填写任务目录或留空后在文件选择器中选取。选的是含 `completed.json` 的整个任务目录。面元探针、七个测站及逐面 CSV 使用所选步长的保存结果。**读取所选步长的归档** 则打开随版本发布的对照结果。

工作台和命令行复用 `outputs/short65_inputs/14D` 中已核对的体网、专用程序与稳态基流。每次新算在独立任务目录保存初始扰动、三步体场、独立释放反力和逐面结果包。保留任务和缓存可分别用于窗口重开与从原场重算。20D 对照的现有归档和运行依赖见[外域对照](../v335/RESULTS.md)。

## 工具链和命令行

在 MATLAB 命令行运行 `!python -c "import numpy, scipy, psutil; print(numpy.__version__, scipy.__version__, psutil.__version__)"` 可核对实际解释器。要在已有的 64 位 Python 3.7 环境复现本机包版本，使用 `python -m pip install numpy==1.21.6 scipy==1.7.3 psutil==5.9.8`。

首次运行会将 Gmsh 4.15.2 和固定 SU2 8.6.0、MinGW、Meson、Eigen、MEL 安装在本项目的 `tools/runtime` 下。SU2 两个静态库采用两路编译。若已下载固定的上游包，可先运行以下命令安装内核，省去重复下载。目录内使用 `compiler.7z.exe`、`su2.zip`、`meson.zip`、`eigen.zip` 和 `mel.zip` 五个文件名，程序会逐个校验内容。[上游地址](../../../scripts/prepareSU2ThermalBuild.py) · [固定大小和 SHA-256](../../../scripts/installShort65Kernel.py) · [Gmsh 来源](../../../scripts/installGmshRuntime.py)

```powershell
python scripts/installShort65Kernel.py --downloads E:/RocketAeroHeat/downloads
```

从命令行直接新算可使用下面第一条命令。第二条从第一次计算的实际前后时间层和独立释放壁行重算逐面输出，输出目录须是新路径。

```powershell
python scripts/shortSupportTime65.py outputs/short65_first
python scripts/shortSupportTime65.py outputs/short65_reread --completed outputs/short65_first/run
```

已存在完整 v332/v333 运行输入的老安装仍可使用历史基流。首次准备的运行依赖限于该银球工况，支管几何和给定温界的适用范围见[实验对照](../v333/RESULTS.md)。
