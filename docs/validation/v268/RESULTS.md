# 等温微球在一个工况下的全表面压力与热流

半径 1 mm 的等温球在缓慢流动的液体中向外放热。单工况入口为 2,608 个原三角面全部求出了压力、入壁热流和首层流体温度。上游角带的入壁热流为 −358.15 W/m²，下游为 −330.22 W/m²。正值约定为热量进入球壁。自绘箭体本次新增双变量工程覆盖 **0 m²**，仍为 0.846754 / 1.964629 m²（43.10%）。[原面审核](evidence/audit_result.json) · [自绘算例基线](../v263/RESULTS.md)

## 逐量方法与结果

| 量 | 求法 | 细网结果 |
|---|---|---|
| 壁压 | 球外 Stokes 近区用 MATLAB `bvp4c` 配点，远区用满足无穷远来流的衰减模态。在切面匹配速度与应力，由动量方程恢复壁压 | 远场 101325 Pa，球面压力扰动模态幅值 0.15 Pa |
| 入壁热流 | 近区和热浴中的有限体积离散保留对流与导热，两侧共同求解切面温度。从固壁取 Fourier 热通量 | 按球面角带积分 −4.329225896 mW；独立整壳配点解 −4.329243394 mW，相差 0.000404% |
| 工作台总热功率 | 按每个原面中心的方向取样，再用原三角形平面面积求和 | −4.318942359 mW。弦面总面积 1.25365760×10⁻⁵ m²，球面面积 1.25663706×10⁻⁵ m² |
| 温度 | 球壁给定温度与首层流体单元解各自保存 | 球壁 302 K，热浴 300 K，首层流体 301.98276–301.98411 K。求解固体壁温、外缘静温和恢复温度保留未知 |

本算例的直径雷诺数为 7.6×10⁻⁶，半径 Péclet 数为 0.19。动量模型取恒物性、无体力的蠕动流极限。流动远场延伸至无穷远，温度边界位于 8 倍半径。物理条件、网格及事前门槛见[运行计划](PLAN.md)。速度、压力扰动和阻力的公开解析参照采用 Mandre 第 9.3.7 节，并转换至静止球的参考系。[解析解与方程位置](evidence/public_reference.json)

## 穿流、收敛与独立对照

半径 2 mm 的人工交界两半球各有 7.46128255×10⁻⁸ kg/s 的质量流入和流出。细网交界净向外对流热功率 0.00721221 mW，导热 4.32201369 mW。二者之和与球壁向外功率一致，整个热域相对能量不平衡为 5.99×10⁻¹⁵。人工交界排除在固壁面积之外。[MATLAB 结果](evidence/sphere_result.json)

独立 Python 核对使用整域解析 Stokes 速度和 Chebyshev–Legendre 配点法。以相同角带平均热流相减，再除以导热尺度 `k|Tw−Tb|/a = 300 W/m²`，粗、中、细三档的最大差依次为 1.2780×10⁻⁴、3.2781×10⁻⁵、8.2956×10⁻⁶。把切面移到 3 倍半径，细网角带变化为 7.7647×10⁻⁶。按原面中心比较点值时最大差为 7.5786×10⁻⁴，包含角带取样误差。静止球有限热浴解析导热、等温零热流、超域与坏网格拒绝均通过。同源导出已逐面回读，浏览器交互本轮未复验。[独立结果](evidence/audit_result.json) · [细网角带](evidence/wall_bands/fine.csv)

以上为解析与独立数值核对。本算例没有匹配的实验测量。已有 NASA 圆柱—扩张段的压力和热流实验对照按其原测站与工况单列。[公开实验结果](../v263/TMR_RESOLVED.md)

## 运行时间与恢复

开发基线 `f6952d101311bf85b3ff313acbffd87901f86139`。最终 MATLAB 批次 28.797 s，含六个球体设置、导出及旧双间隙保存终态的兼容回读。MATLAB 等待进程退出码为 0，监督结束时无遗留子进程，采样进程组峰值 713,023,488 B。独立校核 0.891 s、峰值 52,523,008 B。两项均为各自进程耗时。前两次 MATLAB 批次为 26.015、23.969 s，首次独立校核为 1.031 s。阅读、开发及归档耗时未单独计量。[MATLAB 回执](evidence/matlab_process.json) · [独立回执](evidence/audit_process.json) · [旧入口回读](evidence/legacy_result.json)

从仓库根目录在 MATLAB 中运行以下代码，会在新目录生成六个设置的初场、直接稳态终场、细网原面表和离线工作台。

```matlab
addpath('tests')
folder = tempname;
mkdir(folder)
testHybridCreepingSphere(folder)
disp(fullfile(folder, 'fine', 'workbench', 'index.html'))
```

现有受监督运行的复核命令为 `python docs/validation/v268/check_sphere.py outputs/v268_validation/sphere03 新报告.json`，需 NumPy 与 SciPy，报告使用新文件名。独立审核同时检查该目录的真实进程回执。新机器若要从零生成可审核批次，可在仓库根目录用下列 Python 代码启动单线程 MATLAB，完成后将打印的目录传给上述审核命令。MATLAB 需在 PATH 中，监督器另需 psutil。

```python
from pathlib import Path
import shutil, sys, tempfile

root = Path.cwd()
sys.path.insert(0, str(root / "scripts"))
from runSU2NativeDebugger import observe, sha

folder = Path(tempfile.mkdtemp())
names = ["createCreepingSphereHybridCase.m", "runAeroHeatSingleCase.m",
         "runHybridSurfaceCase.m", "buildCanonicalBody.m", "tests/testHybridCreepingSphere.m",
         "buildAerothermalViewerData.m", "exportAerothermalWeb.m", "docs/validation/v268/PLAN.md",
         "web/index.html", "web/viewer.js", "web/section.js", "web/export.js", "web/viewer.css"]
locks = {str(root / name): sha(root / name) for name in names}
def matlab_path(path):
    return path.as_posix().replace("'", "''")
batch = f"cd('{matlab_path(root)}'); addpath('tests'); testHybridCreepingSphere('{matlab_path(folder)}');"
receipt = observe([shutil.which("matlab.exe"), "-wait", "-singleCompThread", "-batch", batch],
                  folder, locks, 120, 32 * 1024**2)
assert receipt["debuggerExitCode"] == 0 and not receipt["stopReason"]
print(folder)
```

本机最终场、原面 CSV、页面及早期回执的位置与 SHA-256 见[本地文件清单](evidence/local_files.json)。大场和完整页面仅存 E 盘，未在本次 GitHub 发布中备份。终场数据文件可由上述程序重建。已生成页面使用了当时工作区中并发修改的两个前端脚本，哈希随回执保留，发布仅包含本次求解器、数据导出修改和小型证据。[本机离线工作台](../../../outputs/v268_validation/sphere03/fine/workbench/index.html)
