# 一次计算能交付多少原面结果

常用工程入口与混合入口现在都返回本次壁压和热流同时有解的原面数、面积、剩余面积，以及已算面和整表面的有符号热功率。五步工作台在接受本次计算后显示这组数。离线三维页面把它放在结果侧栏开头，逐变量色标、同面探针和各自的覆盖率仍可分别查看。[汇总函数](../../../summarizeComputedWall.m) · [本机微球工作台](../../../outputs/v269_validation/check02/workbench/index.html)

| 输入 | 同面已算 / 全部原面 | 同面面积 / 全面积 | 已算面热功率 | 整表面热功率 |
|---|---:|---:|---:|---:|
| 自绘模型，同一工况存档 | 9,725 / 19,752 | 0.846754 / 1.964629 m²（43.10%） | +121.667 kW | 未知 |
| v268 低速微球，同一工况存档 | 2,608 / 2,608 | 1.253658×10⁻⁵ m²（100%） | −4.318942 mW | −4.318942 mW |
| 标准球，常用工程入口实际执行 | 81 / 168 | 0.0582645 / 0.121107 m²（48.11%） | +669.369 W | 未知 |

自绘模型仍有 1.117875 m² 缺压力、热流或两项。已有热流的面积为 1.307709 m²，面积大于双量交集，因此这些面的积分可以报告。微球按显示的原三角面面积积分，与 v268 的原面结果一致。风洞对照按源壁环带面积积分，复得 24.642781 kW。三种工况各有独立分母。本轮新增自绘模型求解面积和物理验证面积均为 **0 m²**。[逐工况结果](evidence/result.json) · [独立读表](evidence/export_check.json)

## 调用与验收

在 `[current,fs,report]=runAeroHeatSingleCase(mesh,cfg)` 中，`report.jointComputedFaceCount` 是原面交集的面数，`jointComputedAreaM2` 和 `jointMissingAreaM2` 是对应的已算与缺项面积。只有全部原面都同时算出两量时，`fullSurfaceComputed` 才为真。`calculatedHeatPowerW` 为已算面功率，`wholeSurfaceHeatPowerW` 需要每个原面都有计算热流。真实零热流、相抵得到的零功率以及负值都予以保留。两量分布错开、给定压力或全部缺热的情况都已经测试。[接口测试](../../../tests/testComputedWallCoverage.m)

在 MATLAB 仓库根目录运行 `addpath('tests'); testComputedWallCoverage()` 可复验公共入口。要回读本机存档并生成页面，运行 `addpath('docs/validation/v269'); out=tempname; mkdir(out); check_coverage(out)`。测试核对几何、各场数值、状态和逐面算法，并对比原三角面积与源环带面积。静态页面的摘要绑定和位置已通过文件回读核对，本轮未实际操作浏览器。[存档复核程序](check_coverage.m) · [计划及输入](PLAN.md)

基线为 v268 的 `e12e1ba0`。首轮 40.984 s，最终 MATLAB 批次 19.079 s，均含启动和测试；最终采样进程组峰值约 883 MiB。两次监督回执及结果完整，最终受监督的 22 个输入在运行期间一致，退出时无遗留子进程。独立读表与小文件归档 0.25 s，编辑耗时未单列。本机页面沿用了工作区当时的前端脚本，清单记录页面各文件 SHA-256；大型存档和页面留在 E 盘。GitHub 包含本次源码和小型证据，原有并发脚本独立保留。[最终回执](evidence/matlab_process.json) · [首次回执](evidence/first_process.json) · [本机页面校验值](evidence/export_check.json)
