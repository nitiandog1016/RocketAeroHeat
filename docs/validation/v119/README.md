# v119实际证据

总览见[技术说明](../../SU2_ACTUAL_CONTINUATION_V119.md)，唯一续算的门槛在[运行前计划](PLAN.md)。物理未合格、新增资格面积0，数值门槛有一项未通过。

- `cases/continuation01/`：实际配置、源谱系、终态、历史、150/299表面场；完整体场仅本地。
- `processes/`：控件/适配器与结果/物性两项MATLAB实际终态和日志。
- `checks/`：原面/CSV/EOS独立审核、壁面单位高度估计、物性对照、输入拒绝、文件身份及存储清单。
- `research/`：公开DNS原文件、获取来源与输入检查。M8按拟用γ1.4/R287的头部来流核对失败；不放宽容差。
- `failed_reference_check/`：首次参考审核失败源码和记录；没有流场重算。
- `source_versions/`：实际参与的项目源码和本地固定SU2源码依据；`delivery_manifest.json`登记逐项SHA。

两项MATLAB、唯一SU2均终态，全部执行期间来源稳定。组件核对原19752面和92段截面；本段299、累计599分列。没有新主五步骤全链或原生截图验收。

本地只读复开：

```matlab
d = load('outputs/v119_validation/review01/result.mat','r');
fig = uifigure('Name','RocketAeroHeat · v119研究结果','Position',[50 50 1450 850]);
panel = uipanel(fig,'Units','normalized','Position',[0 0 1 1]);
ui = createSU2ResearchView(panel,d.r);
```

只持有GitHub文件时，可将第一行改为本目录`checks/review01/result.mat`。该文件可以显示已保存原面研究结果，不能替代未上传的大体场以重做EOS或启动引擎。源路径/回执为历史身份，不伪造迁移后的计算记录。继续计算从工作台选择仍在本地的原目录，须保留一致网格/配置并通过独立准备器；本版不能自动进入下一次续算。
