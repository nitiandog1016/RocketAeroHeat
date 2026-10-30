# 耦合激波层数值参考

`axis_reference.json`保存v100实际M5/30km上游、给定球形激波、两种物性7至19节点摘要和19节点轴线剖面，以及缩放/真实原面元的未通过项。源完整MAT和实际进程的SHA保存在文件内。

该文件是约120KB的数值参考，不能代替完整场或赋予真实箭体/热流资格。完整复现输入沿用`tests/data/shock_layer_transport/parabolic_shock_layer_evidence.mat`中的原几何、自由流和鼻部种子；不读取其旧壁面场。

复现及当前限制见`docs/validation/v100/README.md`。本次没有实验验证，正式新增热流面积为0。
