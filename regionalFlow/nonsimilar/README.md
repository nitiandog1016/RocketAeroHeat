# 非相似层流条件研究模块

这是从v61已测草稿固定迁出的数值模块，未接入默认生产链。沿流向保留真实物理n网格上的上游u/T剖面，求解TP-3271层流极限。所有入口/边界、温度和来源合同见两个核心完整文件。

```matlab
addpath(projectRoot,fullfile(projectRoot,'regionalFlow','nonsimilar'));
domain=struct('axisymmetric',true,'angleOfAttackDeg',0,'sideslipDeg',0, ...
    'flowRegime','continuum','layerState','laminar', ...
    'crossflowPresent',false,'separationPresent',false,'thinLayerVerified',true, ...
    'inletIdentity','填写真实入口文件/算法/条件身份', ...
    'wallTemperatureHistoryIdentity','填写完整Tw(s)历史身份');
result=solveNonsimilarLaminarStation(initial,history,gas,domain);
```

`thinLayerVerified=true`等字段是**调用者声明**，本函数不具备几何、转捩或分离自动验证能力。其输出`domainDeclarationAutomaticallyVerified=false`。数值收敛只说明单站离散求解，`physicalSolution=false`和`eligibleForSurfaceCoverage=false`恒定，失败的q和tau保持NaN。

`nonsimilarLaminarResidual`与`advanceNonsimilarLaminarStation`保留原草稿数值实现，低层接口允许显式制造源以做数值验证；公共接口不接收制造源。二者忽略O(delta/R)法向曲率修正，使用点式非保守离散，必须另查逐站保守积分能量账。原NASA7焓零点、物性范围80–2500K、给定Tw与未求壁温身份不变。

运行定向测试：

```matlab
addpath(projectRoot,fullfile(projectRoot,'tests'));
testNonsimilarLaminarModule(evidenceFolder);
```

证据保存源码SHA256、MAT和JSON。Couette闭式/制造解及输入拒绝用于检验数值合同，不是球鼻实验验证，也没有新增物理解面积。

## 冻结条件段、真实面元与离线页面

`loadNonsimilarSphericalStudy`读取`tests/data/nonsimilar_sphere/`的固定MAT及SHA清单，并用当前正式残差重放80个保存站。它不重跑长段Newton，不制造旧运行开始前的源码锁定。40°是给定局部相似入口，接受段为40.25°–60°；工况M4/20km、零攻角、球鼻半径0.1m、全程Tw=300K及NASA7物性固定。

`mapNonsimilarSphericalStudy`要求当前几何、实际自由流、物性、完整Tw身份一致；检查完整平面三角形的极角界，按真实质心线性插值，不外推。`sphericalTrianglePolarBounds`计算边内部极小值和轴线穿过情形，避免顶点检查漏掉越界弦段。`nonsimilarStudyDigest`用于检测缓存/站表变更，不是安全签名。

`exportNonsimilarSphericalStudy`单独生成条件研究页面，复用现有WebGL/探针/截面/CSV合同，不改变main默认流程或生产覆盖率。恢复温度和求解壁温保持NaN。面元常值导致曲线阶跃/锯齿，不能解读为横流物理分布。

```matlab
% 完整定向示例：构造匹配球鼻网格、重放、映射、输出独立页面和证据。
addpath(projectRoot,fullfile(projectRoot,'tests'));
testNonsimilarSphericalMapping(fullfile(projectRoot,'outputs','v61_validation','mapping'));
```

测试的152面仅为该几何的结果，接口不硬编码面数。正式工程包中的`example_nonsimilar/web/index.html`可直接离线打开；重放证据在`validation/v61/mapping/`。定向测试18项拒绝，公共站接口29项拒绝；逐站能量/动量独立审计见`docs/V61_ENERGY_AUDIT.md`。当前状态是**已测试的条件研究模块**，未获得生产区域解或实验精度资格。
