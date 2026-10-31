# SU2 有界标量守恒 Jacobian 受控修正

固定上游：SU2 v8.6.0 / 7048fd35f56d4bf52e586429d6746a36a10255f8。上游许可证为 LGPL-2.1-or-later，源码文件原版权声明保留。此目录提供可审查补丁，不发布未标记的“官方”二进制。

SST 隐式未知量是 rho*k 与 rho*omega。有界质量散度修正残差为 -phi*div(massFlux)，在冻结流动状态下，对守恒未知量的导数为 -div(massFlux)/rho。本补丁修正着色边两端、归约边和边界函数四个对角表达式；Conservative=false 分支保持原值。残差、生产耗散、扩散、原下限和松弛规则不改。它不保证全 SST 更新正值或收敛；SA 和其他标量模型不在当前验证范围。

复现入口：
- scripts/verifySU2BoundedChainRule.py：提取准确原生边界函数，编译差分核对和独立纯对流反例。
- scripts/buildSU2BoundedRuntime.py：核对固定源码与 Eigen/MEL，复用现有构建目录但强制重编 CFD 的全部对象；新运行包独立安装，最后还原临时源码。仅保留三个既有编译包含路径修复，不带旧热边界/冻结能量钩子。需要既有固定工具链与 Meson 构建图，来源见 scripts/prepareSU2ThermalBuild.py。
- scripts/reviewSU2BoundedBuild.py：独立核对实际构建树差异与安装字节。
- scripts/checkSU2BoundedPackage.py：官方/新包原对流兼容与新包有界接口，各同初态两步。
- scripts/prepareSU2SSTBoundedStudy.py：只有前置检查通过才准备一次原整箭100更新比较。原参考终态不作新初态。

各运行入口拒绝覆盖已有回执；失败必须保留并决定下一步，不能直接重跑。结果与限制见 docs/validation/v131/ 和对应版本说明。运行时/大场只存本地，源码和小型证据同步 GitHub。
