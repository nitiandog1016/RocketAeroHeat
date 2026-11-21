# B 型热敏珠 6 m/s 启动工况输入

`geometry.json` 记录名义尺寸、SI 工况和构造初场的来源。`base6.cfg` 和 `relax.cfg` 分别用于 6 m/s 扇区基流及全周松弛。`halfstep.cfg` 用于原生第一半步，第二步只切换重启和物理时间层编号。三个配置的文件字节分别与 v288、v289、v294 的实际运行一致。

`CIntegration.cpp` 是本项目在 SU2 8.6.0 的对应源文件上修改的专用观察器。它在 172,482 节点的热敏珠全周网格、单线程、第一或第二个 BDF1 物理步的指定装配阶段记录壁反力，并在附加读数时保留正式终场。作者和版权声明见文件头。本文件依 GNU LGPL 2.1 或更高版本分发，许可证全文见 [COPYING](COPYING)。固定的 SU2 源码提交为 `7048fd35f56d4bf52e586429d6746a36a10255f8`。[上游源文件](https://github.com/su2code/SU2/blob/7048fd35f56d4bf52e586429d6746a36a10255f8/SU2_CFD/src/integration/CIntegration.cpp)
