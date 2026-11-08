# 温变输运的原生注册适配

固定 OpenFOAM 13 提交 `18870c24d21c6b982e2cdec27b2f59738cca5f90`，采用 GPLv3；完整许可复用 [OPENFOAM_COPYING](../openfoamPhaseAwareAir/OPENFOAM_COPYING)。有界反算使用 Boost 1.83 TOMS748，许可见 [BOOST_LICENSE_1_0.txt](../openfoamPhaseAwareAir/BOOST_LICENSE_1_0.txt)。上游类型及本项目适配来源分别保留。

本目录只注册 `icoTabulatedTransport + rahPhaseAwareInternalEnergy + hIcoTabulatedThermo + perfectGas`，没有新流动/时间推进内核。构建入口将既有 `rahPhaseAwareInternalEnergy.H`、`rahPhaseDomain.H` 和两项许可复制到独立构建目录；原始头文件位于相邻 `openfoamPhaseAwareAir`。量热、气相/低密度域、60—2000 K 范围、参考焓和严禁外推规则均不变。

`scripts/phaseAwareViscousAir.py` 用已有 Lemmon 2004 稀薄极限输运公式在3882节点生成成熟线性mu/kappa表，Pr随温度变化。输运存在不代表该流场符合连续介质假设，也不自动获得湍流/分离热资格。

实际库位于 E盘 `outputs/v207_validation/native01/library/librahPhaseAwareViscousAir.so`，SHA256为 `75f6204a37b68772750c3b5b626e7c91dd01689459774c4d65e58f0a5224cc8e`。该独立注册库与v206最终状态守卫库同时加载；守卫本身不重复注册在本库。实际构建源及许可留在v207证据中。

仅已验证32单元层流Couette的两项8步控制。`scripts/runViscousThermalControl.py` 的CLI已关闭，不能直接重复运行历史任务。复现或新区域构建/CFD必须另登记独立预算和目录。新读取器只支持已测层流Fourier、fixedValue壁温、准确气体字典；原生wallHeatFlux取负得到项目正向入壁热流，负热流不能裁零。
