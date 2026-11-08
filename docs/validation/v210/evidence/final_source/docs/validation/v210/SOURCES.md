# v210来源、分发和证据范围

上游Eilmer固定提交`f53f4609a0331d48efee69a4e4f3c3598378cc03`，原生GasModel/factory/Brent/complex接口沿用v209逐字审核；GPLv3或以后，许可证见[前轮归档](../v209/evidence/GDTK_COPYING)。factory覆盖文件保留上游代码并只增加具名模型分支，补丁明确改动；运行库不随本轮提交。

物理数据和Lemmon来源沿用v206/v207，完整原始数据已在`data/thermo/phase_aware_dry_air_v1.json`，本轮不重复拷贝或声称新实验验证。只验证新原生接口/实复数导数/联合域拒绝，不能把有限状态通过外推成全域物理准确性。

准确执行代码在各任务`executed_entry.py`/`executed_prepare.py`/`executed_runner.py`，失败结果保持；final_source为交付时源，不替代历史运行源。prepared02仅替换helper导入，气体/factory/数据不变；query_inputs02按原生行号复用，audit02只修导数参考分段，没有重算CFD、裁T或放宽容差。归档清单明确本地二进制和未上传资产。
