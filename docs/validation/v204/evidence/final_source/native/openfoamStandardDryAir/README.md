# OpenFOAM标准干空气联合域映射

独立注册`rahStandardDryAirInternalEnergy`，使用v204标准干空气Cp表、相对hs/s参考及工程低密度域。沿用Boost1.83 TOMS748；未改流动内核/原NIST映射。实际状态必须满足100—2000K与`0<p<=min(1e6,1e4*(T/100)^3)`Pa。

先拒绝真能量超域，再在固定p允许温区内求根，所有返回路径检查联合域/正密度。最多8个相邻浮点数定位区间下界，不修改输入或允许域。能量求值上限64次，包含3次表/联合域端点、60次内部及1次末态；反算成功不等于流场收敛、已验证外缘或壁热。详见`docs/validation/v204/NATIVE_DOMAIN_PLAN.md`。

Boost使用Boost Software License1.0（随附）。OpenFOAM注册接口继续遵守其GPLv3；项目旧物性和默认模型保持。
