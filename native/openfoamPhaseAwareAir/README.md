# 相态有界标准干空气适配

此目录是固定 OpenFOAM 13（18870c24d21c6b982e2cdec27b2f59738cca5f90）的 GPLv3 适配库，并保留 OpenFOAM 与 Boost 许可。不包含通量、推进或回滚内核。Boost 1.83 TOMS748 用于严格60—2000K内的热内能反算；Cp采用已有 hIcoTabulated 的线性表及解析积分。生成表的具体来源、组分、适用域和证据见 `docs/validation/v206/PLAN.md`。

`rahPhaseAwareInternalEnergy.H` 检查真实目标能量和初始种子，只计算温区内的根；求解器此时尚未更新压力，不能将旧压力和新温度组合为最终物理状态。`rahPhaseDomain.H` 定义最终p/T联合域。`rahPhaseStateGuard.C` 是只读 functionObject，对每个完成时刻的全部内部和边界p/T/rho/e/U检查有限性、正密度、EOS和联合域，并核对timeIndex连续。它不修改状态，不将非法值填成零。

守卫不逐单元审核e(T)残差或完整守恒方程；这些仍由保存场独立读取和空间残差工具承担。标量查询通过不能证明foamRun中的实际分派/调度，须看相应真实任务报告。无黏Euler算例里的mu=0/Pr=.72只是原生无黏字典，不能当作工程热链使用的Lemmon2004输运。

新MATLAB完整入口为 `phaseAwareAirState.m` 和 `phaseAwareAirReferenceState.m`。仅显式新模型调用使用，旧NASA7/NIST/OTH及既有工程默认保持。恢复因子必须由区域方法提供；接口不自行给出热流、固体壁温、连续流或已认证外缘资格。60K来源下界、气相距离和局部稀薄判据均不能靠裁温/改门槛绕过。

首次隔离构建入口 `scripts/buildPhaseAwareAir.py` 需要numpy/scipy和已安装OpenFOAM/Boost；v206的WSL缺少scipy导致构建后审核中断，已通过 `scripts/auditPhaseAwareNative.py` 在Windows只读核对后继续未运行的拒绝任务。两个原生构建与已有正向查询没有重复。一次性目录保护禁止重复消耗同一预算。
