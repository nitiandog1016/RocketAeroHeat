# v163 原生物性异常定位与高精度保存

## 完整域结果

实际任务real02复用v160局部149，初场、网格和配置分别与v161控制逐字相等。原19752面、335726节点、705490单元，M5、20°、给定Tw300K。运行包为rah-ideal-air-observed-v1，数值基线为phase-v2，新增C++终止记录与17位ASCII。

原生69.282秒后以3221225501退出，父监督86.172秒、峰值进程树RSS约2.24GB，时间/内存限制均未触发。终止记录包含std::domain_error及“inverse target outside temperature range”，随后调用原终止处理器。完整帧局部0、最后有限历史16分别对应系列151次、167次更新，异常瞬间场未保存。[实际终态](validation/v163/evidence/outputs/v163_validation/real02/engine_process.json)

调用路径为OpenMP迭代、流动预处理、CNSSolver::SetPrimitive_Variables、CNSVariable::SetPrimVar、CSU2IdealAir::SetTDState_rhoe、rah_air::temperature_from。固定二进制反汇编确认返回地址0x1404804a3对应首次当前态查询，原生合法性检查和Solution_Old回退尚未执行。已确定的阻断为物性越域异常绕过原生回退。异常节点、rho/e/k数值、越过温区的方向和候选状态越域的具体动力学原因仍待后续改动记录。[调用栈复核](validation/v163/evidence/outputs/v163_validation/native_review.json) [调用点复核](validation/v163/evidence/outputs/v163_validation/callsite_review.json)

## 保存精度与物理资格

| 检查 | 结果 |
|---|---|
| C++相邻9温度逐位往返 | 16位5/9，17位9/9 |
| 最长标量文本含分隔符 | 26字节，容量预留32字节 |
| 小型0/2/3全字段相对差最大 | 4.35e−16 |
| 真实新旧首帧全字段相对差 | 5.72e−16 |
| 真实首帧独立阶段能量误差 | 2.00e−13 |
| 真实首帧EOS误差 | 4.43e−16 |
| 真实首帧热梯度误差 | 7.36e−15，通过原2e−10 |
| 原残差−8门槛 | 未通过 |
| 首帧y+99 | 2.07121，未通过 |
| 新物理资格 | 0 |

首帧温区88.687—1678.033K，壁流体270.675—382.151K，给定Tw300K。完整首帧接口已通过，进程整体仍为异常终态。v161原始低精度首帧失败保持，旧内存值未恢复。[首帧独立审核](validation/v163/evidence/outputs/v163_validation/real02/audit0/report.json)

## 构建与失败记录

独立原生夹具7.250秒完成编译，普通异常、OpenMP异常、非std异常和实际文本精度四项通过。引擎范围为三个翻译单元，CEulerSolver对象逐字恢复phase-v2，其余295继承对象逐字保持。累计编译29.751秒，唯一链接1.953秒，共享源恢复通过。

main前两次编译因Windows的byte及宏名称冲突失败，第三次采用受限Windows头并在SU2/Eigen头之后包含观察器，实际通过。小型首次准备因错误历史快照路径拒绝，恢复Git准确源后唯一4更新0.859秒exit0。真实首次准备因误填49输出间隔触发容量预检，恢复原50后启动唯一真实任务。四次失败均保留准确源和回执。[事前计划](validation/v163/PLAN.md) [恢复记录](validation/v163/RECOVERY_PLAN.md)

## 下一项

主五步骤只读验收79.359秒exit0、1733源输入稳定，新包方法、8原面、4段周向导出、未知求解Tw及旧phase-v2读取通过。独立原面压力/热流/温度CSV最大相对差3.30e−15，压力曲线与原面一致。首版独立曲线审核误用热流列，准确失败源保留，按实际字段身份修正后只读通过。异常整箭任务继续由正常结果入口拒绝。此次采用已有WebGL实现，无新原生鼠标或截图验收。[工作台检查](validation/v163/evidence/outputs/v163_validation/main01/report.json) [导出独立复核](validation/v163/evidence/outputs/v163_validation/export_review.json)

限定变比热模型，将首次物性拒绝接入已有旧态恢复，完整重算相同k下的热状态。旧态仍越域时记录并终止。下一项先做准确函数小检查，再独立数值合同和有限完整域验证。[实施决策](validation/v163/NEXT_DECISION.md)

小源码、完整MATLAB替换文件、目标和172项小证据同步公开GitHub，证据原文件合计10,471,290字节。源码散列表超归档自设大小限制后，通过逐字校验恢复归档，原失败脚本保留。大场及运行包仅本地，无工程包、无独有文件删除，本轮清理0字节。总目标继续进行。
