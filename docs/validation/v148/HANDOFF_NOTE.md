# 接续核对说明

接续时先检查实际进程和终态。AUSM+-up2原生进程已完成，其热状态审核也已由先前监督流程完成，exit0。

随后一次尝试调用同名审核监督器，在创建子进程之前被 `assert not folder.exists()` 拒绝。没有启动第二次审核或数值求解，没有覆盖已有回执。原命令为 `outputs/v148_validation/run_bounded.py outputs/v148_validation/up2_audit_request.json`，异常位于该监督脚本第11行。后续直接使用并检查已有终态。

本轮整箭数值运行仍只有控制/候选各一次；不把本次准备拒绝记作第三次数值计算。
