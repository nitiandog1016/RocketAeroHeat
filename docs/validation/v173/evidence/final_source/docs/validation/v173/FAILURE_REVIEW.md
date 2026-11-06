# v173 实际候选异常终止：已记录证据与未知边界

2026-11-06。唯一 `real01` 原生运行 91.281 秒、退出码 134；父监督 108.640 秒、采样峰值 1,057,988,608 字节，85 项输入稳定。不是超时、取消或输出容量终止。原生容量记录 245,313,841 字节；完整保存 local 0 / series 1502，最后有限历史 local 56 / series 1558。末历史对应完整场和异常瞬间场没有保存，不能补造。

**本次退出的实际原因、异常类型、调用点、节点和温区方向均未知。** 不能因退出码 134、上一轮错误消息或本轮恢复事件，把最终退出写成 `domain_error`、某一边界查询或某一节点失效。

## 小范围只读核对

`scripts/reviewNACASpatialStepTermination.py` 唯一审核进程上限 15 秒，实际监督 **0.219 秒、exit 0**，脚本内部 0.078 秒。17 项输入 SHA 稳定，4 份准确源码、2 份小日志、4 条事件核对完成；读取大场、CFD、构建均为 0。准确源码与回执见 `outputs/v173_validation/termination_review01`。

已扫描 `engine_stdout.log` 的全部 **12,791 字节**，包含原始十六进制尾部检查；全部为 ASCII，无 NUL，末尾是：

```text
Simulation Run using the Single-zone Driver

```

全文件没有 `RAH_NATIVE_TERMINATION`、`RAH_PROPERTY_QUERY`、`exception_type=`、`domain_error` 或温度反演错误文本。父任务 `native_real01.log` 共 185 字节，也只有终态 JSON，没有异常消息。末尾不存在被文本显示隐藏的二进制信息。

任务目录只找到这一份文本日志；准确执行的 `runEmbeddedSU2.py` 将标准错误以 `stderr=STDOUT` 合并到同一文件，没有另行指定的 stderr 文件。这里仅确认任务目录与监督器已知日志，未检索操作系统事件日志，不能声称所有外部日志渠道均不存在。

## 准确终止处理器仍是旧版本

本轮 runtime 及准确 `source_tree.json` 已绑定，源清单 SHA256 为 `5e60fbaaccd412ced994f2831c9785ca9a7e246900101177996b2dfbaed6ada8`。`SU2_CFD.cpp` 和 `NativeTermination.hpp` 与 v171 空间壁温父包逐字相同；主对象也相同，SHA256 为 `7f5786d58fcd5b66b14428eff95eb2bb9f8c72e057392d47c538776c3c718254`。

主函数仍调用 `rah_termination::install()`，由 `std::set_terminate` 注册 C++ 终止处理器。其准确限制是：

1. `NativeTermination.hpp` 第 19 行在输出任何终止文字前检查并发/递归标志；已经进入记录时，后续进入直接 `std::_Exit(134)`。该分支不等待记录完成。
2. 首次进入先输出 BEGIN、异常信息和各栈帧，第 51 行才显式 `fflush(stderr)`，之后委托原处理器并保留 `abort()` 后备路径。它没有保证所有终止路径都完成记录和刷新。
3. 这是 C++ terminate 处理器，不是所有 Windows 异常、信号、显式进程退出或外部终止的通用捕获器。

因此，准确源码确实存在可绕过记录完成的路径；但**本次并发/递归进入、未刷新或其他退出路径均未被实测证明**。缺少消息既不能证明处理器未安装，也不能证明究竟经过哪一条路径。退出码 134 不能单独作为异常类型或 POSIX 信号诊断。

## 四条事件只证明已恢复的单元试探

`conservative_steps.csv` 4 条记录均为 `node.trial / rhoe`，`old_valid=1`、`theta=0.5`、`attempts=1`，局部历史为 46、56、56、56。记录明确给出局部 RCM 与全局编号；本复审不把两者混同，也不重复主核验器的几何/EOS审核。

记录的试探热内能依次为 164,136.497006、1,829,781.584547、1,857,901.956282、1,831,179.554023 J/kg。这些数值和各自 `reason=inverse target outside temperature range` 属于**已捕获并回溯恢复的事件**，不是最终异常的输入。不能笼统把四项都写成同一温区方向，更不能把最后一条事件当作终止位置。

准确 `ConservativeStepTrace.hpp` 在每条事件后检查写错误并调用 `fflush(trace.file)`，因此事件文件与尚未产生完整终止文字并不矛盾。四条完整事件没有赋予未保存终态任何信息。

不为补栈再运行、构建或扫描当前候选。原异常退出、缺失终态、原 −8/热流/近壁和实验门槛保持；本次只读复审不新增物理资格，也不把日志完整性审核当作计算成功。
