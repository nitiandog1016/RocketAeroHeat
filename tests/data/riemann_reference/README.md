# 独立 Riemann 参考的小型原文

此目录用于独立数值校核，不是外部 CFD 场，也不替代 MATLAB 求解器。

- `raw/clawpack_exact_euler.py`：2026-10-24 获取的 Clawpack `riemann_book/exact_solvers/euler.py` 原始代码，17157字节，SHA256 `2a926d59d4935807b93c425aa617bf64d64858e65231fd556b6ae6202d0ec2ef`。来源：[上游代码](https://github.com/clawpack/riemann_book/blob/master/exact_solvers/euler.py)。代码采用 BSD 3-Clause，完整声明随附 `raw/LICENSE-CODE`；作者为 David I. Ketcheson、Randall J. LeVeque 和 Mauricio J. del Razo。
- `raw/probe_exact_riemann.py`：本项目2026-10-24保存的独立 Brent/二分采样脚本原文，7640字节，SHA256 `402ffeee97aed006c6784fa6b46d3276a3ab428897d77d07921308d9708c067a`。它用上述公开压力曲线作额外核对，不是上游原文件。

`tests/riemannReservoirReference.py` 先核对这两个文件 SHA，再只提取公开压力曲线和本地 `exact` 函数。不会运行原脚本的文件写出、绘图或交互代码。原文字节由 `.gitattributes` 的 `tests/data/**/raw/** -text -eol` 保持。BSD 声明来自上游 `LICENSE-CODE`，2026-10-28通过已连接 GitHub读取，Git blob为 `9fb3b0b009fbd584e5afe744b657c68a29ba2c2a`。

复现 v91 小域独立核对：先在 MATLAB 中调用 `testRocketCellInterfaceLookahead('outputs/v91_validation/new_unit')`，然后运行 `python tests/auditRocketCellInterfaceLookahead.py outputs/v91_validation/new_unit`。输入制造场与两单元夹具也已随库保存，无需下载大型工程包。该核对属于解析/数值参考验证，不是实验热流验证。
