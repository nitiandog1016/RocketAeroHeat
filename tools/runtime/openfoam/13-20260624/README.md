# 固定OpenFOAM13原生运行入口

运行库安装于E盘承载的 `RocketAeroHeat-Ubuntu-24.04` 内 `/opt/openfoam13`，Debian `openfoam13=20260624`，GPL-3.0-or-later。Git只保存本说明和 `runtime.json`，不分发二进制。官方来源：https://openfoam.org/download/13-ubuntu/ 。

`runtime.json`绑定实际二进制/库/关键源码和三个准确上游blob。全部发行源的Git逐项比较未执行，不能由build前缀推断全部一致。调用通过具名进程加载发行包环境，不修改用户shell启动文件。

NASA7与现有工程物性相容合同见 `scripts/openFoamAirThermo.py`；原生RR常数舍入差通过具名EOS参数显式处理。此安装包并不自动意味着热流/三维箭体求解已验证。v198 LTS小控制已因温域警告关闭，后续按目标与NEXT_DECISION另行准入。
