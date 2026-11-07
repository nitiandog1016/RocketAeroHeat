# Eilmer CEA空气回调来源与许可

`scripts/eilmerAirThermo.py`将固定GDTk/Eilmer提交的空气CEA9曲线、系数混合、能量/熵参考及Sutherland输运适配为逐面Python回调。本适配沿用GPL-3.0-or-later；本目录保留上游COPYING与完整GPL文本。此说明不声称第三方作者审核本项目结果。

- 上游：[gdtk-uq/gdtk](https://github.com/gdtk-uq/gdtk/tree/f53f4609a0331d48efee69a4e4f3c3598378cc03)
- Copyright (C) 2015-2016 Peter A. Jacobs, Rowan J. Gollan, and other contributors；CEA曲线源署名Rowan G. and Peter J.
- 准确来源：`src/gas/thermo/cea_thermo_curves.d`、`therm_perf_gas_mix.d`、`src/gas/species-database/air.lua`、`physical_constants.d`、`composite_gas.d`。
- 本项目改动：Python向量化、固定生成气体身份检查、200—6000 K表外拒绝、显式区别原生混合区cp与dh/dT；未改变求解器或原生热物性。原生查询及源码SHA见`docs/validation/v191/evidence/gas01/report.json`。
- 本次不分发Eilmer二进制；准确本地运行包来源和构建信息沿用v187。
