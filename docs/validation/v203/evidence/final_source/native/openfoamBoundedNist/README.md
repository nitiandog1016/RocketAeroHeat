# 具名NIST有界内能反算

`rahBoundedNistInternalEnergy`继承OpenFOAM原生能量映射，只把The的Newton反算替换为Boost.Math TOMS748。物性表、es/hf/熵、流动推进及原生边界算法保持。输入es超出表对应范围直接拒绝，不能裁剪或外推；Pressure/T/seed合同与Python端必须一致。

固定Foundation OpenFOAM13提交18870c24d21c6b982e2cdec27b2f59738cca5f90（GPL-3.0），Boost1.83包1.83.0-2.1ubuntu3.2（Boost Software License 1.0）。上游算法源码：https://github.com/boostorg/math/blob/boost-1.83.0/include/boost/math/tools/toms748_solve.hpp 。实际包源码/许可证/SHA和所有失败见docs/validation/v203；本目录不重写TOMS748或CFD内核。

适配库为独立新类型，不覆盖v202原sensibleInternalEnergy。原生标量和一步分派已通过，真实候选触发100000Pa压力合同并关闭，尚无真实正向发展场/热资格。不得通过修改本头文件阈值绕过当前准入；接续见docs/validation/v203/NEXT_DECISION.md。
