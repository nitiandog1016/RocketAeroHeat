# NIST气体原生注册适配

本目录只注册OpenFOAM Foundation 13已有模板组合，热模型、插值积分、状态方程与求逆仍使用其原生实现。固定源码提交为18870c24d21c6b982e2cdec27b2f59738cca5f90，官方源码和许可证见https://github.com/OpenFOAM/OpenFOAM-13 （GPL-3.0）。不修改共享安装，不在Git中分发编译库；构建输出留项目E盘outputs。

调用方式、误差、已知失败和未取得的物理资格见docs/HYBRID_NIST_AIR_V202.md。不能将模板注册成功当作任意输入能量/温度反算种子的稳健性证明。
