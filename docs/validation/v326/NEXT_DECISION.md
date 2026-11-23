# 下一项只检验双支管试件的守恒扩散通量修正

采用同一 Short Test 65 混合网、原冷态、逐节点定温和 v326 的内部质量通量。v326 在第 8 次更新的最低值流体节点上，测得非正交扩散修正大于沿边扩散与对流项；冻结旧场的逐边幅值限制足以改变这一行的残差方向。下一项据此固定一个候选系数 0.5，检验它在真实推进中能否维持全域温差范围，并向原收敛门槛推进。[该点的原始分解](RESULTS.md)

参考 OpenFOAM Foundation 13、提交 `18870c24d21c6b982e2cdec27b2f59738cca5f90` 的 `limitedSnGrad`（GPL-3.0-or-later）。该实现按每个面的两点项和原非正交修正确定限制系数，参数 0.5 将修正幅值限制在两点项幅值内。本机核过的 `.C`、`.H` 的 SHA-256 分别为 `58a1ceb12433fce37b180066c0610cc75e9dc26ddfcae86dacb09d3c1e371569`、`e1a2044cbb8a47822a78fd0355e1f9b02776769e358578b01c3ee0164bd7c4e9`。OpenFOAM 用单元中心和面，当前 SU2 用节点和对偶面，两者的两点系数及切向贡献定义有别。采用的是逐面限制思想，具体适配只在本试件的诊断驱动中验证。[上游实现](https://github.com/OpenFOAM/OpenFOAM-13/blob/18870c24d21c6b982e2cdec27b2f59738cca5f90/src/finiteVolume/finiteVolume/snGradSchemes/limitedSnGrad/limitedSnGrad.C) · [上游定义](https://github.com/OpenFOAM/OpenFOAM-13/blob/18870c24d21c6b982e2cdec27b2f59738cca5f90/src/finiteVolume/finiteVolume/snGradSchemes/limitedSnGrad/limitedSnGrad.H)

先用常量、正交线性场与斜网仿射场核对通量，分别报告原修正、受限修正、对角处理、解析法向通量和两侧等量反向性。说明斜网限幅引起的一致性误差及网格细化趋势；事前登记准入精度与失败条件。保留切向扩散、原物性、已登记的标量范围与残差要求。任何制造场一致性或内部守恒门槛失败，即停止该候选并记录失败位置。

通过通量核对后，从原冷态启动一个独立 200 次稳态阶段，4 线程、300 s、4 GiB、450 MiB。保留周期、终帧和首次越界更新，核对全域上下界、六分量原始残差与同状态的边界通量。若再次越界，结束这项对照。达标后沿用既有的连续热稳定、同阶段壁反力及热量账，再查空间、近壁与外域误差，按五区返回原面。实验对照仅使用测量定义相符的外露球面量热和七站读数。下一阶段启动前另登记输入和计算预算。
