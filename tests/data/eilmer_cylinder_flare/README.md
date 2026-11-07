# CUBRC Run 14 圆柱—外扩段公开参考

本目录只保存固定GDTK v5.0.0源码树中13个小型原文件，共105,928字节。版本为
`f53f4609a0331d48efee69a4e4f3c3598378cc03`；下载URL、字节数、Git blob SHA1和SHA256见
`source_manifest.json`。原文件逐字保留，不执行其中的Lua/D源码。

来源为[官方算例](https://github.com/gdtk-uq/gdtk/tree/f53f4609a0331d48efee69a4e4f3c3598378cc03/examples/eilmer/2D/cylinder-flare/sg)。
代码许可与作者说明原样保存在`raw/COPYING`、`raw/gpl.txt`、`raw/contributors.txt`；上游声明GPL-3.0或以后版本。
两份实验表保留CUBRC测量身份，不把算例数值解称作实验。原论文只引用
[MacLean与Holden，AIAA 2004-0529](https://chimeracfd.com/professional/public/aiaa2004_0529.pdf)，没有复制整篇论文。

离线准备命令：

```text
python scripts/prepareEilmerCylinderFlareReference.py --output outputs/v177_validation/reference01
```

程序拒绝覆盖既存输出；生成SI热流/压力测站CSV、`conditions.json`、原始输入锁和执行回执。
测站X严格用原表4.004英寸；原表未明确Btu变体，项目明确采用Btu_IT并保留原数值。
54条活跃热流和24条活跃压力全部保留，另保留注释行。P30/P32超出官方示例220 mm尾端，标为未映射，不能裁到末面；HT56/P31缺失值保持空白。

此候选是具名轴对称层流理想N2控制，不是已运行/已验证的计算，也不授真实20度三维RANS、底面、突起或全面积热环境资格。
固定源码实际使用CEA黏度，低于200 K返回端点常值；导热因Prandtl选项优先而使用`cp*mu/0.72`，并非CEA导热拟合。
实验Tv=2467 K、p/T/rho不完全相容、输运范围及逐点不确定度仍须明确处理；准备程序不猜补，也不改变默认求解器。
