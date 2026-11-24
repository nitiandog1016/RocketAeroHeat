# DLR RFZ 多工况表面热流公开数据核查

2026-11-24。通过 DataCite 找到 Bykerk 的论文配套软件存档，再沿存档指向的 GitHub 仓库核对实际文件。RFZ-REC 的 STEP 和八个跨声速工况的表面网格、压力、热流均已下载并回读。该套导出覆盖作者的八分之一对称计算域；读到全部公开壁面块，不等同于已经证实项目完整 CAD 的每个原面都能映射。这里记录的数值来自 DLR 的计算结果，本项目求解覆盖和实验验证面积均不因下载而增加。[作者仓库（固定提交）](https://github.com/rfz-model/RFZ_data/tree/3a174ee3d0b0626d6029e9608d862b674f3b266b) · [2023 年配套数据 DOI](https://doi.org/10.5281/zenodo.8146567)

## 可用于同一表面网格对照的八个工况

- AFMC 2024 的 RFZ 一级滑翔数据共有八个马赫数，0.75、0.85、0.95、1.05、1.15、1.25、1.35、1.45。每个工况一个 VTM 索引和 11 个 VTU 表面块（B_4 至 B_14）。实际读取全部 88 个块的坐标、单元连接、节点压力及热流。每组 196,401 个表面单元，102,294 条分块节点记录（块间接缝可能重复计点）。每块所有节点都被单元引用、压力和热流均为有限数；每块均有非零热流。八组坐标、连接、单元类型和偏移逐块相同，八份热流场各不相同。[八组结果目录](https://github.com/rfz-model/RFZ_data/tree/3a174ee3d0b0626d6029e9608d862b674f3b266b/Aerodynamics/Stage_1/AFMC_2024/Transonic)
- 热流位于 VTK 的节点数据 `PointData/heat_flux`，压力是 `PointData/pressure`。可直接对同一公开网格的节点做跨工况比较。与项目的面平均值比较时须统一节点到面元的口径，按面积对照，并核对单位、热流正方向、壁温、部件及对称面。论文的稳态计算利用八分之一对称域；复制扇区只表示对称展开，独立参考仍来自原扇区。论文还描述另一项全域非定常计算，两者须分别识别。[AFMC 原论文（表 1 与计算域）](https://elib.dlr.de/204233/) · [本地 0.75 工况入口](../../outputs/rfz_reference/author_data_20261124/files/Aerodynamics/Stage_1/AFMC_2024/Transonic/vtk_0.75/glide.surface.pval.29338.vtm)
- [RFZ-REC 三维几何](../../outputs/rfz_reference/author_data_20261124/files/CAD/stp/RFZ-REC.step) 已与固定提交的 Git 对象核对，SHA-256 为 `0442f31ef5cc0ffb0676879c17299da68245b35a3ac104d82d7d46419f2033db`。[作者 CAD 下载目录](https://github.com/rfz-model/RFZ_data/tree/3a174ee3d0b0626d6029e9608d862b674f3b266b/CAD/stp)

## 与当前 v342 首选论文的关系

EUCASS 2023 另有六组点云文件，作者脚本注明压力用 Pa、热流用 W/m²。该发布物给出 XYZ、压力、热流等数组，没有给出点与表面单元之间的连接表。它支持论文工况的点位对照。AFMC 2024 的网格和工况应单列，不能拿来给 EUCASS 点云补连接。[2023 原始数据](https://github.com/rfz-model/RFZ_data/tree/3a174ee3d0b0626d6029e9608d862b674f3b266b/Aerodynamics/Stage_1/EUCASS_2023) · [作者字段说明](https://github.com/rfz-model/RFZ_data/blob/3a174ee3d0b0626d6029e9608d862b674f3b266b/Aerodynamics/Stage_1/EUCASS_2023/visualise_data.py)

| 原目录 | 逐点检查 |
|---|---|
| `pfd10deg180`、`pfd10deg185` | 每组 70,984 条。两组 XYZ 逐项相同，热流不同；可作同位置、不同工况的点值对照。需匹配该舵面构型。 |
| `deg185`、`deg190` | 每组 72,706 条，各自长度一致。跨组 X、Y 相同，Z 逐值相反。须核实姿态与坐标约定后再配对。 |
| `deg180` | X、Y、压力、热流各 18,701 条；Z 却有 72,706 条。按原文件无法逐点配对。 |
| `pfd10deg190` | X、Y、压力、热流各 70,984 条；Z 有 72,706 条。同样暂停逐点使用。 |

另核对了 RFZ-ST2 上面级两个工况的整流罩、贮箱和喷管三块文件，均包含坐标、连接、压力和热流。仓库还列出出口边界，未计入以上核对。上面级作为单独构型留存。[作者的上面级数据目录](https://github.com/rfz-model/RFZ_data/tree/3a174ee3d0b0626d6029e9608d862b674f3b266b/Aerodynamics/Stage_2/HISST_2024)

## 其他科研档案与复核记录

- RETALT2 给出多个工况的整箭热流图和分区均值；所核对的 Zenodo 下载清单是气动力数据库，未核实同工况完整逐面热流原件。[RETALT2 论文](https://www.eucass.eu/doi/EUCASS2023-323.pdf) · [下载清单](https://zenodo.org/records/7081540)
- ESA 的 IXV 项目有整表面数据库研究；hAEdb 开发方说明其收存的 IXV 原始 CFD 数据可申请提取。此次未取得公开匿名下载包。NASA 航天飞机历史档案公开了检索指南，另有多工况风洞实测报告，但全球测热报告明确给出的区域为迎风面和左舷。[IXV 数据库论文](https://doi.org/10.1177/0954410012469493) · [hAEdb 访问说明](https://smr.ch/projects/haedb/index.html) · [NASA 存档指南](https://ntrs.nasa.gov/citations/19950013208) · [NASA 测热范围](https://ntrs.nasa.gov/citations/20030005486)
- 在项目本地 `outputs/rfz_reference/author_data_20261124/` 留存 149 份源文件（共约 247 MB）及文件清单，固定提交为 `3a174ee3d0b0626d6029e9608d862b674f3b266b`。独立检查脚本在读取前重算每个 Git blob；科学数据只读，检查未修改原数组。此前论文页和资料库简介不足以判断数据是否存在；作者的独立软件发布补足了证据。[文件核查明细](../../outputs/rfz_reference/author_data_20261124/field_audit.json) · [AFMC 原件校验清单](../../outputs/rfz_reference/author_data_20261124/transonic_all_sources.json) · [核查脚本](../../outputs/rfz_reference/author_data_20261124/check_fields.py)
