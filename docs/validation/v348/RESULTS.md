# deg190 的 Mach 标注与原始字段口径

DLR 的 EUCASS 2023 目录说明将这批结果称为 **Mach 5.35**。同一发布包的工况工作簿在 **J4** 列出 **Mach 5.38**。两处标注均予保留。工况按 RFZ-REC、滑翔阶段、作者攻角 190°、舵偏 0° 和发动机状态匹配。按这些条件筛选，工作簿只有一行匹配。论文表 4 的编号为 **3**，工作簿 `Aerodynamic Coefficients` 表第 4 行的 Case 为 **2**。同为 190° 的工作簿第 7 行有 10° 舵偏，属于另一组。[发布说明](https://github.com/rfz-model/RFZ_data/blob/3a174ee3d0b0626d6029e9608d862b674f3b266b/Aerodynamics/Stage_1/EUCASS_2023/readme) · [作者工作簿](https://github.com/rfz-model/RFZ_data/blob/3a174ee3d0b0626d6029e9608d862b674f3b266b/Aerodynamics/Stage_1/EUCASS_2023/aerodynamic_coefficients.xlsx) · [论文表 4，第 7 页](https://doi.org/10.13009/EUCASS2023-018)

| 比对项 | 论文表 4，编号 3 | 工作簿第 4 行 |
|---|---:|---:|
| 来流速度（m/s） | 1689.7 | 1689.72 |
| 来流静温（K） | 245.4 | 245.36 |
| 来流静压（Pa） | 368.1 | 368.1 |
| 作者攻角 / 舵偏（°） | 190 / 0 | 190 / 0 |
| Mach | 表内未列 | 5.38 |

工作簿 B4 的 422 s 标识轨迹取点时间。H4 给出的来流动压为 7,464.66 Pa，用 D4 速度除以 I4 声速得到 5.3811，与 J4 一致。发布说明与工作簿采用不同 Mach 标注的原因仍未知。冻结工况按论文表 4 的 SI 输入保存。[机器可读来源记录](../../../cases/rfz_rec/deg190_source_record.json) · [原冻结输入](../../../cases/rfz_rec/full_cad_input.json)

## 压力零点与热流正向

作者绘图脚本给出的单位是壁面压力 Pa、热流 W/m²。脚本把选定的原数组按色标上下限线性着色。虽然脚本从列出的速度、密度算得约 7,137.72 Pa 的动压变量，后续没有使用该变量。论文图 13 展示底板压力系数及热流，工作簿 F4、H4 记录来流压力及动压。这些来源没有明确给出 `p_data.npy` 的压力零点或 `q_data.npy` 的正法向。因而样点仍以作者原单位和符号展示，压力系数转换及有符号误差尚缺字段定义。[作者的对应绘图脚本](https://github.com/rfz-model/RFZ_data/blob/3a174ee3d0b0626d6029e9608d862b674f3b266b/Aerodynamics/Stage_1/EUCASS_2023/deg190/visualise_data.py)

相对 v347，项目新增同时有独立壁压、热流的原面为 **0 面、0 m²**。本版解决了对照时只凭 Mach 标签或 Case 编号选错工况的问题。RFZ 项目独立解和误差仍待取得。完整原壁分母为 130,044 面、740.729003 m²。

## 使用与复核

工作台的参考图题列出两处 Mach 标注，侧栏列出论文编号与工作簿行。导出的 PNG 同时写入当前变量的零点或符号核实状态。在 **计算 → 查看作者压力与热流样点 → 导出参考样点（CSV / ZIP）** 中，选择 ZIP 可同时取得原始样点 CSV 和 `source.json`。该说明包含数据 DOI、固定提交、工况来源、单位、原始文件摘要以及包内 CSV 的 SHA-256。[工作台导出示例](evidence/pressure.png) · [实窗回执](evidence/workbench.json)

2023 年发布标签和项目固定提交中的六个样点文件及对应工作簿具有相同 Git 对象。首次检查误将两版绘图脚本也按字节判等。回读原件确认两者的行结束符不同，统一换行后内容相同，原始失败日志已保存。独立解包后，72,706 行的 436,236 个场量与坐标值按原顺序验证通过。包内 CSV 与使用原入口单独导出的文件逐字节相同，单个压力值被修改的测试表在打包前被拒收。[数据发布 DOI](https://doi.org/10.5281/zenodo.8146567) · [独立解包结果](evidence/bundle.json)

在项目根目录，安装有 NumPy、openpyxl 和 pypdf 的 Python 可运行以下核查。第一条将缺少的一手资料下载到指定目录，已有资料校验后复用。第二条复核从工作台导出的 ZIP。

```powershell
python docs/validation/v348/verify_sources.py --cache outputs/rfz_source_check
python docs/validation/v348/verify_bundle.py outputs/v348_validation/ui01/reference.zip
```

基线为 v347 `9728dd96910a2aae4211239ad3d623c4e44e0f9b`。资料恢复测试复用了已核对的论文和六个 NPY，其余源文件按固定版本重新取得，用时 7.59 秒。离线复读 0.84 秒。实窗及异常数据验收 27.67 秒，进程组峰值 1.99 GiB。独立解包复核 0.69 秒。开发与文档用时未单计，流场求解用时为零。[事前计划](PLAN.md) · [运行记录](evidence/runs.json) · [下一项](NEXT_DECISION.md)
