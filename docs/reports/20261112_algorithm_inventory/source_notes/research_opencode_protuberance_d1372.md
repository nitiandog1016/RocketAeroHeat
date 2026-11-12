# NASA-TN-D-1372 原文核查（任务三）

日期：2026-10-30  
范围：核对 NASA-TN-D-1372 的正式书目身份、可获取的原文全文、试验条件、表格与图号；只转录原文明确写出的条件，不把图读值或 OCR 猜测写成测量表。  
结论性质：**来源核查，不是实验验证，不是求解器对比。** 本文不写“实验验证通过”。

## 1. 结论摘要

- NTRS 正式记录和官方全文 OCR 文本已取得；全文文本 1,601,061 字节，SHA-256 为 `665C95E4852B4B8845CDD6FCA95AF9804123D1E50DCD364BB14D979E09F6F218`。
- 原文为 Burbank、Newlander、Collins 的 **NASA TN D-1372**，1962 年 12 月，Langley Research Center；不是部分二手引文误写的“D-1732”。
- 原文摘要称测试了 **20 种构型**，来流 Mach 为 **2.65、3.51、4.44**，单位长度 Reynolds 数约 **1.3×10⁶～4.7×10⁶ /ft**，湍流边界层厚度变化约 **10:1**。
- 三张热流表按边界层厚度分组：Table III 为 **δ=0.70 in**（p.47 起），Table IV 为 **δ=1.50 in**（p.115 起），Table V 为 **δ=6.00 in**（p.241 起）。Table VI 是无限长圆柱的理论热流系数（p.344 起），不是测量表。
- 表 III/IV/V 的逐点列包含测点坐标、周向角、`Te/Tt`、`Tw`、`h` 和 `h/h0` 等信息；但当前只有 OCR 文本，表格列序在 OCR 中被拆散，**不能可靠还原逐点对应关系**。因此本次只交付条件登记表，不转录逐点热流、壁温或压力数值。
- 原文的 `h` 由瞬态壁温法反算，含“蒙皮等温、忽略横向导热、忽略向模型内部传热、忽略辐射”等假设；绝对 `h` 还存在因驻点温度探针位置造成的相位差（原文称约 15%～20%），`h/h0` 比值才是原文强调的可比量。
- 未运行任何 CFD、热流、壁温、v100、Q41 或大型 NS 求解；新增正式合格面积 **0 m²**，未新增物理解面积。

## 2. 实际检索与获取记录

| 项目 | 地址/命令 | 实际结果 | 说明 |
|---|---|---|---|
| NTRS 记录页 | `https://ntrs.nasa.gov/citations/19630001117` | 浏览器式 `webfetch` 返回 403 | 无 `User-Agent` 时被拒；不代表记录不存在 |
| NTRS API | `https://ntrs.nasa.gov/api/citations/19630001117` | HTTP 200，JSON 3,069 字节 | 取得题名、作者、报告号、摘要、下载链接 |
| 全文文本 | `https://ntrs.nasa.gov/api/citations/19630001117/downloads/19630001117.txt` | HTTP 200，1,601,061 字节 | 官方全文 OCR；SHA-256 见上 |
| PDF HEAD | `https://ntrs.nasa.gov/api/citations/19630001117/downloads/19630001117.pdf` | HTTP 200，`Content-Length: 226,898,018` | 扫描 PDF 很大；HEAD 只能证明对象存在 |
| PDF 全文 GET | 同一地址，浏览器式 `User-Agent` | 本次任务窗口内未完成 | 传输很慢；未把未完成的局部文件当作原文证据 |
| 纯文本 API | 同一 NTRS API | 403 | 未使用 |
| 二手摘要 | 未使用 | — | 题名、条件、表号均以 NTRS 全文和原文页为准 |

补充说明：

- NTRS 页面在没有浏览器式 `User-Agent` 时返回 403；API 和全文下载在带 `User-Agent` 时可用。这是本次实际遇到的访问条件，不是原文缺失。
- 全文 OCR 的页脚/页眉印刷页码可读，表 III、IV、V、VI 的起始页分别由 p.47、p.115、p.241、p.344 的页脚直接核对；Table I 起始页 p.33、Figure 1 起始页 p.348 也由页脚核对。
- 本文件中的页码均为**印刷页码**，不是 PDF 物理页序号。由于 PDF 全量下载未完成，未做“印刷页码 ↔ PDF 物理页”的逐页对照；引用时以表号/图号为主，页码为辅。

## 3. 书目身份

| 字段 | 值 | 来源 |
|---|---|---|
| 题名 | Heat-Transfer And Pressure Measurements on a Flat-Plate Surface and Heat-Transfer Measurements on Attached Protuberances in a Supersonic Turbulent Boundary Layer at Mach Numbers Of 2.65, 3.51, and 4.44 | NTRS API；原文封面 |
| 作者 | Paige B. Burbank, Robert A. Newlander, Ida K. Collins | 原文封面；NTRS 作者字段 |
| 单位 | NASA Langley Research Center, Hampton, VA | 原文封面；NTRS |
| 报告号 | NASA TN D-1372 | 原文封面；NTRS `otherReportNumbers` |
| 日期 | 1962 年 12 月 | 原文封面；NTRS publication date |
| NTRS ID | 19630001117 | NTRS |
| 旧登录号 | 63N10991 | NTRS `legacyMeta.accessionNumber` |
| 分发 | PUBLIC | NTRS |
| 主题词 | altitude test; turbulent boundary layer; flat plate; pressure distribution; aerodynamic heat transfer | NTRS keywords |

**编号核对：** 原文封面、页眉和 NTRS 记录均为 **D-1372**。部分二手引文写成 D-1732，属于编号错误；本文件不采用该错误编号。

## 4. 原文结构（印刷页码）

| 部分 | 页码 | 内容 |
|---|---|---|
| Summary | pp.1–2 | 20 种构型、Mach/Reynolds/δ 范围、主要结论 |
| Introduction | pp.2–3 | 外露管路/突起背景、两类模型、试验范围 |
| Symbols | pp.3–4 | `h`、`h0`、`h/h0`、`Tw`、`Te`、`Tt`、`δ` 等定义 |
| Apparatus, Models, and Instrumentation | pp.5–10 | 风洞、平板、支撑、模型、热电偶、压力耙 |
| Accuracy | pp.10–11 | 压力、Mach、温度、热流系数精度 |
| Data Reduction | pp.11–12 | 压力与热流换算、瞬态壁温法假设 |
| Discussion of Results | pp.12–30 | 流场显示、压力、热流分布与关联 |
| Index to Tabular Data | pp.31–32 | 表 I–VI 的页码索引 |
| Tables I–VI | pp.33–347 | 测点位置、压力、三组热流、理论表 |
| Figures 1–71 | pp.348–418 | 试验装置、模型照片、油流、压力、热流图 |

页码口径：

- Summary、Introduction、Symbols、Apparatus、Accuracy、Data Reduction、Discussion 的页码由 OCR 中的印刷页码标记核对。
- Table I 起始 p.33、Table III 起始 p.47、Table IV 起始 p.115、Table V 起始 p.241、Table VI 起始 p.344、Figure 1 起始 p.348，均由页脚数字直接核对。
- Table I 的结束页、Table II 的起始页、Table VI 的结束页、Figure 71 的结束页由相邻部分的起始页和索引推定；已在表中标为范围，不把推定写成逐页核对。

## 5. 原文明确写出的试验条件

### 5.1 总体条件

| 项目 | 原文值 | 页码/位置 |
|---|---|---|
| 构型数 | 20 种构型 | Summary, p.1 |
| 来流 Mach | 2.65、3.51、4.44 | Summary/Introduction, pp.1–3 |
| 单位长度 Reynolds 数 | 约 1.3×10⁶～4.7×10⁶ /ft | Introduction, p.3 |
| 边界层厚度变化 | 约 10:1 | Summary, p.1；Introduction, p.3 |
| 风洞试验段 Mach 范围 | 2.3～4.65 | Apparatus, p.5 |
| 平板尺寸 | 4 ft × 5 ft 仪器化平板 | Apparatus, p.5 |
| 热电偶 | 84 支 | Apparatus/Table I, p.5 |
| 静压孔 | 10 个 | Apparatus/Table I, p.5 |

### 5.2 边界层厚度

| 位置 | 名义 δ | 页码 |
|---|---:|---|
| 前位置，x=40.47 in | 0.7 in | Results, p.16 |
| 后位置，x=116.47 in | 1.5 in | Results, p.16 |
| 风洞侧壁（最大边界层） | 6 in | Results, p.16 |

原文说明这些名义值与参考 11 的预测式在“良好一致”范围内；δ 随 M 和 R 变化。Table III/IV/V 分别按 0.70、1.50、6.00 in 分组。

### 5.3 三张热流表的条件

| 表 | δ | 起始页 | 构型节数 | 内容 |
|---|---:|---:|---:|---|
| Table III | 0.70 in | p.47 | 9 | 逐点热流与 `h/h0` |
| Table IV | 1.50 in | p.115 | 17 | 逐点热流与 `h/h0` |
| Table V | 6.00 in | p.241 | 12 | 逐点热流与 `h/h0` |
| Table VI | — | p.344 | — | 无限长圆柱理论热流系数，不是测量 |

每个表节内部按 Mach/Reynolds 组分组，组头给出 `M` 和 `R`；不同构型的 Mach/Reynolds 组合并不完全相同。由于 OCR 把表格列拆散，本文件没有逐组转录这些组头，也没有把组头与逐点行强行配对。

### 5.4 精度

原文 p.11 给出热流系数精度：

| `h` 范围（Btu/ft²-sec-°R） | 精度 |
|---|---:|
| > 0.0150 | ±10% |
| 0.0010～0.0150 | ±15% |
| < 0.0010 | ±20% |

温度测量精度为 **±2°F**（热电偶丝校准、冷端和记录器分辨率，p.11）。

### 5.5 热流定义与单位

- `h` 的单位为 **Btu/ft²-sec-°R**（Symbols, p.4；表 III 表注）。
- `h/h0` 定义为同一热电偶位置“有突起时的局部热流系数”除以“无突起平板的局部热流系数”（pp.17–18）。
- 表 III/IV/V 以 `h/h0` 形式给出结果；逐点绝对 `h` 也列在表中。

### 5.6 壁温与驻点温度

- 原文没有给出一个适用于所有测点的固定壁温。`Tw` 是每个测点的实测列；瞬态法用 `Tw(t)` 反算 `h`（Data Reduction, pp.11–12；表列头含 `Tw`）。
- 驻点温度 `Tt` 由探针测量；第一阶段探针位于支撑结构上方，后续阶段改到试验表面同侧。原文称这造成平板绝对热流系数的相位差：0.70-in 边界层约 **15%**，1.50-in 边界层约 **20%**；`h/h0` 比值仍有效（pp.17–18）。
- 因此，使用这张数据集时不能把不同阶段的绝对 `h` 直接混在一起；应优先用同一阶段的 `h/h0`，或在复算时严格匹配 `Tt` 探针位置与 `Tw` 历程。

### 5.7 热流换算假设

原文 p.12 明确列出瞬态壁温法假设：

- 蒙皮厚度方向温度均匀；
- 横向导热可忽略；
- 向模型内部的传热可忽略；
- 辐射热损失可忽略。

这些假设是原文数据不确定度的一部分，不应在项目中被当作“直接测量的壁面热流”。

## 6. 表清单与构型节

### 6.1 表清单

| 表 | 页码 | 内容 |
|---|---|---|
| Table I | p.33 起 | 热电偶和静压孔位置、突起相对位置；84 支热电偶、10 个静压孔 |
| Table II | p.35 起（索引） | 平板与边界层耙压力；II-1 局部 Mach 数，II-2 压力系数 |
| Table III | pp.47–114 | δ=0.70 in 的逐点热流与 `h/h0` |
| Table IV | pp.115–240 | δ=1.50 in 的逐点热流与 `h/h0` |
| Table V | pp.241–343 | δ=6.00 in 的逐点热流与 `h/h0` |
| Table VI | pp.344–347 | 无限长圆柱理论热流系数 |

Table II 的起始页来自索引（p.31），结束页由 Table III 起始 p.47 反推；没有逐页核对 Table II 的最后一页。

### 6.2 Table III（δ=0.70 in，p.47 起）

| 节 | 构型 |
|---:|---|
| 1 | Flat plate alone |
| 2 | 1-inch by 2-inch rectangular stiffener |
| 3 | 2.8-inch-diameter instrumented cylinder |
| 4 | 2.8-inch-diameter instrumented cylinder and cap plate |
| 5 | 1.4-inch-diameter instrumented cylinder |
| 6 | 1.4-inch-diameter instrumented cylinder and dummy cylinder 6.4 diameters upstream, offset 26.5° |
| 7 | 2.8-inch-diameter cylinder swept forward 45° |
| 8 | 2.8-inch-diameter cylinder swept back 45° |
| 9 | 30° cableway |

### 6.3 Table IV（δ=1.50 in，p.115 起）

| 节 | 构型 |
|---:|---|
| 1 | Flat plate alone |
| 2 | 1-inch by 2-inch rectangular stiffener |
| 3 | 2.8-inch-diameter instrumented cylinder |
| 4 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 3.2 diameters upstream |
| 5 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 6.4 diameters upstream |
| 6 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 3.2 diameters upstream, offset 26.5° |
| 7 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 3.2 diameters upstream, offset 45° |
| 8 | 2.8-inch-diameter instrumented cylinder and cap plate |
| 9 | 1.4-inch-diameter instrumented cylinder |
| 10 | 1.4-inch-diameter instrumented cylinder and dummy cylinder 6.4 diameters upstream |
| 11 | 1.4-inch-diameter instrumented cylinder and dummy cylinder 6.4 diameters upstream, offset 26.5° |
| 12 | 2.8-inch-diameter cylinder swept forward 45° |
| 13 | 2.8-inch-diameter cylinder swept back 45° |
| 14 | 30° cableway |
| 15 | 26° fairing |
| 16 | 45° fairing |
| 17 | Lox line |

### 6.4 Table V（δ=6.00 in，p.241 起）

| 节 | 构型 |
|---:|---|
| 1 | Flat plate alone |
| 2 | 2-inch by 4-inch rectangular stiffener |
| 3 | 2-inch by 4-inch stiffener with 1/4 round fairing |
| 4 | 2-inch by 4-inch stiffener with 30° wedge |
| 5 | 1-inch by 2-inch rectangular stiffener |
| 6 | 2.8-inch-diameter instrumented cylinder |
| 7 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 3.2 diameters upstream |
| 8 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 6.4 diameters upstream |
| 9 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 3.2 diameters upstream, offset 26.5° |
| 10 | 2.8-inch-diameter instrumented cylinder and dummy cylinder 3.2 diameters upstream, offset 45° |
| 11 | 2.8-inch-diameter cylinder swept forward 45° |
| 12 | 2.8-inch-diameter cylinder swept back 45° |

**20 与 9+17+12 的关系：** 摘要的“20 种构型”指不同的几何/硬件构型；三张表按三个 δ 重复测试，所以表节总数为 9+17+12=38，不等于 38 种独立构型。不能把 38 个表节当成 38 个新实验。

## 7. 图清单（关键图号）

| 图号 | 内容 | 页码/说明 |
|---|---|---|
| Figure 1 | 试验平板、热电偶/静压孔位置、典型模型相对位置 | p.348 起 |
| Figures 2–5 | 层压平板、翼型支撑、密封、试验段安装 | pp.348–350 附近 |
| Figures 6–7 | 可换平板面板、热电偶与突起位置 | 图题可见，页码未逐页核对 |
| Figure 8 | 四种整流/凸条模型及热电偶位置 | 图题；模型照片见 Figs.15–18 |
| Figure 9 | 2.8 in 直径正圆柱模型 | 图题 |
| Figure 10 | 2.8 in 圆柱加端板 | 图题 |
| Figure 11 | 1.4 in 直径圆柱 | 图题 |
| Figure 12 | 2.8 in 圆柱 45° 后掠 | 图题 |
| Figure 13 | 30° cableway 模型 | 图题 |
| Figure 14 | 外部管路/整流罩模型 | 图题 |
| Figures 15–31 | 典型突起照片 | 图题 |
| Figures 32–34 | 小/大边界层耙 | 图题 |
| Figure 35 | 整流/凸条模型油流 | p.380 |
| Figure 36 | 2.8 in 圆柱油流，M=2.65/3.51/4.44 | p.380 附近 |
| Figures 37–40 | 串列/偏置圆柱油流 | pp.384–388 |
| Figures 41–42 | 后掠/前掠圆柱油流 | p.388 |
| Figure 43 | 典型边界层剖面 | p.388 |
| Figure 44 | 平板压力分布 | p.391 |
| Figure 45 | 整流/凸条上游压力 | p.392 |
| Figures 46–48 | 圆柱/后掠压力与 Cp | 图题 |
| Figures 49–50 | 平板热流与 Van Driest 理论比较（δ=0.70/1.50 in） | 图题 |
| Figures 51–56 | 整流/凸条热流比、油流关联、轴向分布 | 图题 |
| Figures 57–60 | 单圆柱、边界层影响、驻点理论与实验 | 图题 |
| Figures 61–65 | 串列/偏置圆柱热流比与关联 | 图题 |
| Figures 66–70 | 前掠/后掠圆柱热流比与驻点理论 | 图题 |
| Figure 71 | 特定模型平板上的等 `h/h0` 线 | 图题，pp.416–418 |

图号是原文图题，页码只在能由页脚直接核对时给出；“图题可见、页码未逐页核对”不是缺失图，而是没有把推定页码写成事实。

## 8. 条件 CSV：本次转录了什么、没转录什么

新增 `docs/research/d1372_verified_conditions.csv`，共 28 行，字段为：

```text
record_id, scope, quantity, value, unit, source_page, source_table_or_figure, derivation, note
```

它只登记**原文明确写出的条件/精度/单位**，例如：

- 总体 Mach 2.65/3.51/4.44、Re/ft 1.3×10⁶～4.7×10⁶、δ 变化 10:1；
- 平板 4 ft×5 ft、84 支热电偶、10 个静压孔；
- δ=0.7/1.5/6 in 三个位置；
- Table III/IV/V 的 δ=0.70/1.50/6.00 in；
- 热流系数精度 ±10%/±15%/±20%、温度 ±2°F；
- 油流显示 M=3.51、Re/ft=3×10⁶、δ=6 in；
- `h` 单位 Btu/ft²-sec-°R；
- `Tw` 为逐点实测、没有单一固定壁温；
- Table VI 是理论表；原文共 71 幅图。

**没有转录：**

- Table III/IV/V 的逐点 `x, y, z, θ, Te/Tt, Tw, h, h/h0` 数值；
- Table II 的逐点压力/局部 Mach 数；
- Figure 49–71 的图读热流值；
- 任何由 OCR 猜测的组头与测点配对。

原因：

1. 当前取得的是官方全文 OCR，不是逐页图像；表格列在 OCR 中被按列拆散，同一行的 `x/y/z/θ/Tw/h` 不能可靠对齐。
2. 原文表格很大（Table III 68 页、Table IV 126 页、Table V 103 页），逐点 OCR 误读会把错误数值带进项目。
3. 图 49–71 是等值线/曲线图，不是原始测量表；从图上读数不能冒充“原文表格数据”。
4. 任务要求“只有原文条件明确时才转录 CSV”；条件明确，逐点数值的对应关系不明确，所以只转录条件。

如果后续要建逐点验证集，建议的转录协议：

- 用 PDF 页图或清晰的表格扫描，不用当前 OCR 列序；
- 每个点记录表号、节号、页码、`M`、`R`、`δ`、`x/y/z/θ`、`Te/Tt`、`Tw`、`h`、`h/h0`；
- 至少两人/两次独立转录并逐字段比对；
- 先做 1 个表节的小样本，确认列序和单位，再扩展；
- 图读值单独标为“figure-digitized”，不得与表列值混用。

## 9. 对项目当前里程碑的可用性判断

**可用之处：**

- 这是一份公开、可追溯、包含 20 种突起/硬件构型的超声速湍流边界层热流与压力数据源，适合作为“突起干扰热流”方法的候选验证数据。
- 它给出了明确的 M、Re/ft、δ、模型几何、测点位置、精度和 `h/h0` 定义；比只给一条经验倍数曲线更适合做逐点比较。
- Table III/IV/V 覆盖 δ=0.70/1.50/6.00 in，可用于检查“突起高度/边界层厚度”这个控制参数；Figure 8/9/13/14 给出模型类型。

**不能直接宣称：**

- 不能宣称“已经用 D-1372 验证通过”。本次没有做任何求解器对比，也没有转录逐点测量值。
- 不能把 `h/h0` 的图读峰值当成原始测量表；也不能把不同δ、不同M/R的表节点混成同一工况。
- 不能用这份数据直接覆盖当前 20° 箭体、真实飞行、高超声速、分离/再附、真实气体等区域；原文是平板+附着突起的超声速风洞试验，M 2.65～4.44。
- 不能忽略原文的瞬态壁温法假设和 15%～20% 的绝对 `h` 相位差；做定量比较时必须匹配 `Tt` 探针位置、`Tw` 历程和同一阶段平板基线。
- 本任务没有新增物理解面积；新增正式合格面积 0 m²。

**建议的下一步（未执行）：**

1. 先取 PDF 的清晰页图，核对 Table III 的一个小表节（例如 flat plate alone 或 1-inch×2-inch stiffener）的列序和单位。
2. 选择与当前几何/工况最接近的 1～2 个构型，重建模型几何和测点坐标；只对原文实际测量的位置评价误差。
3. 将 `h/h0` 作为主比较量，绝对 `h` 作为辅助；分别报告 M、Re/ft、δ、`Tw`/`Tt` 的匹配情况。
4. 如果无法匹配 `Tw` 或 `Tt`，只能做趋势/量级比较，不能写“验证通过”。

## 10. 修改文件、命令与证据

新增：

- `docs/research/opencode_protuberance_d1372.md`：本文件。
- `docs/research/d1372_verified_conditions.csv`：28 行条件/精度/单位登记表。
- `docs/research/ntrs_19630001117_metadata.json`：NTRS 记录的精简元数据（题名、作者、报告号、日期、下载链接、访问日期）。
- `docs/research/evidence_manifest.json`：本任务文件的字节数和 SHA-256，以及全文文本的 SHA-256。

未修改：

- 任何求解器、网格器、生成器、v100/Q41/v99 数据；
- `outputs/`；
- 任务一/任务二的代码和证据。

实际命令（关键步骤）：

```powershell
# NTRS API（浏览器式 User-Agent）
Invoke-WebRequest -Uri 'https://ntrs.nasa.gov/api/citations/19630001117' `
  -Headers @{'User-Agent'='Mozilla/5.0 ...';'Accept'='application/json'} `
  -UseBasicParsing -TimeoutSec 60

# 官方全文 OCR
Invoke-WebRequest -Uri 'https://ntrs.nasa.gov/api/citations/19630001117/downloads/19630001117.txt' `
  -Headers @{'User-Agent'='Mozilla/5.0 ...'} -OutFile $env:TEMP\19630001117.txt

# PDF 只做 HEAD，确认大小
Invoke-WebRequest -Uri 'https://ntrs.nasa.gov/api/citations/19630001117/downloads/19630001117.pdf' `
  -Method Head -Headers @{'User-Agent'='Mozilla/5.0 ...'} -UseBasicParsing -TimeoutSec 60
```

本任务没有运行 MATLAB、CFD、热流或壁温求解；没有重跑 v99/Q41/v100；没有清理 `outputs`。

## 11. 可追加到 PROJECT_GOALS.md 的执行记录

```markdown
- **2026-10-30 NASA-TN-D-1372 原文核查（已测试的检索与条件登记；未做实验验证）**：取得 NTRS 正式记录与官方全文 OCR（1,601,061 字节，SHA-256 `665C95E4852B4B8845CDD6FCA95AF9804123D1E50DCD364BB14D979E09F6F218`）；PDF HEAD 确认 226,898,018 字节，全文 GET 因传输慢未在本任务窗口完成，未把局部文件当证据。核对原文为 Burbank/Newlander/Collins，NASA TN D-1372，1962-12，Langley；20 种构型，M=2.65/3.51/4.44，Re/ft≈1.3×10⁶～4.7×10⁶，δ 变化约 10:1。三张热流表按 δ=0.70/1.50/6.00 in 分组（Table III p.47 起、IV p.115 起、V p.241 起），Table VI p.344 起为无限长圆柱理论表；Figures 1–71 在 pp.348–418。新增 `docs/research/d1372_verified_conditions.csv` 28 行，只登记原文明确条件/精度/单位；逐点 `h/Tw/p` 因 OCR 列序不可靠未转录，图读值不计入原始测量表。原文热流系数精度 ±10%/±15%/±20%，温度 ±2°F；瞬态壁温法含等温蒙皮、忽略横向导热/内部传热/辐射假设，绝对 h 的探针位置相位差为 0.70-in 边界层约 15%、1.50-in 约 20%，`h/h0` 为原文强调的可比量。未运行任何 CFD/热流/壁温求解；新增正式合格面积 0 m²；不写“实验验证通过”。证据 `docs/research/`。
```

## 12. 明确未验证/未交付

- 未完成 PDF 全量下载和逐页图像核对；所有页码来自官方 OCR 的印刷页码标记，不是 PDF 物理页序号。
- 未逐页核对 Table II 的结束页、Table VI 的结束页、Figure 71 的结束页；这些范围由相邻部分推定。
- 未转录 Table III/IV/V 的逐点数值，未建立逐点验证集。
- 未把原文数据与当前求解器、20° 箭体或任何 CFD 结果比较。
- 未证明原文数据适用于当前飞行器、当前物性、真实气体或高超声速条件。
- 未新增物理解面积、未新增正式合格面积。
