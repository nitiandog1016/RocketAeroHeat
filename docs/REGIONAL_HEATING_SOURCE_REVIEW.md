# 区域气动热来源审查与待实现工程方法（v54，2026-10-14）

本文件整理已核验的公开公式原页、球鼻无量纲散点草稿及底热下一份一手论文线索。本项研究没有实现新的生产热流模型，没有把代表底热映射到底面，也没有完成实验比较。下面“可实施”表示公式及输入已识别；不表示适用性或精度已通过项目测试。其他 v54 计算修改的状态由 `docs/PROJECT_GOALS.md` 单独记录。

固定数据和标定脚本位于 `tests/data/experimental_hemisphere/`；原 PDF、关键页、点击标注图及研究日志保留在本地验证证据目录 `outputs/v54_validation/heat_source/`。本文件中的验证证据相对路径以项目根目录为起点，完整工程包的相应源页是研究资料而不是生产外流输入。

## 1. 球鼻：可作为独立工程对照的 Lees/Kemp 分布

官方来源为 Michael E. Tauber, *A Review of High-Speed, Convective, Heat-Transfer Computation Methods*, NASA TM-101077, 1989，NTRS 19890017745。

- 官方全文：https://ntrs.nasa.gov/api/citations/19890017745/downloads/19890017745.pdf
- 打印7页/PDF11：Eq42、43，局部热流与驻点热流比及压力梯度参数。
- 打印8页/PDF12：Γ 的定义，明确 Γ=1 同时用于 Eq42/43 就得到 Lees 对照；Fig5 是半球柱热流分布；Eq44 为冷壁近似相关式。
- 打印9页/PDF13：平头转角局部相似与测量比较，非相似/转角问题不能以平滑补值解决。

原页 PNG 和相同页的文本已保存至 `outputs/v54_validation/heat_source/19890017745_pdfpage_11.png`、`...12.png`、`...13.png`。以下是保留所有热力变换参数的公式转写草稿，实施前还需核对各符号的原始定义：

```text
q_w/q_ws = sqrt[(rho_w mu_w (dU_e/ds)/(beta t_e)) /
                    (rho_ws mu_ws a/(beta_s t_s))] Gamma

beta = [2 t_s (dU_e/ds)/(rho_w mu_w U_e^2 (r/L)^(2j) Gamma^2)]
       * integral_0^s[rho_w mu_w U_e (r/L)^(2j) Gamma^2 ds]

j=1 为轴对称；a=(dU_e/ds)_stagnation。
```

这里 `rho_w mu_w` 是**当地壁面状态**的密度黏度乘积，不是项目 v53 的 `rho_e mu_e`。`t_s、t_e` 也是该报告定义的热力变换参数，不可从 OCR 下标猜测后直接省略。原页的 Eq42/43 可以实施独立 Lees 对照，但必须先将 Tauber 的全部变换定义与现有 Levy–Lees 内核逐项核对；不能将省略这些参数的简式直接当成一般高焓球鼻公式接入。

冷壁工程法可用当前内部 `p_e(s)、U_e(s)`、核验的球面弧长/半径、给定均匀壁温和独立驻点热流作输入，输出独立方法的 `q(s)`。它不包含当前 BVP 的热壁绝热根，不能借此声称 BVP 已解；也不能从一条冷壁 `q/q_s` 关系伪造恢复温度。原源适用的是薄、附着、轴对称零攻角冷壁层流；热壁、分离、跨鼻部/后体接口和沿程熵层吞入均须另外证明。

### 图5来源链与草稿数据

Tauber Fig5 重印 Kemp, Rose, Detra, *Laminar Heat Transfer Around Blunt Bodies in Dissociated Air*, Journal of the Aerospace Sciences 26(7), 421–430 (1959), DOI https://doi.org/10.2514/8.8128，原 Fig3。

24个清晰可辨的圆点已手工选取，保存为 `tests/data/experimental_hemisphere/tauber_fig5_circle_points_draft.csv`；同目录 `tauber_fig5_points_metadata.json` 记录来源链、待补输入和数据身份。`tests/data/experimental_hemisphere/digitize_tauber_fig5.py` 保留点击坐标、原页坐标、线性角度轴与对数热流轴标定；`outputs/v54_validation/heat_source/tauber_fig5_selected_points.png` 可逐点复核。只选择圆点，没有把理论曲线采样当实验。程序实际运行通过结构/域检查，**尚未完成实验热流精度比较**。

标定脚本默认校验 CSV/元数据与固定原始像素选点一致，不修改数据：

```powershell
python tests/data/experimental_hemisphere/digitize_tauber_fig5.py
```

需要重新生成固定数据时使用 `--write`。源图存在时，可使用 `--render-overlay` 输出点击复核图；可用 `--source-dir` 指定本地源页目录和 `--overlay-output` 指定输出位置。源图仅用于研究复核，不是求解所需的外部 CFD 数据。

名义点击半宽是原渲染页2像素，约0.442°及1.04%热流比；它仅描述这次选点。仪表覆盖约5–10°的弧段平均、模型连接误差及原实验测量不确定性均大于此选点量，未重建。

原论文官方全文返回403；可读的镜像 OCR 提供了来源链的工况线索，但未把它视为经过原页逐项核验的实验输入：

- 镜像：https://www.scribd.com/document/418715814/Laminar-Heat-Transfer-Around-Blunt-Bodies-in-Dissociated-Air
- 原实验有两组传播激波 Mach 范围：约10.9–13.5及7.65–8.5，初压分别约1 cmHg及0.5 cmHg；这些是 `Ms`，**不是模型自由流 Mach**。
- 名义计算例 `Ms=12`，驻点压力约2.242 MPa、温度约6945 K。热化学状态涉及高焓解离。
- 圆点系列身份、每点具体运行、壁温、入口气温、表面催化性等尚未通过原图及实验报告核验。CSV 的系列身份故意保留未知。

这些散点目前只能做图形/来源准备；不能作为 M4、20 km、Tw=300 K 冻结 NASA7 算例的已完成独立实验验证。下一步有意义的比较是先复原原 shock-tube 工况/壁温/测点平均，然后运行可支持相同热化学假设的模型；不将实验角度插到当前工况里计算误差。

## 2. 底面：找到径向梯度原论文，不能继续只依赖代表 q

已取得 Crossref 官方书目信息：Bruce M. Bulmer, *Radial Base Heat-Transfer Gradients in Turbulent Flow*, Journal of Spacecraft and Rockets 11(11), 798–800 (1974), DOI **https://doi.org/10.2514/3.62179**。`outputs/v54_validation/heat_source/radial_base_crossref.json` 保存原返回记录。

它直接讨论底面径向梯度，应优先于再搜一般底热综述。AIAA官方本轮403，未读正文，尚不清楚它采用何种 r/r_b 函数、传感器位置、前体类别、壁温、是否可配合 SAND93 代表热量归一化。不能凭题名发明多项式、把 SAND 的代表 q 当中心/平均值，或把另一几何的径向曲线缩放到当前箭体。

下一次应取得这三页原文，至少核对：

1. `q(r)/q_ref` 中 q_ref 是中心值、指定测点值还是面积平均；
2. 径向函数是否给出全圆底、底唇和 r=0 极限，或者仅报告几处梯度；
3. 前体/底半径、肩部边界层类型与厚度、Mach/Re、攻角、喷流/支撑、壁温；
4. 原数据能否与 SAND93 的热流定义对应，归一化是否涉及总焓差；
5. 传感器足迹及不确定性。只有这些明确后才可建立面积平均/径向积分一致性测试和逐面元 q 来源。

## 3. 已核实的分离区公式：能给均值或峰值，不能代替完整热流分布

### NACA TN-3792，Chapman1956

官方：https://ntrs.nasa.gov/api/citations/19930084514/downloads/19930084514.pdf

打印23页/PDF25明确其能量积分只给总体平均热率，不能计算局部壁面热流；局部分布依赖回流运动和壁形。附录C打印34页/PDF36给轴对称分离混合层的总体热率及表面平均。不能把这里的轴对称平均误读成平圆底径向公式。`outputs/v54_validation/heat_source/NACA_TN3792_pdfpage_25.png`、`...36.png` 已核对。

可推进的独立内核是：给定上游层流边界层/自由剪切层状态和分离区域几何，解附录给出的混合层能量积分得到总体 Q/平均q及恢复因子，并为后续局部回流求解提供整体守恒约束。它没有足够信息直接涂满当前圆底，也不适用于把湍流底区当层流。

### Tauber1989，打印16–18页/PDF20–22

有来源的工程关系包括：

```text
Eq60: C_h,peak/C_h,upstream = (p_peak/p_upstream)^0.85
       C_h = q/(h_aw-h_w)；不同实验约0.73--0.85指数，受壁温/流态影响。
Eq63: St_max = 0.664 (Re'_x)^(-0.5) (Pr')^(-0.667), laminar/laminar shock interaction
Eq64: St_max = 0.0576 (Re'_x)^(-0.2) (Pr')^(-0.667), laminar/turbulent or turbulent/turbulent
Eq65: St_ave = 0.48 (Re_L)^(-0.4) (L/H)^0.2, planar turbulent cavity average
Eq66: St = 2.07/Re_delta, laminar separated/cavity leeward heating
Eq67: (Nu/Pr)_b/(Nu/Pr)_c =35.5 (Re_b/Re_c)^2.2, cone turbulent base representative
```

Eq66 的原来源是 Lamb1980, *Convective Heat Transfer Correlations for Planar, Supersonic, Separated Flows*, DOI https://doi.org/10.1115/1.3244287。它是平面分离/腔体关联，不是轴对称圆底径向求解。Eq63/64 是**峰值**，并不定义再附点位置或上游分离长度，也没有全区域空间形状。

所需内部输入是上游已求边界层状态/厚度、流态，压力跳跃、已确定再附位置、当地参考焓/温度/黏度及壁温。当前项目缺少闭合的分离长度/回流/再附位置；不应重新播种流线越过未解带，然后套峰值当全区 q。可先做受支持压缩坡/入射激波的独立峰值模块与原图比较，再实现实际位置/区间的流场闭合。

## 4. 真正含空间底热及温度的下一内部求解器依据

J.C.S. Meng, *A New Technique for Calculating Reentry Base Heating*, **NASA TN-D-7508**, 1973，官方全文：https://ntrs.nasa.gov/api/citations/19740005592/downloads/19740005592.pdf。

报告给近尾区完整 Navier–Stokes、其余区边界层方程，Telenin条带离散及多重打靶；需要由内部前体给出入口热边界层，以及下游近尾迹剖面。它输出局部 q(y)、底区温度/压力；Fig12（PDF38）是真实空间热系数曲线，可以成为未来内部尾迹求解器的数值校核。它是二维特定几何数值结果，**不能作为轴对称箭体的一条通用工程径向函数**。

复算基准原文：零攻角，M_inf=5.31，T_inf=217.43 K，Tw=255 K，Re约0.8e5--1.2e5，Pr=1，线性黏度/导热系数。官方元数据及全文保存在 `outputs/v54_validation/heat_source/19740005592_metadata.json` 和 `19740005592.pdf`。这套基准是数值方法复算，不是独立实验验证。

若径向梯度原论文不足以形成闭合，物理上直接推进完整底热的路径是内部轴对称回流/自由剪切层/壁面能量求解，使用当前内部肩部状态、热边界层剖面、给定Tw作为入口/壁面边界，输出实际 p、T、q并做网格/下游边界收敛与守恒。二维 TN-D-7508 可先用于独立方法校核，轴对称扩展仍需方程/基准；不能把其2D精度继承到真实3D。

## 本轮状态与明确下一步

已完成官方全文/页图保存、24个球鼻实验圆点草稿及可复现结构检查、底面径向原论文 DOI 定位。没有运行生产模型热流比较，没有导入外部CFD作为程序求解前提，没有新增全面积覆盖。

建议下一批实现依次是：核对球鼻BVP热壁附着支并做现有工况测试；建立严格可选的冷壁 Lees/Kemp 对照；取得 Bulmer1974 三页径向原文后决定底热关联能否实装；若资料不支持全底面空间分布，就实现内部轴对称回流能量闭合。峰值/平均分离关联可提供验证约束，不用其代替实际全区热场。
