# v92验证与当前状态

索引为`progress_repreview.json`，方法和限制见`../../TEMPERATURE_REPREVIEW_V92.md`。

当前真实检查点Q40，累计18次接受/78次主NS；另计本轮1次准备NS和2次独立探针NS。1/8步实际接受，RMS下降12.4336%，仍未稳态，新增正式合格热流面积为0。整箭本次仅需一次受限重算；多轮新条件补全由小域制造案例验证。

验证分开记录：小域稠密方向/集合/公开温度导数；真实C/W与投影支撑；逐轮全部CV/界面输入分区及实际公开抽样；全部19752开放储库面、墙/Fourier/五量和历史；原MAT到面元云图、曲线/CSV/方法/给定与未知身份；浏览器实际操作。内部Riemann公开参考只抽查266/360面，不能写成全部内部解独立重算。浏览器初始化错误照实记录，不宣称零错误。

`raw/`保留失败日志、旧测试/合并器源字节及旧整数scale小输入，均非新增物理试算；最终空集合身份检查单独记录。图像包含布局修复前和缓存未更新时的证据，最终界面图为`q40_workbench.jpg`、`q40_heat.jpg`和`q40_diagnostics.jpg`。

六项Q39派生矩阵在全部实际审计结束后回收，实收约2.10GB。删除前后SHA、清单与独立核对见`obsolete_cache_*.json`。保留Q39/Q40、几何、方向、检查点及结果；旧完整方向审计不能直接加载已清理的缓存，历史重放须按对应源码和保留输入在新目录重建研究及新来源绑定，不保证历史MAT文件SHA相同。

小域复现不需要这些大矩阵：在MATLAB运行`addpath('tests'); testRocketTemperatureRepreview('outputs/repreview_new'); testRocketTemperatureLookaheadIntegration('outputs/lookahead_legacy_new');`，随后运行`python tests/auditRocketTemperatureRepreview.py outputs/repreview_new`。输出目录须为新目录，输入随仓库保留。GitHub交付完整源码及必要小资料，不默认生成本地工程ZIP；总目标仍在进行中。
