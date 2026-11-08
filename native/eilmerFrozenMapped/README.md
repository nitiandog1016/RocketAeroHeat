# 冻结mapped边界的隔离原生控制

固定上游：[GDTk/Eilmer](https://github.com/gdtk-uq/gdtk/tree/f53f4609a0331d48efee69a4e4f3c3598378cc03)，GPLv3或以后；依赖组件各自许可。本目录适配沿用该授权，许可证原文见项目已有v185来源证据。保留原生流动算法，不分发第三方运行库。

`lmr/bc/ghost_cell_effect/rah_frozen_mapped.d`负责独立捕获供体数据、受体ghost对应和原生六组交换。`probe.d`负责真实近壁小片的梯度、原生通量/残差及原始量方向导数控制；不包含时间推进。使用v210相态空气factory覆盖和v196已证几何修正，旧包/共享源不改。

v213已有基准与复数响应等价，但中心差分失败；正式NK守恒变量Jv未测。本接口尚未接生产NPZ/factory，母ID/轮次/完整实复数及未用槽守卫未齐，禁止直接用于区域生产求解。历史运行入口已关闭；新任务必须按目标文档事前登记，不能直接重跑历史构建。

完整范围与失败见[原生对照](../../docs/HYBRID_FROZEN_NATIVE_OPERATOR_V213.md)，下一项见[执行点](../../docs/validation/v213/NEXT_DECISION.md)。
