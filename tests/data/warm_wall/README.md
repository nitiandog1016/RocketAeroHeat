# 暖壁附着支的独立数值参考

`independent_inputs.json` 保存 M4、20 km、零攻角球鼻外缘中的四个局部站（40°、45°、50°、54.797931°），与原局部相似内核输入相同。气体常数 R=287.05287 J/(kg K)、gamma=1.4；NASA7干空气物性、Sutherland黏度、Pr=0.72。输入并不是实验测量，外缘压力使用球面修正牛顿工程近似。

`independent_adiabatic_summary.json` 是独立 SciPy `solve_bvp` 的直接绝热边界解：壁面直接施加 `(l/Pr)*g_eta+omega=0`（壁面F=0故n=0），求出壁面总焓后反解Taw。没有调用MATLAB生产内核、MATLAB的Tw热流根搜索或诊断副本。相同输入/物性与NASA TP-3271局部相似方程是有意共用的物理模型，不是独立实验或独立外流证明。

产生脚本和原页证据在阶段交付包的 `validation/nose_branch/independent_adiabatic_bvp.py`、`NASA_TP3271.pdf`、`NACA_TN3325.pdf` 中；开发树位置为 `outputs/v54_validation/nose_branch`。脚本使用 SciPy1.7.3/NumPy1.21.6，容差1e-9，etaMax=12/20，记录原始过冲、导数、残差、温度、剪应力及Pr1控制。

经典 NACA TN-3325 印页17–18、22解释暖壁有利压梯下的粘性附着支可以速度过冲；其 Pr1/线性黏性物性条件不能作为这里的Pr0.72/NASA7数值热流基准。这个来源仅支持“速度全域单调不是通用附着条件”。
