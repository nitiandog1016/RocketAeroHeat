# 暖壁附着支速度过冲修正

v53球鼻局部相似模型在45°之后的绝热边界未解，根源是剖面内核把“速度全域单调且不超过外缘”当作附着支条件。300 K冷壁热流剖面仍可解，较热的绝热试探剖面则被该条件拒绝。v54修正保留方程、外缘、质量积分和各物理筛选，只更正附着支的错误定义；尚未完成球鼻/后体接口、斜流、全区域求解或实验验证。

## 原方程及来源核对

[NASA TP-3271 (1992)](https://ntrs.nasa.gov/api/citations/19930006245/downloads/19930006245.pdf) 印页4（PDF页6）的式19、21–25、27及边界条件给出下列层流局部相似极限。沿xi非相似项在本模块省略，f_eta=F，l=rho*mu/(rhoe*mue)，D=rhoe/rho=T/Te，m=l*F_eta，n=(l/Pr)*g_eta+omega：

```
m_eta = -f*F_eta + beta*(F^2-D)
omega = (1-1/Pr)*(Ue^2/H0)*F*m
n_eta = -f*g_eta
f(0)=F(0)=0; g(0)=gw; F(infinity)=g(infinity)=1
```

原页不要求F<=1或F_eta>=0处处成立。速度出现内部局部最大值时m=0，式19允许D>F²，从而F_etaeta<0；温暖低密度流层可在有利压梯作用下比外缘加速更多，随后返回外缘。

[Cohen与Reshotko，NACA TN-3325 (1955)](https://ntrs.nasa.gov/citations/19930084102) 印页17–18明确说明：暖壁有利压梯的速度剖面可以从大于1的值返回外缘；这一粘性支须与逆压梯下非唯一的无粘支区分。印页22结论2、3重复确认。该报告假定Pr=1及线性黏性温度关系，仅用于核对支类型，不作为本NASA7/Pr0.72数值热流的定量基准。

[NASA SP-8062 §2.3.2](https://ntrs.nasa.gov/api/citations/19710021703/downloads/19710021703.pdf) 说明局部相似在中等压力梯度有用，在强有利压梯/小曲率区或趋近分离时准确性下降。因此本修正仍属有范围说明的工程近似，并未由“BVP有解”证明实际球鼻精度。

## 保留的判据及诊断

`solveAxisymmetricLaminarMassProfile.m` 保留原物性80–2500 K、0<=beta<=1、输入能量/动量兼容、ODE与边界残差、正壁面剪应力及无反向流判据。外缘仍同时核对速度、总焓、静温、密度、剪切/热梯度及数值边界余量；原有限剖面积分及厚度检查也保留。球鼻模块的几何、零攻角、Kn、Re及delta/R门限未变。

删除错误的`max(F)<=1+tol`及全域`min(F_eta)>=-tol`附着门限。过冲和导数不被裁剪，新增`velocityAudit`记录完整数值域内的最小/最大速度比、过冲量、最小导数和壁面支判据。旧`massProfile.monotonicityTolerance`配置名兼容保留，现仅用于无逆流的最小速度数值容差。

## 数值证据

四站为原M4/20 km/零攻角半球头柱体数据，NASA7气体、Pr0.72。独立SciPy解直接施加绝热n(0)=0，不经过MATLAB内核或MATLAB壁温热流求根。

| 球面角度 | 独立Taw (K) | 速度过冲 (ppm) | 壁面剪应力 (Pa) |
|---:|---:|---:|---:|
| 40° | 862.789028159 | 0.297522 | 181.626271689 |
| 45° | 857.253565311 | 1.275719 | 192.495699857 |
| 50° | 850.827746105 | 5.062306 | 199.526676876 |
| 54.797931° | 843.719649046 | 17.627896 | 202.400759505 |

独立直接BVP与MATLAB诊断副本的Taw差<4e-9 K、剪应力相对差<6e-11；etaMax12到20的Taw差<3e-10 K。Pr1、Tw=T0的同站绝热控制得到g=1、q=0、最大F=1。过冲稳定发生在eta约4–5处，远离数值外界，原物性及无逆流条件均满足；属于同方程独立数值校核，不是实验验证。

2026-10-14生产修正回归已通过：`tests/testWarmWallAttachedMassProfile.m`核对四站的直接基准、实际热流根、冷壁加热/暖壁冷却、物理壁坐标焓导数、eta及BVP容差收敛、Pr1控制及真实暖壁Pr1过冲；旧`testAxisymmetricLaminarMassProfile`的Blasius/Falkner-Skan/Crocco极限和原拒绝校核也通过。开发日志 `outputs/v54_validation/nose_branch/warm_kernel_regression.log`包含较早试跑断言及最终通过，不能把整份当纯通过日志；发布复测记录为 `outputs/v54_validation/warm_kernel_release_final.log`，其终态由项目执行记录及归档清单保存。逐项量化结果在`warm_kernel_regression_summary.json`。主球鼻整链路由单独测试证明，不能从本局部站回归推断。
