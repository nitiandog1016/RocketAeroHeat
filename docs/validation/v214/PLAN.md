# v214事前计划：实际NK矩阵向量积与同阶独立控制

上一目标轮为progress：v213原生边界与真实梯度/结果等价已实现并推送，原中心差分失败及正式入口缺项明确。总目标active。继续原全面积/科研UI路线，不改成仅接口目标。

一次≤60秒定向复用准确newtonkrylovsolver、已编译几何覆盖的encode/decode、原nk-config及v213逐量文件，核对实际Jv调用链、设置恢复、时间层和变量尺度。原生入口是evalJacobianVectorProduct→evalComplexMatVecProd→setJacobianEvalSettings→U[1]守恒量扰动/decode→evalResidual(1)→zed，再恢复二阶残差设置。当前Jacobian=1/Residual=2、两冻结限幅false、无blending/smoothing。不能把一阶近似Jv当二阶精确导数，也不能从v213非正式入口失败断定生产Jv失败。

## 唯一控制与预算

固定v213 control02的网格、初场和配置，无新网格/prep/物性helper；复用真实13列/416棱柱/32层/2微米、原气体和边界。新增隔离nk_probe，不改成熟内核/已有probe及失败包。

沿用预定三种速度。方向由原p/T/U/SA方向经原生encode_conserved转换为每单元守恒量切向，另用解析乘积/能量式独立审核。调用实际evalJacobianVectorProduct(sigma=1e-30,phase=0,steps=0,false)，不自己替代该函数；只分配此单算子所需zed向量，未调用完整Krylov/ILU或伪时间对角，不能授线性求解/区域峰值资格。

独立实数基准及±1e-4使用相同setJacobianEvalSettings与U[1]/decode/evalResidual(1)，受体变、供体冻结。驻留与快照两路径均采基准、正式Jv与±实扰动。保留每原生六通量/六残差/有效梯度和母单元ID，核对接口缩放差≤1e-11、正式zed与dUdt切向≤1e-11、同阶中心差分≤1e-5。原v213二阶失败不追改；本轮不扫扰动/通量/CFL/初猜或放宽门槛。未过则依据完整分量复审具体组合，不无界追加控制。

源/输入核对与准备≤60秒/1GiB/32MiB；一个完整helper构建≤150秒/父180秒/6GiB/512MiB/2CPU，禁止第三个完整包或自动重试；若纯编译接口错误，保留并在同150秒总预算内最多一次确定性修正。一次查询≤60秒/父90秒/2GiB/32MiB/1CPU；无时间/牛顿发展准入。所有进程输入SHA、组清空、采样RSS明确，旧包不覆盖。

接口成立后尽快生产文件桥接/母ID/未用槽守卫、物理终端/反馈与完整区域资源；工程与UI独立最高，并发Web保持。原面/面积/责任/32层/历史19层尾缘/全部门槛与29保护保持。不新代理/computer-use，无包/独有删除，小源码证据公开、大资产E盘。
