# 固定组分理想空气模块

显式模型为 `Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1`。`IdealAir.hpp`提供无外部库依赖的C++17状态核，MATLAB入口位于项目根目录的`idealAirThermodynamicState.m`。原生验证查询器另链接固定CoolProp库。

## 使用

```matlab
state = idealAirThermodynamicState('PT',101325,300);
restored = idealAirThermodynamicState('rhoe',state.rho,state.e);
```

```cpp
#include "IdealAir.hpp"
const auto state = rah_air::state("PT", 101325.0, 300.0);
const auto transport = rah_air::transport(state, 0.0);
```

六种状态对为`rhoT`、`PT`、`rhoe`、`Prho`、`hs`、`Ps`，参数顺序与名称相同。全部使用SI，能量和焓单位J/kg，熵和比热J/(kg K)，密度kg/m³，压力Pa，温度K，黏度Pa s，导热系数W/(m K)。MATLAB支持同尺寸数组与标量扩展。

`e`为热内能。流动调用者按`e = rhoE/rho - |u|²/2 - k`恢复SST热状态。能量与熵采用固定CoolProp Air的OTH参考，重启动转换须先明确原模型参考。反解温区60—2000K，越界报错。仅`Prho`的派生温度允许8eps×2000K以内的端点舍入恢复。

热力学输出包含cp、cv、e、h、s、声速以及固定其他自变量的压力、温度、焓、熵偏导。工程输运采用Sutherland、Pr=0.72和Prt=0.9，导热系数导数包含dcp/dT，并冻结给定湍黏度。输运的实验适用范围须随工况另外核对。

## 数据与验证

`air_ideal.json`保留固定来源理想项和OTH偏移，MATLAB直接读取。C++常数与该数据通过原生库及跨语言对照验证。气体常数为8.31451/0.02896546 J/(kg K)，缩减温度132.6312K，缩减密度302.622436442kg/m³。

来源为CoolProp7.2.0提交`98b3523d5daa98454618d381d2ae53f7471d216b`的Air.json，原文件SHA256为`c5e0d064d83fc91ed35d7354cea1fa38fea832a09d2e85faf4336e0017c33616`。附原CoolProp许可。理想项源于Lemmon等2000，式18、式24和Table12。物性范围内仍需核对固定组分、气相与理想气体假设，实际压力效应及化学反应另行建模。[原论文](https://doi.org/10.1063/1.1285884)

本模块已通过v153代数、原生库对照和MATLAB入口检查，完整引擎接入待完成。测试记录和性能口径见项目文档`docs/FAST_IDEAL_AIR_V153.md`。
