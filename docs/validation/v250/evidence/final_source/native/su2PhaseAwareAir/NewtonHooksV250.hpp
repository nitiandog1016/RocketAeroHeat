#pragma once
// 模块职责：只读暴露原生NK的阶段和方向；在已声明CNewtonIntegration后包含，不改变数值对象布局。
namespace rah_newton_observer {
using Direction = CSysVector<CNewtonIntegration::Scalar>;
using Callback = void (*)(const char*, CGeometry*, CSolver**, CConfig*, const Direction*, double, double);
inline Callback event = nullptr;
}
