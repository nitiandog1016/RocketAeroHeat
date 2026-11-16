// 模块职责：只读层间观察接口，不改变原生类布局和数值公式。
#pragma once
class CGeometry; class CSolver; class CConfig;
namespace rah_fas_observer {
using Callback=void(*)(const char*,CGeometry*,CSolver**,CGeometry*,CSolver**,CConfig*);
inline Callback transfer=nullptr;
}
