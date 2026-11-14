// 模块职责：只读观察回调的窄接口；不增加原生类成员，不参与数值决策。
#pragma once
class CGeometry; class CSolver; class CConfig; class CNSVariable;
namespace rah_single_observer {
using StageCallback=void(*)(const char*,CGeometry*,CSolver**,CConfig*,unsigned short);
using RejectCallback=void(*)(CNSVariable*,unsigned long,const double*,const double*,bool,const char*);
inline StageCallback stage=nullptr;
inline RejectCallback rejection=nullptr;
}
