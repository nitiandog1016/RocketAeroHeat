// 模块职责：在实际一般EOS边界提交通量前调用只读出口审核，不修改原生数值公式。
#pragma once
namespace rah_outlet_v240 {
using Hook=void(*)(unsigned long,unsigned long,double,double,const double*,const double*,const double*);
inline Hook hook=nullptr;
}
