#pragma once
// 模块用途：只给本理想空气模型的物性域拒绝附上查询身份，不改变合法状态或异常恢复规则。
#include "IdealAir.hpp"
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>

namespace rah_query {
// 模块一：线程局部作用域。局部RCM编号与全局编号分开，未知编号不得冒充0号节点。
struct Context {
  const char* phase="unscoped";
  bool localKnown=false, globalKnown=false, markerKnown=false, innerKnown=false;
  unsigned long local=0, global=0, marker=0, inner=0;
};
inline Context& context() {
  static thread_local Context current;
  return current;
}
struct Guard {
  Context previous;
  explicit Guard(const Context& value):previous(context()) { context()=value; }
  Guard(const Guard&)=delete;
  Guard& operator=(const Guard&)=delete;
  ~Guard() { context()=previous; }
  void phase(const char* value) { context().phase=value; }
};

// 模块二：成功路径直接返回同一内核结果；仅domain_error加说明，其余异常原样传播。
// 不打印每个可恢复试探，不产生无限日志；调用者原有的回退/终止继续决定后续行为。
inline rah_air::State evaluate(const std::string& mode,double first,double second) {
  try {
    return rah_air::state(mode,first,second);
  } catch(const std::domain_error& error) {
    const auto& c=context();
    std::ostringstream message;
    message.imbue(std::locale::classic());
    message<<std::setprecision(17)<<"RAH_PROPERTY_QUERY mode="<<mode<<" phase="<<c.phase
      <<" localKnown="<<c.localKnown<<" local="<<c.local
      <<" globalKnown="<<c.globalKnown<<" global="<<c.global
      <<" markerKnown="<<c.markerKnown<<" marker="<<c.marker<<" innerKnown="<<c.innerKnown<<" inner="<<c.inner
      <<" input1="<<first<<" input2="<<second<<" reason="<<error.what();
    throw std::domain_error(message.str());
  }
}
} // namespace rah_query
