// 模块职责：原生配置解析器识别稳定反算选择，值由 COptionBase 保存并按 SU2 规则继承。
#pragma once

class CStableCaloricOption final : public COptionBase {
public:
  std::string SetValue(const std::vector<std::string>& tokens) override {
    if (tokens.size() != 1 || tokens.front() != "STABLE_PIECEWISE") {
      return "RAH_CALORIC_INVERSE expects STABLE_PIECEWISE. Omit this option for TOMS748.";
    }
    return COptionBase::SetValue(tokens);
  }
  void SetDefault() override {}
};
