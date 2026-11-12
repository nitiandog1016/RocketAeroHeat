# 工程算法报告交付说明

2026年11月12日。报告正文为`docs/ENGINEERING_ALGORITHMS_REPORT_2026-11-12.md`，Word版为本目录的`RocketAeroHeat_工程算法与适用边界_20261112.docx`。

报告包含21个方法族、56组原生可编辑Word公式、6张汇总表；最终19页。按v231工程源码并纳入编制期间主线新增v232记录、明确标注的独立工作区资料整理，区分来源条件、实际代码门限、项目验证状态、条件研究和未实现候选。没有运行新的MATLAB、CFD、引擎构建、物性验证或浏览器验收；没有改变求解器、参数、原面面积或任何物理资格。

`source_manifest.json`绑定102份本地源码/资料的长度、SHA、构建时Git身份与56组LaTeX公式。`source_notes/`只冻结4个小型来源说明/条件表，便于确认未提交工作区材料的引用身份；没有复制候选源码、工程包或大型结果。候选自身的“通过”文字仍为该来源的报告，不自动成为生产验收。

清单SHA对应读取时本地文件原始字节，不是Git对象ID；普通源码在另一系统检出时可能因换行策略而改变文件SHA。四份冻结笔记按原始字节保留，并使用本目录.gitattributes禁止其换行规范化。

`document_qa.json`记录Word原生公式数、分页、逐页目视和纸面安全区检查。内部PDF和PNG只用于排版检查，位于被忽略的`outputs/engineering_algorithm_report_20261112/render/`，不作为额外交付或物理证据上传。

## 重建方式

`scripts/buildEngineeringAlgorithmReport.py`从正文生成Word及来源清单，需要python-docx、lxml、latex2mathml与Microsoft Office的MML2OMML.XSL。可用`--math-package-dir`指定latex2mathml所在目录、`--math-xsl`指定转换表；本次使用Codex已配置的Python运行时和本机Office转换表。重建会更新文档容器时间与SHA，不承诺旧DOCX容器逐字节相同。

```text
python scripts/buildEngineeringAlgorithmReport.py --math-package-dir outputs/engineering_algorithm_report_20261112/runtime
```

当前打包渲染脚本因未安装LibreOffice而明确失败，随后使用本机Microsoft Word只读打开DOCX、后台导出PDF，再由Poppler生成全部19页PNG。逐页查看后修复首版孤立表头；最终没有观察到公式缺字、剪裁、表格重叠或越出纸面安全区。报告使用原生OMML，未把公式栅格化。

首次数学转换因转换器不支持aligned环境而失败；转换器适配成等价matrix行列后成功。这个排版修正不改变公式和数值。全部源语言LaTeX保留于正文和公式索引，便于独立核对。

科学边界：旧工程66.5626%是有计算值面积；新相态鼻部0.1790893%是候选面积；均不能从本报告获得新的物理认证。主线执行点继续由PROJECT_GOALS及v232/NEXT_DECISION维护。
