# 工程算法报告交付说明

原版编制于2026年11月12日，2026年11月16日补充解释。报告正文为`docs/ENGINEERING_ALGORITHMS_REPORT_2026-11-12.md`，Word版为本目录的`RocketAeroHeat_工程算法与适用边界_20261112.docx`，保留原文件名便于继续引用。

补充版共26页，包含21个方法族、56组原生可编辑Word公式和7张表。新增第14章按A01至A21说明用途、原理、输入输出及适用边界，并解释整条计算链的分工，覆盖NASA7、外流波系、流管、边界层、换热和局部关联。

原第1至13章的科学证据仍按2026年11月12日编制基线理解：v231工程源码、编制期间新增v232记录和明确标注的独立工作区资料。补充版未重做历史来源或数值验证，没有运行新的MATLAB、CFD、引擎构建、物性验证或浏览器验收；没有改变求解器、参数、原面面积或物理资格。最新项目进展以`docs/PROJECT_GOALS.md`为准。

`source_manifest.json`原样保留，绑定原版102份本地源码/资料的长度、SHA、构建时Git身份与56组LaTeX公式；其中DOCX的SHA属于原版，不是补充版。`revision_20261116.json`记录本次正文、构建器和补充版DOCX的SHA，逐项确认56组公式未改，并列出A01至A21解释编号；不声称102份当前文件重新核验通过。`source_notes/`原样保留4个小型来源说明/条件表，没有复制候选源码、工程包或大型结果。候选自身的“通过”文字仍为该来源的报告，不自动成为生产验收。

清单SHA对应读取时本地文件原始字节，不是Git对象ID；普通源码在另一系统检出时可能因换行策略而改变文件SHA。四份冻结笔记按原始字节保留，并使用本目录.gitattributes禁止其换行规范化。

`document_qa.json`记录补充版Word的公式数、分页、逐页目视和纸面安全区检查；原版QA另存于`document_qa_20261112.json`。内部PDF和PNG只用于排版检查，补充版位于被忽略的`outputs/engineering_algorithm_report_20261116/render/`，不作为额外交付或物理证据上传。

## 重建方式

`scripts/buildEngineeringAlgorithmReport.py`从正文生成Word，需要python-docx、lxml、latex2mathml与Microsoft Office的MML2OMML.XSL。使用`--explanatory-revision`保留原版科学清单并生成补充清单，不能省略该参数而重建历史科学基线。可用`--math-package-dir`指定latex2mathml所在目录、`--math-xsl`指定转换表；本次使用Codex已配置的Python运行时和本机Office转换表。重建会更新文档容器时间与SHA，不承诺旧DOCX容器逐字节相同。

```text
python scripts/buildEngineeringAlgorithmReport.py --math-package-dir outputs/engineering_algorithm_report_20261112/runtime --explanatory-revision
```

打包渲染脚本因未安装LibreOffice而明确失败，随后使用本机Microsoft Word只读打开DOCX、后台导出PDF，再由Poppler生成全部26页PNG。逐页查看并将附录完整安排于新页；最终未观察到公式缺字、剪裁、表格重叠或越出纸面安全区。报告使用原生OMML，未把公式栅格化。原版孤立表头修正与本次渲染失败记录分别保留。

首次数学转换因转换器不支持aligned环境而失败；转换器适配成等价matrix行列后成功。这个排版修正不改变公式和数值。全部源语言LaTeX保留于正文和公式索引，便于独立核对。

科学边界：原版中旧工程66.5626%是有计算值面积；新相态鼻部0.1790893%是候选面积；均不能从解释补充获得新的物理认证。主线执行点继续由最新PROJECT_GOALS及其具名NEXT_DECISION维护。
