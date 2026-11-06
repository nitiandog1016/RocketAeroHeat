# v171 原输运表有界补证

2026-11-06。本次只补资料证据，原始开始时间09:23:26 UTC、结束时间09:29:00 UTC，实际耗时334秒，超出300秒预算34秒，时间预算未遵守。最多3组定向检索的数量上限已遵守，三组已完成并停止。此前称“5分钟内”不准确，本处及回执明确纠正，原时间和记录保留；本次纠正没有再检索或提取。未触碰正在运行的候选、旧data、配置、门槛或其他文件，无CFD/构建/Git操作。

## 结论

**未取得1950年原表2.39/2.42的可核对数值页，150–350 K原表数值新增0；原输运差异仍未闭合。** 未用现代数据库、现代Sutherland参数或二手教科书数据冒充原表。没有生成物性CSV，也没有据此修改本次候选。

已定位[NIST官方NBS Circular 564书目](https://www.nist.gov/publications/circular-bureau-standards-no-564tables-thermal-properties-gases-comprising-tables)，1955年出版，作者包括Ralph L. Nuttall，内容范围包括空气的热力学与输运性质。其[DOI](https://doi.org/10.6028/NBS.CIRC.564)重定向至[官方PDF](https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular564.pdf)。这是未来可追溯的**官方后续汇编线索**，本轮未核验表页和编修说明，不能称其逐字重印或数值等同TN3344的1950参考26。

## 三组检索与证据界限

1. 原作者/表号定向检索：Morey、2.39、NBS-NACA及Nuttall、2.42、1950。找到NASA原始报告对表号的引用，但没有原散页表。
2. 官方域及准确题名检索：NIST官方书目找到Circular 564；NASA报告再次确认两表的题名和1950年身份。
3. 官方PDF路径及564与原表号关联检索：未得到能证明原表数值保留的表页或编修说明。没有继续扩展第四组检索。

TN3344本地文本参考26仍为Morey的Dry Air Coefficients of Viscosity、Nuttall的Dry Air Thermal Conductivity，NBS1950。旁证[NACA TN2843](https://ntrs.nasa.gov/api/citations/19930083593/downloads/19930083593.pdf)参考18/19列出2.39/2.42和December1950；[NACA TN3163](https://ntrs.nasa.gov/api/citations/19930083885/downloads/19930083885.pdf)参考15/19写July1950。这些仅是检索返回的书目信息，未当作物性数值证据；月份差异也不能在未见原散页时自行消解。

## 获取限制和停止条件

官方PDF的HTTP响应为200，Content-Length为53,394,377字节。浏览工具因文档大小拒绝读取。限时本地恢复依次遇到fitz/pypdf模块缺失及指定位置的pdftotext.exe不可用；失败记录保留在`tests/data/naca_tn3344/transport_v171/search_receipt.json`。没有为此安装依赖、重复检索或扩大文献支线，未保存完整大PDF，也未宣称已目视检查官方表页。

未来若另行进行取表工作，应先核验1955汇编关于1950散页的编修说明，再逐页确认空气黏度、导热系数的温标、归一化基准、单位和150–350 K数据；仅可把未证明与1950相同的数值标成1955独立对照。当前`originalTransportTableEquivalent=false`保持；v170物性归约比算子只证明给定两模型时的代数换算，不证明原表已取得。
