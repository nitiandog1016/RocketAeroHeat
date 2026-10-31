# v104 验证与来源记录

对应 [v104 方法说明](../../KINNEY_ACREAGE_V104.md)。这是公开公式对照和单时刻数据／界面检查，不是实验精度验收；新增物理解面积 0。

## 实际运行

| MATLAB任务 | 终态 | 身份 |
| --- | --- | --- |
| engineering_baseline01 | exit1 | 测试入口的整型角度错误，保留失败 |
| engineering_baseline02 | exit0 | 原工程链最终基线 |
| kinney_tests01 | exit0 | 较早方法元数据版本 |
| kinney_actual01 | exit0 | 较早方法元数据版本 |
| kinney_workbench01 | exit0 | 较早五步骤控件检查 |
| kinney_regression01 | 4294967295 | 主动中断长历史，不算通过、不重启 |
| kinney_tests02 | exit1 | 导出结构字段不一致 |
| kinney_actual02 | exit1 | MAT已保存但导出失败，不作最终结果 |
| kinney_tests03 | exit0 | 最终公式、图、经典体与温度身份检查 |
| kinney_actual03 | exit0 | 最终物理数值，原对照标签修复前容器 |
| kinney_workbench02 | exit0 | 五步骤实际控件、显式两次求解和默认回放 |
| kinney_export04 | exit0 | 温度合同、统一色标，原标签修复前导出 |
| kinney_provenance01 | exit0 | 逐面原来源修复、物理逐位不变、最终导出 |

`processes/` 和 `logs/` 保留当时原字节，包括主动中断说明。`provenance.json` 绑定 13 项终态、139 个确切源码版本及输入 SHA；`sourceStableAfter` 是各文件布尔检查结果，不能当作散列值。

源码闭包是工程函数名的保守词元闭包，包含注释、字符串和显式 `evalc` 函数名；它不是运行时调用轨迹，也不声称所有被快照的无关文件均归档。`source_versions/*.raw` 用完整 SHA 核验；初次失败入口从修复前文本重建后，与当时记录的完整 SHA 相同才接受。Git 属性禁止这些证据自动转换换行。

## 数值与展示证据

- `independent_final_physics_audit.json`：保存最终数值的独立公式／图约束检查，242,140 项；早期 `independent_audit.json` 保留其原身份。
- `summaries/`：基线、最终定向算例、真实对照、五步骤和标签修复的小型原摘要。
- `baseline_provenance.mat`：标签修复保存的原方法来源，不能代替完整物理 MAT。
- `viewer_provenance_contract.json`：最终资源数据探针／截面检查；较早 `viewer_contract.json` 不覆盖修复后的版本。
- `published_viewer_provenance_assets.json`：最终 payload SHA 与五个静态资源的 Git 提交／SHA。`published_viewer_assets.json`、`published_viewer_final_assets.json` 保留早期资源身份。
- `browser/final_workbench.png`、`nose_original_source.png`：真实浏览器截图；AX 记录显示来源、滑块和未知壁温。`console.json` 为本次捕获的 warn/error 日志，结果为空。
- `browser/section_minus1_dom_export.csv`：实际点击后，从可见下载链接捕获的生成内容。`browser/csv_audit.json` 将其与保存 MAT 和原三角形交点独立比较。**不是已证实的原生下载文件**；三次下载等待超时如实记录。

最终浏览器资源在本地 `outputs/v104_validation/published_workbench05/`。`case-data.js` 与最终导出硬链接共用，禁止原位改写；应向新目录重新导出。其 SHA 为 `81bc1e610b7875a7d8ba75b8ebe86fecf995d16ab41f251e81c25a94dfc1e1d8`。

`browser/session_summary.json` 汇总实际动作及证据散列。截图证明渲染状态，CSV/Node证明数据，MATLAB控件流程证明计算操作；三者不能互相替代。

## 存储与下一起点

大型 `baseline02`、`actual03`、经典体/工作台 MAT 和研究 PDF 仍只存本地，GitHub同步不等于已备份。没有新建工程包或重复隔离工程；标签修复只保存新增小来源。旧失败与唯一结果未删除。

后续从[现行目标](../../PROJECT_GOALS.md)第4节B开始，复用现有缺口/父链证据，不重跑这些已完成任务、旧球鼻、Q41或长历史。不以方法分歧或标签修复宣称新增面积。
