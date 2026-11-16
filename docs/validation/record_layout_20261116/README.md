# 2026-11-16 记录结构整理

用户要求将长期规则、目标、当前交接和历史记录分开。只整理记录与检查入口，不运行CFD/MATLAB、修改物理算法或增加合格面积，不推送GitHub。

## 新入口

- [根执行约定](../../../AGENTS.md)：长期约定和恢复顺序。
- [项目目标](../../PROJECT_GOALS.md)：主路线、职责和验收；移除版本流水账及过期“下一步”。
- [当前交接](../../CURRENT_STATE.md)：唯一当前状态和行动入口。
- [维护规则](../../RECORDS_GUIDE.md)：何时更新哪个文件、并发写入和交接规则。
- [历史索引](../../EXECUTION_INDEX.md)：原文与版本报告/决定的定位，不复制长摘要。
- [只读检查程序](../../../scripts/checkProjectRecords.py)：大小、当前标记、引用、基础提交和历史原文SHA。不是科学验收。

## 完整保留的原文

| 时点 | 原AGENTS | 原目标 |
|---|---|---|
| 开始整理，科学基线v251 | [AGENTS.before.md](AGENTS.before.md)，310742字节 | [PROJECT_GOALS.before.md](PROJECT_GOALS.before.md)，68160字节 |
| 接入并发v252后的写入基线 | [AGENTS.pre_apply.md](AGENTS.pre_apply.md)，313013字节 | [PROJECT_GOALS.pre_apply.md](PROJECT_GOALS.pre_apply.md)，69020字节 |

均按原字节保存，包括全部用户决定、原错误/失败、存储清理记录和当时有效/已过时的指令。SHA、时间和源路径见 [migration_manifest.json](migration_manifest.json)。快照内部相对路径需按原源文件位置解释，不能把它们当作现行AGENTS或当前目标。较早历史和版本证据仍原位保存。

## 并发接续和语义纠正

最初读取时v252只有计划及未提交实现。准备新入口期间另一个工作流将v252结果写回两份共享文件；比较读取时SHA的检查失败，**没有覆盖新文件**。随后再次按原字节保存v252版本，读取其报告、NEXT_DECISION、最终归档回执及本地提交`78da62aa80f25593a9644c3d4eeb4075d438ae20`，才更新当前交接。

v252最终归档passed不代表原三链MATLAB父监督终态已知：原退出码/耗时/峰值及全程稳定仍未知。新交接保留该缺口、单站软件修复的独立回执、研究候选范围和新增合格面积0，不重跑补证据。所有v252源码/原报告/计划/原始回执由其原任务管理，本次不改。

旧目标第7节仍称v244是下一起点，而首部已推进到v251/v252；该冲突原文保留在快照，现行目标只指向CURRENT_STATE。具体关闭组合与长久规则分开，避免把旧CLI关闭误当永久禁令，也避免由历史计划自动重启计算。

[进程观察](process_observation.json)仅是定向Windows时点查询，不覆盖所有WSL进程/聊天编辑，不从“无匹配进程”推断原成功。并发工作区最初清单见 [git_status_before.txt](git_status_before.txt)。

## 交付核对

检查结果在`verification.json`生成；四份历史SHA、入口大小、当前v252指针、引用和基础提交均需通过后交付。检查并不重新验证科学结果或进程存活。新入口发布前再次核对最新源SHA，旧输入不匹配时停止该次覆盖、重新读入合并；不使用Git复原他人修改。

`AGENTS.next.md`与`PROJECT_GOALS.next.md`是此次准备的准确替换文本，保留用于审查。原输入和替换文本都在本目录，不要求重读整个历史才能继续工作。后续只维护实际入口，勿重新执行此次迁移或历史finalize脚本。

本目录的`.gitattributes`禁用文本换行转换，防止Git把原文CRLF归一化后改变SHA。提交前从Git暂存对象回读四份历史快照，核对字节数与SHA；结果另存`git_archive_verification.json`。该检查只证明归档字节保持，不授科学资格。
