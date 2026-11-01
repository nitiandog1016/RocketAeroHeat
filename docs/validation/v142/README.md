# v142 验证索引

唯一TetGen已有网格重建小例失败，未启动真实候选。独立复核通过只确认失败事实，不授予质量或物理资格。

- [运行前计划和固定参数](plan.md)
- [源码下载来源](runtime/source.json)、[构建与二进制身份](runtime/runtime.json)、[有限源码合同](runtime/source_contract.json)
- [小例失败](hollow_test/report.json)、[原生终态](hollow_test/native/process/process.json)、[日志](hollow_test/native/process/stdout.log)
- [独立只读核对](hollow_test/independent_review.json)
- [准确源/存储清单](storage_and_sources.json)
- [方法、限制和下一起点](../../TETGEN_FIXED_REPAIR_V142.md)

`source_versions/`保存实际执行字节；本轮未修正/重跑数值。运行时状态记录是构建时事实，后续小例失败见本目录，不能将built解释为可用修复器。完整官方源按固定URL及9文件Git blob可取得；原许可在`runtime/upstream_LICENSE`。大运行时未上传，v141真实输入仍只在本地并保持不变。
