# v147 验证索引

结果与限制见[版本说明](../../SST_NATIVE_UPDATE_V147.md)，后续算法与停止条件见[算法决策](ALGORITHM_DECISION.md)。本目录包含固定源、事前计划、实际终态、失败和逐字执行来源；不以观测兼容授予物理资格。

- `PLAN.md`、`NUMERICAL_PREFLIGHT.md`：首轮运行前的问题、预算与门槛。
- `INDEX_CORRECTION.md`、`failure01/`：局部/全局编号错误及纠错前准确来源，原v1包/小例仍本地。
- `failure02/`：首次MATLAB附加字段拒绝及准确来源；只读修正，无数值重跑。
- `official_source/`、`source_basis.json`：固定官方更新、重编号和许可证。归档首次COPYING后缀查询不唯一，改用精确根路径，没有改变原源或数值结果。
- `candidate_source/`、`candidate_source_basis.json`：下一AUSMPLUSUP2候选的源码依据，尚未测试。
- `evidence/`：两次构建、五个小例、唯一真实更新、热审核/更新重放、两项MATLAB终态及独立导出证据。
- `source_versions/`：最终来源与运行期间的不同准确版本；导出复核器仅在MATLAB结束后改读workbench02。
- `delivery_manifest.json`：文件/散列、局部运行次数、存储口径以及未混入提交的既存UI依赖。

完整源文件在项目中直接替换，MATLAB注释说明模块用途。大CSV/MAT、既存源场和编译运行包未上传；仅靠克隆代码仓库不能打开本地已存真实场，需按合同准备相应原模型/初态。禁止覆盖已有运行记录或把历史失败重写为成功。
