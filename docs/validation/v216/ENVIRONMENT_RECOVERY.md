# v216离线环境恢复

首inventory因当前WSL `/usr/bin/python3` 无SciPy在导入阶段失败（监督0.554656秒），未读取母网格。失败保留，不安装新依赖。改用已有 `D:/anacondapython/python.exe`，原SciPy/NumPy直接复用。

事前调整离线监督合同：Windows用psutil对本任务进程树同时采样RSS/虚拟字节，超过4GiB即停止；1CPU、累计120秒、128MiB输出不变。没有硬RLIMIT_AS，不把采样上限冒称硬地址空间保护。只读Python设计不启动原生/网格/CFD，不将该合同转授原生求解。

inventory2成功读取母数据，但摘要器错误要求不存在的pyramids字段而KeyError（.812秒）；准确原v144只有棱柱/四面体，缺字段按空集合处理，与既有auditSU2HybridTopology一致。仅恢复元数据报告，不改变单元类型/连接或重新网格化，累计预算不增加。

design01的NumPy旧版本对Fortran连续数组直接view结构键出现警告并导致边界匹配断言失败（1.891秒）；原整数面计数已通过，但没有生成分区结果。仅对规范键显式ascontiguousarray，再在design02执行同一预定选择规则；不改变原坐标、连接、分区宽度或门槛，首准确任务源码/日志保留。

design02越过面配对，继而暴露“完整底面节点编号可直接取模”的错误（2.031秒）。v133棱边采用两侧法向，底面有额外列节点；v144以sourcePreparedNodeIDs保留其祖先编号。design03改为按准确祖先映射及已绑定fan报告的baseNodeIDs回原列，再检查每面32层及全部2852尾缘扇。仍为同一预定分区，不移动点、不删尾缘、不改变物理门槛；失败不追改。
