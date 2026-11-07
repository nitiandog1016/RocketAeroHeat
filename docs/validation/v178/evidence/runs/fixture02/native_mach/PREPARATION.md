# Eilmer 完整域输入：仅准备，原生未验证

固定 Eilmer v5.0.0 / f53f4609a0331d48efee69a4e4f3c3598378cc03。
本目录不包含求解结果；不表示原生准备/收敛/热流/实验验证已通过。
引用原 SU2 文件，禁止改动；没有复制大网格。实际 WSL 路径尚须验证。
原生读取必须保留 expect_gmsh_order_for_wedges=true。边界标签和全部原壁面保留。
使用 Wilcox k-omega，不能称为 SST，不能接续旧 SST 保存场。
官方伪空气热力表 200–6000 K；表外采用端点定cp延拓h/s，不裁温，表外物理资格未知。
不能把有限值或数值收敛当成全面积物理求解资格。RANS近壁/壁热和分离区须独立验证。
loads 的 q_total=q_cond+q_diff 以内部面法向定义，必须连同outsign核对入壁正号；
原生非结构y+使用L_min/2估距，不能直接充当项目严格法向距离验收。

安装/编译与系统变更另行处理，本准备器没有启动它们。固定运行时就绪且输入SHA复核后，
官方命令顺序是 prep-gas gas.inp rah-air.gas、lmr prep-grids --job=grid.lua、lmr prep-sim --job=job.lua，
然后才可在已登记预算与监督下执行 lmr run。本轮这些命令全部未执行。
原始面顺序的结果映射、q符号/面积积分与稀疏输出读取均待独立实现，不能直接装载研究结果。
