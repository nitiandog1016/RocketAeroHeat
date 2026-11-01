# MMG局部四面体重网格

固定MMG5.8.0（LGPL-3.0-or-later），源码来自https://github.com/MmgTools/mmg/tree/v5.8.0 。独立入口`local_tetra_remesh.cpp`同许可；源码和许可证见`docs/validation/v138/official/`。上游数值内核未修改，CMake仅显式调用Windows已有sh。无SCOTCH/VTK/ELAS。

准备/构建：`scripts/prepareMMGLocalRuntime.py`、`scripts/buildMMGLocalRuntime.py`。运行时位于被Git忽略的`tools/runtime/mmg/5.8.0-local/`，保留源档、静态库与可重链接入口源；不覆盖已有安装。`repairMMGLocalWrapper.py`记录本轮首入口参数错误的单独重链接恢复，不是新安装的必需步骤。当前入口源已包含修正。

`mmgLocalInterface.py`冻结全部边界三角形和未选单元，按明确reference恢复固定节点，仅容许预先约定的缩放回写舍入；新自由节点保留原生坐标，禁止最近邻匹配。`testMMGLocalRemeshing.py`验证带孔域，`remeshSU2WithMMG.py`为本轮唯一真实候选。它们拒绝覆盖旧证据目录，重现实验须显式另选输出目录，不删除旧结果。

小例和真实读写/固定接口通过，但真实仍有201个低质四面体，未获计算资格。不得把MMG成功退出或此接口当成收敛气动热解。完整记录见`docs/MMG_LOCAL_REMESH_V138.md`。
