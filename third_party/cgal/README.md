# CGAL几何审核运行时

v136通过独立Python子进程使用CGAL 6.0.1.post202410241521 Windows cp311 wheel的Polygon Mesh Processing自交检查。源项目：https://github.com/CGAL/cgal-swig-bindings 。发布元数据：https://pypi.org/pypi/cgal/6.0.1.post202410241521/json 。wheel采用GPL-3.0-or-later/商业双许可，本项目记录和使用开源GPL路径；wheel内原许可随验证证据保存，未购买或声称商业许可。

本地Python3.11.9嵌入式运行时来自 https://www.python.org/ftp/python/3.11.9/python-3.11.9-embed-amd64.zip ，PSF许可原文保留。运行时位置`tools/runtime/cgal/6.0.1-py311/`被Git忽略，二进制未随项目上传。`scripts/installCGALGeometryRuntime.py`记录来源、发布SHA（wheel）、本地SHA、文件清单和真实导入状态，拒绝覆盖已有安装。

此Python包只测试`CGAL_Polygon_mesh_processing`、Polyhedron及Kernel的几何读取/自交路径；不声称所有wheel模块及其可选依赖均可用。v137新增的独立C++模块见下节，不属于此Python接口。

## v137 原生局部重网格

新增`local_tetra_remesh.cpp`，采用GPL-3.0-or-later，链接CGAL6.0.1头文件和Boost1.87.0；完整许可证原文见`docs/validation/v137/runtime/licenses/`。准备与构建脚本为`scripts/prepareCGALRemeshingRuntime.py`和`scripts/buildCGALRemeshingRuntime.py`，在本地`tools/runtime/cgal/6.0.1-native/`生成独立可执行文件，不修改系统环境或旧自交运行时。

只支持有效凸包三角剖分输入，原生导入与局部/边界冻结经两个小例验证。当前真实箭体带孔外域及既定局部选区不满足直接导入前提，不能用于当前候选或授予网格/物理资格。详情见`docs/CGAL_LOCAL_CONTRACT_V137.md`。安装脚本拒绝覆盖已有运行时；重现时应先核对现有清单，不删除仍被引用的文件。
