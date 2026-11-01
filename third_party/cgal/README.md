# CGAL几何审核运行时

v136通过独立Python子进程使用CGAL 6.0.1.post202410241521 Windows cp311 wheel的Polygon Mesh Processing自交检查。源项目：https://github.com/CGAL/cgal-swig-bindings 。发布元数据：https://pypi.org/pypi/cgal/6.0.1.post202410241521/json 。wheel采用GPL-3.0-or-later/商业双许可，本项目记录和使用开源GPL路径；wheel内原许可随验证证据保存，未购买或声称商业许可。

本地Python3.11.9嵌入式运行时来自 https://www.python.org/ftp/python/3.11.9/python-3.11.9-embed-amd64.zip ，PSF许可原文保留。运行时位置`tools/runtime/cgal/6.0.1-py311/`被Git忽略，二进制未随项目上传。`scripts/installCGALGeometryRuntime.py`记录来源、发布SHA（wheel）、本地SHA、文件清单和真实导入状态，拒绝覆盖已有安装。

当前只测试`CGAL_Polygon_mesh_processing`、Polyhedron及Kernel的几何读取/自交路径；不声称所有wheel模块及其可选依赖均可用。CGAL四面体重网格的C++集成仍是后续工作，不包含在此运行时接口的已测能力内。
