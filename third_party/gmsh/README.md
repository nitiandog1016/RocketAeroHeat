# 项目本地Gmsh运行时

固定官方Python发行包4.15.2，来源与完整SHA见`source_lock.json`，许可证及官方例外原文见`LICENSE.txt`。本仓库不包含DLL或上游源码副本；官方完整源码可从[gmsh.info](https://gmsh.info/)及[4.15.2源包](https://gmsh.info/src/gmsh-4.15.2-source.tgz)取得。

在项目根目录以已有Python运行`scripts/installGmshRuntime.py`。安装到被Git忽略的`tools/runtime/gmsh/4.15.2`，不修改系统PATH；重复运行核对已安装文件。计算仍由MATLAB管理，通过独立Python进程调用此网格器，不需要外部CFD场。

安装与调用测试采用Windows x64、Python3.7.7、NumPy1.21.6、SciPy1.7.3；没有验证其他平台。下载包约42.24MB，提取的Python/DLL/许可证约89.67MB，运行时容量不计入每次研究数据预算。已校验安装后删除下载副本，保留必要运行文件。库名或安装成功不授予体网格与气动热资格。
