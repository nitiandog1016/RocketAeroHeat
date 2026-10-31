# 独立审核启动失败记录

第一次用`python -c`读取并compile执行`tests/auditActualRocketResult.py`时，启动包装没有提供`__file__`，在第7行计算项目根目录时发生`NameError: name '__file__' is not defined`，退出1。该失败发生在读取数值和写审核结果之前。

随后仅为执行命名空间提供真实脚本路径与`__name__='__main__'`，源码没有修改，重新执行独立审核退出0，生成`checks/independent_audit.json`。未重跑MATLAB、Gmsh或SU2；不将启动失败隐藏为首次通过。
