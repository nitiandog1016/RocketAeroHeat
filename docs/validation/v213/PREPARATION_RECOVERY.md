# 原生准备的包布局修复

prep_sim01在读取prepsim.lua前失败：官方prepsim入口按实际可执行文件的`../lib/prepsim.lua`定位，不只使用DGD环境。原binary已构建通过，未执行流动；首失败保留。

将同一binary以硬链接放到独立包bin/probe，包lib只复制准确官方prepsim.lua；不复制大二进制、不改共享安装，不重建。唯一prep_sim02仍在原两次准备中prep-sim累计30秒预算内，自动追加0，输入包括原失败回执及准确Lua SHA。该布局恢复不改变模型或数值门槛。
