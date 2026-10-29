# 共同三维激波入口：小型输入与冻结证据

- `actual_rocket_geometry_and_freestream.mat`仅含原箭体V/F/面积及M5、30km、20°均匀自由来流。无旧压力、热流或边界层外缘；用它即可独立重算新入口与几何接口。原大文件的值/来源绑定见`input_binding.json`。
- `shock_surface_inlet_evidence.mat`为共同候选激波的定比热/NASA7状态、导数、平面激波和几何验证原MAT，373470字节。给定候选形状没有通过箭体匹配，不能解释为真实激波层或外缘。
- `shock_body_geometry_evidence.mat`保存候选激波内法线与真实原三角面的最近几何交点。不赋壁面压力、温度或热流。

完整运行命令、原任务终态、独立审核及限制见`docs/validation/v97/README.md`。本组输入没有大型矩阵或完整工程副本；不声称重新生成MAT容器必然得到历史SHA。
