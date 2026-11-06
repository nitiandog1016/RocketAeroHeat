# NACA TN3344 公开半球—圆柱热交换数据

选定Table I最后一项热交换工况：M1.97、Re_D2.28e6、总温570°R，直径4英寸。NASA报告：[原始来源](https://ntrs.nasa.gov/citations/19930084172)。

`heat_parameter.csv`为Figure 4(b)带中心点圆圈的20个手工图读点，量是局部Nu_s/√Re_s；不是直接测得的W/m²。`wall_temperature.csv`为Figure 3(a)上方加热壁三角符号22点，温度比使用同一工况总温换算。横坐标为从鼻尖起的母线弧长s/D，不是轴向距离。误差列只含图读误差，原文实验Nu约±15%等不确定度另见`contract.json`。

原图、像素中心、三点标定和第四网格点检查均保留。图读小数来自像素线性换算，不代表同位数实验精度。鼻尖温度点按图中轴线归为s/D=0，原像素映射偏差保留。热流从加热壁流向气体，与程序“入壁为正”约定相反。

下载原PDF的SHA及渲染图SHA见`source_manifest.json`；大PDF只在本地`outputs/v168_validation/sources`，未备份到GitHub。历史文件带436字节前缀，保留原字节并经恢复读取/Poppler渲染核对。图读准备和相容来流已经测试；层流轴线/Tw分布/原文归一化适配和数值求解尚未完成，物理资格为0。底面及支撑件不在实验热测量范围。
