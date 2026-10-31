function scope=standaloneAirThermalScope(mode)
%STANDALONEAIRTHERMALSCOPE 独立研究内核的实际热物性范围说明。
% 模块作用：与生产侧壁外流的回退策略分开描述，避免独立内核拒绝
% 域外状态时却在来源元数据中写成“退回定比热”或“占位填值”。
if strcmp(mode,'thermallyPerfect')
    scope=['固定组分干空气NASA7，80~2500 K；80~150 K与150 K连续的' ...
        '振动冻结近似，150~200 K少量外推；若使用Sutherland，' ...
        '80~150 K黏度为工程外推；不含解离/电离；' ...
        '独立内核最终域外拒绝，不切换定比热、不占位填物理场'];
elseif strcmp(mode,'caloricallyPerfect')
    scope=['明确选用定比热理想气体cp=gamma*R/(gamma-1)，' ...
        '高温比热误差未修正；各模块另外检查自身适用温度及范围；' ...
        '不是NASA7失效后的静默回退'];
else
    error('RocketAeroHeat:AirThermal:Mode','未知独立内核热物性模式。');
end
end
