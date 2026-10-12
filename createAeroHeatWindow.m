function [fig,host] = createAeroHeatWindow()
%CREATEAEROHEATWINDOW 创建统一应用框架，不执行网格或物理计算。
% 顶部阶段指示保持可见，内容区依次容纳网格、飞行条件和压力结果。
% 通过 setAeroHeatStage 更新阶段；宿主面板不随页面切换而销毁。
fig = uifigure('Name','RocketAeroHeat · 气动与气动热分析', ...
    'Position',[30 30 1580 940],'Color',[0.97 0.97 0.97]);
layout = uigridlayout(fig,[2 1]);
layout.RowHeight = {46,'1x'}; layout.Padding = [6 6 6 6];
steps = uilabel(layout,'Text','1  几何朝向与网格   →   2  飞行条件   →   3  压力与气动热', ...
    'FontSize',17,'FontWeight','bold');
host = uipanel(layout,'BorderType','none');
fig.UserData = struct('stageLabel',steps,'host',host);
end
