function setAeroHeatStage(fig,stage)
%SETAEROHEATSTAGE 更新统一框架的阶段说明，不改变计算状态。
if isgraphics(fig) && isstruct(fig.UserData) && isfield(fig.UserData,'stageLabel')
    labels = {'【1 几何朝向与网格】 → 2 飞行条件 → 3 压力与气动热', ...
        '1 几何与网格已确认 → 【2 飞行条件】 → 3 压力与气动热', ...
        '1 几何与网格已确认 → 2 飞行条件已确认 → 【3 压力与气动热】'};
    fig.UserData.stageLabel.Text = labels{stage};
    drawnow;
end
end
