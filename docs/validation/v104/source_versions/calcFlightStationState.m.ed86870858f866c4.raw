function current = calcFlightStationState(mesh,condition,cfg)
%CALCFLIGHTSTATIONSTATE 从同一几何独立建立一个飞行时间站的外流与热边界。
% 模块作用：供初始壁温历程预计算和瞬态耦合共用，避免两条入口采用
% 不同的轴向波系/准流线流程。这里只建立给定工况的准定常快照，
% 不推进材料温度；壁温迭代由调用方在此快照上进行。
cfg.flight=condition;
% 每站基准统一采用初始壁温，不继承调用方旧的非均匀热流快照。
if isfield(cfg,'heating') && isfield(cfg.heating,'wallTemperatureByFace')
    cfg.heating=rmfield(cfg.heating,'wallTemperatureByFace');
end
fs=calcFreestream(condition);
current=modifiedNewtonianPressure(mesh,fs,cfg);
current=findStagnationRegion(current,fs,cfg);
current=solveExternalFlow(current,fs,cfg);
current=marchAxisymmetricShockExpansion(current,fs,cfg);
current=traceSurfaceStreamlines(current,cfg);
current=calcSurfaceHeating(current,fs,cfg);
end
