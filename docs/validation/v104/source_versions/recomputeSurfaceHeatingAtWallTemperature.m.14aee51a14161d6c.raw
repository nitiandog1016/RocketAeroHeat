function updated = recomputeSurfaceHeatingAtWallTemperature(mesh,freestream,cfg,temperature)
%RECOMPUTESURFACEHEATINGATWALLTEMPERATURE 给定非均匀壁温，重算整条热边界。
% 模块一：校验面元顺序与壁温；不改变外流、几何和准流线。
% 模块二：由每个面元壁温重算 Eckert 参考黏度、Thwaites/Head 沿程
% 状态、分离筛选、Reynolds 比拟和球鼻局部焓差。上游壁温可以影响
% 下游热流，即使下游自身壁温没有改变；不是逐点乘固定幂指数。
% 模块三：清除旧壁温历程，保存算法和输入壁温场。
% 这是给定壁温的准定常热边界评估接口，尚不独立推进材料温度，
% 也不重算粘性-无粘耦合外流；未来隐式热响应迭代可反复调用它。

if ~isfield(mesh,'surfaceStreamlines') || ~isfield(mesh,'externalFlow')
    error('RocketAeroHeat:WallField:ExternalState', ...
        '沿程壁温重算需要当前外流快照和三维准流线路径。');
end
cfg.heating.wallTemperatureByFace=temperature;
temperature=surfaceWallTemperature(size(mesh.F,1),cfg);
updated=calcSurfaceHeating(mesh,freestream,cfg);
updated.heating.wallTemperatureByFace=temperature;
updated.heating.wallCouplingMethod= ...
    '给定逐面元壁温，重积分 Thwaites/Head 状态并重判适用域';
updated.heating.wallCouplingScope= ...
    '固定几何、外流和准流线的准定常热边界评估；不等于瞬态材料温度求解';
end
