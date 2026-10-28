function policy=cutCellWallThermalPolicy(boundary)
%CUTCELLWALLTHERMALPOLICY 显式解析曲面Dirichlet温度梯度策略，默认保留旧计算。
policy='legacy';hasPolicy=isfield(boundary,'thermalGradient');hasGradient=isfield(boundary,'temperatureSurfaceGradient');
if ~hasPolicy&&~hasGradient,return;end
assert(isfield(boundary,'type')&&strcmp(boundary.type,'wall')&&hasPolicy, ...
    'CutNS:WallThermalBoundary','切向温度数据只属于显式指定策略的wall边界。');
value=boundary.thermalGradient;
assert(((ischar(value)&&isrow(value))||(isstring(value)&&isscalar(value)&&~ismissing(value)))&& ...
    ismember(char(value),{'legacy','dirichlet-tangent'}), ...
    'CutNS:WallThermalBoundary','thermalGradient只支持legacy或dirichlet-tangent。');
policy=char(value);
if strcmp(policy,'legacy')
    assert(~hasGradient,'CutNS:WallThermalBoundary','legacy不消费所给切向导数，不能静默忽略。');
else
    assert(isfield(boundary,'temperature')&&(~isa(boundary.temperature,'function_handle')||hasGradient), ...
        'CutNS:WallThermalBoundary','函数墙温必须显式提供temperatureSurfaceGradient，单位K/m。');
end
end
