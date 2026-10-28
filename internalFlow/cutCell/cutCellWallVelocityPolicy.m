function policy=cutCellWallVelocityPolicy(boundary)
%CUTCELLWALLVELOCITYPOLICY 草稿：显式解析原曲面静止墙的速度切向梯度策略。
% 模块一：缺省/legacy保留旧路径，不静默消费未知表面导数。
% 模块二：dirichlet-tangent只支持已给定常数静止[0 0 0]，两条切向导数有确切来源。
% 在对应测试和完整NS集成通过前，本草稿不属于v82交付库或可用物理解。
assert(isstruct(boundary)&&isscalar(boundary),'CutNS:WallVelocityBoundary','墙边界必须为单一结构。');
policy='legacy';hasPolicy=isfield(boundary,'velocityGradient');
hasDerivative=any(isfield(boundary,{'velocityTangentialGradient','velocitySurfaceGradient'}));
if ~hasPolicy&&~hasDerivative,return;end
assert(isfield(boundary,'type')&&strcmp(boundary.type,'wall')&&hasPolicy, ...
    'CutNS:WallVelocityBoundary','速度梯度策略属于wall；导数字段不能单独出现。');
value=boundary.velocityGradient;
assert(((ischar(value)&&isrow(value))||(isstring(value)&&isscalar(value)&&~ismissing(value)))&& ...
    ismember(char(value),{'legacy','dirichlet-tangent'}),'CutNS:WallVelocityBoundary', ...
    'velocityGradient只支持legacy或dirichlet-tangent。');
assert(~hasDerivative,'CutNS:WallVelocityBoundary', ...
    '当前静止壁不消费额外速度表面导数；非零/函数壁速须另建完整导数合同。');
policy=char(value);if strcmp(policy,'legacy'),return;end
assert(isfield(boundary,'velocity')&&isnumeric(boundary.velocity)&&isreal(boundary.velocity)&& ...
    numel(boundary.velocity)==3&&all(isfinite(boundary.velocity(:)))&&all(boundary.velocity(:)==0), ...
    'CutNS:WallVelocityBoundary','当前dirichlet-tangent仅支持常数静止velocity=[0 0 0]。');
end
