function item=createRocketAxisymmetricWallCase(options)
%CREATEROCKETAXISYMMETRICWALLCASE 同一完整箭体的静止壁/内面黏性功受控工况。
% 模块一：完整沿用v68几何、飞行、气体、远场、Tw、SA输入和均匀Q初猜。
% 模块二：墙速度梯度和内面黏性功速度是两个独立旋钮，不导入旧流场。
% 模块三：新工况工厂单独标识，用新源码/新目录形成真实新链，不冒称v68续算。
% 本工况仍是零攻角无突起的数值候选；未验证稳态、物理精度或新增生产面积。
if nargin<1,options=struct;end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options), ...
    {'flight','geometry','domain','nuTildeRatio','wallVelocityGradient','viscousWorkVelocity'})), ...
    'RocketWallCase:Options','只接受原工况及两个显式离散对照选项。');
wallPolicy='dirichlet-tangent';workPolicy='legacy-wls';baseOptions=options;
if isfield(options,'wallVelocityGradient'),wallPolicy=options.wallVelocityGradient;baseOptions=rmfield(baseOptions,'wallVelocityGradient');end
if isfield(options,'viscousWorkVelocity'),workPolicy=options.viscousWorkVelocity;baseOptions=rmfield(baseOptions,'viscousWorkVelocity');end
assert(isText(wallPolicy)&&ismember(char(wallPolicy),{'legacy','dirichlet-tangent'})&& ...
    isText(workPolicy)&&ismember(char(workPolicy),{'legacy-wls','arithmetic'}), ...
    'RocketWallCase:Policy','墙梯度只允许legacy/dirichlet-tangent；内面功速度只允许legacy-wls/arithmetic。');
wallPolicy=char(wallPolicy);workPolicy=char(workPolicy);
item=createRocketAxisymmetricCase(baseOptions);
assert(isequal(item.boundary.jMin.velocity,[0 0]), ...
    'RocketWallCase:Boundary','本对照只支持显式常数静止墙速度。');
if ~strcmp(wallPolicy,'legacy'),item.boundary.jMin.velocityGradient=wallPolicy;end
if ~strcmp(workPolicy,'legacy-wls'),item.boundary.viscousWorkVelocity=workPolicy;end
item.caseFactory='createRocketAxisymmetricWallCase';
item.options=options;
item.inputContract.wallVelocityGradientPolicy=wallPolicy;
item.inputContract.viscousWorkVelocityPolicy=workPolicy;
item.inputContract.wallVelocityGradientBoundaryMeaning='static Cartesian wall velocity components have prescribed zero wall-tangential derivative';
item.inputContract.changedFromV68=['controlled wall gradient=' wallPolicy '; interior viscous-work velocity=' workPolicy ...
    '; same initial Q, geometry, gas, Tw, farfield and SA inputs'];
item.inputContract.wallGradientIsPhysicalAccuracyValidation=false;
end

function yes=isText(value)
yes=(ischar(value)&&isrow(value))||(isstring(value)&&isscalar(value)&&~ismissing(value));
end
