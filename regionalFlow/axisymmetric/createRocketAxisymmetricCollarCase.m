function item=createRocketAxisymmetricCollarCase(options)
%CREATEROCKETAXISYMMETRICCOLLARCASE 完整箭体近壁法向层的独立五量SA研究工况。
% 模块一：沿用同一箭体/标准大气/物性/Tw/SA及均匀初猜，不导入旧末态或外部CFD。
% 模块二：默认normal-collar；可显式选同层数star-ray-reference做几何/流动对照。
% 模块三：壁速度梯度与内部功策略独立记录；默认legacy/arithmetic便于对照v69主链。
% 模块四：重新绑定真实体积和初积分，不把新几何、初猜或数值稳态晋升为生产资格。
if nargin<1,options=struct;end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options), ...
    {'flight','geometry','domain','nuTildeRatio','wallVelocityGradient','viscousWorkVelocity','gridPolicy','collar'})), ...
    'RocketCollarCase:Options','只接受已有工况、明确网格策略及近壁层配置。');
policy='normal-collar';baseOptions=options;
if isfield(options,'gridPolicy'),policy=options.gridPolicy;baseOptions=rmfield(baseOptions,'gridPolicy');end
assert(((ischar(policy)&&isrow(policy))||(isstring(policy)&&isscalar(policy)&&~ismissing(policy)))&& ...
    ismember(char(policy),{'normal-collar','star-ray-reference'}), ...
    'RocketCollarCase:Policy','网格策略只允许normal-collar或star-ray-reference。');policy=char(policy);
collar=struct;
if isfield(options,'collar'),collar=options.collar;baseOptions=rmfield(baseOptions,'collar');end
assert(isstruct(collar)&&isscalar(collar),'RocketCollarCase:Options','collar必须为标量结构体。');
assert(strcmp(policy,'normal-collar')||isempty(fieldnames(collar)), ...
    'RocketCollarCase:IgnoredOption','射线参考网格不得静默忽略近壁层配置。');
if ~isfield(baseOptions,'domain'),baseOptions.domain=struct;end
assert(isstruct(baseOptions.domain)&&isscalar(baseOptions.domain),'RocketCollarCase:Options','domain必须为标量结构体。');
if ~isfield(baseOptions.domain,'radialLayers'),baseOptions.domain.radialLayers=24;end
if ~isfield(baseOptions,'wallVelocityGradient'),baseOptions.wallVelocityGradient='legacy';end
if ~isfield(baseOptions,'viscousWorkVelocity'),baseOptions.viscousWorkVelocity='arithmetic';end
item=createRocketAxisymmetricWallCase(baseOptions);
if strcmp(policy,'normal-collar')
    [mesh,domain]=buildRocketAxisymmetricCollarDomain(item.surfaceMesh,item.geometryTruth,baseOptions.domain,collar);
    assert(mesh.cellCount==item.mesh.cellCount&&isequal(mesh.shape,item.mesh.shape), ...
        'RocketCollarCase:Topology','完整近壁层必须与同层数参考域具有同一逻辑控制体数量。');
    item.mesh=mesh;item.domain=domain;effectiveCollar=domain.collar.options;
else
    effectiveCollar=struct;
end
item.caseFactory='createRocketAxisymmetricCollarCase';item.options=options;
item.initialization.initialConservedIntegralPerRadian=sum(item.Q.*item.mesh.volumePerRadian,1);
item.initialization.geometryFactory=policy;
item.inputContract.gridPolicy=policy;item.inputContract.effectiveDomainOptions=item.domain.options;
item.inputContract.effectiveCollarOptions=effectiveCollar;
item.inputContract.firstOffsetIsOwnerWallDistance=false;
item.inputContract.changedFromV68=['complete original wall/farfield preserved; grid=' policy ...
    '; uniformly initialized same freestream five-component state on new control volumes; no imported solution'];
item.inputContract.wallGradientIsPhysicalAccuracyValidation=false;
end
