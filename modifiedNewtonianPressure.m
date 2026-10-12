function mesh = modifiedNewtonianPressure(mesh,freestream,cfg)
%MODIFIEDNEWTONIANPRESSURE 计算修正牛顿面元压力，结果存入 mesh.pressure。
% 输入：最终网格及 localGeometry、当前来流、可选 cfg。
% 迎风公式：Cp=CpMax*max(mu,0)^2，p=pinf+qinf*Cp。
% 背风及切向默认 Cp=0，即自由来流静压；这不是分离流/尾迹压力模型。
% 可设置 cfg.pressure.leewardCp 和 cfg.pressure.aftBaseCp 为常数。
% aftBaseCp 仅覆盖非迎风尾部平封盖；不覆盖反向流中的迎风尾盖。
% 方法主要用于高超声速迎风面工程估算；M>=1 只是公式定义域。
% 不处理凹形表面遮挡、激波干扰、粘性、真实气体或气动热。

if nargin < 3, cfg = struct; end
%% 1. 计算参考上限，重新分类以避免使用上一次攻角的分类结果
[CpMax,reference] = calcCpMax(freestream);
mesh = classifyAerothermalRegion(mesh,freestream,cfg);
r = mesh.aerothermalRegion;
leewardCp = 0; aftBaseCp = 0;
if isfield(cfg,'pressure')
    if isfield(cfg.pressure,'leewardCp'), leewardCp = cfg.pressure.leewardCp; end
    if isfield(cfg.pressure,'aftBaseCp'), aftBaseCp = cfg.pressure.aftBaseCp; end
end
validateattributes(leewardCp,{'numeric'},{'real','finite','scalar'});
validateattributes(aftBaseCp,{'numeric'},{'real','finite','scalar'});

%% 2. 使用原始 mu 而不是分类阈值，保留略微迎风面的非零压力增量
mu = r.mu; windward = mu > 0;
Cp = leewardCp*ones(size(mu));
Cp(windward) = CpMax*mu(windward).^2;
base = r.isAftCap & ~windward;
Cp(base) = aftBaseCp;
qinf = reference.dynamicPressure;
pinf = freestream.pressure;
pressure = pinf+qinf*Cp;
if any(~isfinite(pressure) | pressure < 0)
    error('RocketAeroHeat:Pressure:Nonphysical','指定压力系数产生负绝对压力或非有限值。');
end

%% 3. 保存压力、近似来源与计算时的来流，供 GUI 和后续积分读取
result = struct;
result.model = '修正牛顿法（定比热完全气体）';
result.Cp = Cp; result.CpMax = CpMax;
result.staticPressure = pressure; result.pressureIncrement = qinf*Cp;
result.freestreamPressure = pinf; result.dynamicPressure = qinf;
result.reference = reference; result.freestream = freestream;
result.mu = mu; result.isNewtonianFace = windward;
% 保存压力所对应的几何，供整体载荷模块识别过期压力。
result.geometry = struct('normal',mesh.localGeometry.normal, ...
    'center',mesh.localGeometry.center,'area',mesh.localGeometry.area);
result.isPrescribedBase = base; result.isPrescribedLeeward = ~windward & ~base;
result.leewardCp = leewardCp; result.aftBaseCp = aftBaseCp;
result.summary = struct('minimumCp',min(Cp),'maximumCp',max(Cp), ...
    'minimumPressure',min(pressure),'maximumPressure',max(pressure));
mesh.pressure = result;
% 旧载荷不能继续搭配新的压力；调用方需重新积分。
if isfield(mesh,'aerodynamicLoads'), mesh=rmfield(mesh,'aerodynamicLoads'); end
if isfield(mesh,'stagnation'), mesh=rmfield(mesh,'stagnation'); end
if isfield(mesh,'boundaryLayerEdge'), mesh=rmfield(mesh,'boundaryLayerEdge'); end
if isfield(mesh,'boundaryLayer'), mesh=rmfield(mesh,'boundaryLayer'); end
if isfield(mesh,'heating'), mesh=rmfield(mesh,'heating'); end
if isfield(mesh,'heatingCoverage'), mesh=rmfield(mesh,'heatingCoverage'); end
if isfield(mesh,'stagnationHeating'), mesh=rmfield(mesh,'stagnationHeating'); end
if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
if isfield(mesh,'surfaceStreamlines'), mesh=rmfield(mesh,'surfaceStreamlines'); end
if isfield(mesh,'axisymmetricMarch'), mesh=rmfield(mesh,'axisymmetricMarch'); end
if isfield(mesh,'externalFlow')
    saved=mesh.externalFlow.freestream;
    same=true;
    for key={'Mach','alpha','beta','altitude','gamma','R', ...
            'pressure','temperature','density','speed','totalPressure','totalTemperature'}
        name=key{1};
        same=same && isfield(saved,name) && ...
            isequal(saved.(name),freestream.(name));
    end
    if ~same, mesh=rmfield(mesh,'externalFlow'); end
end
end
