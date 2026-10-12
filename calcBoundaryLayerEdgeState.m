function mesh = calcBoundaryLayerEdgeState(mesh,freestream,cfg)
%CALCBOUNDARYLAYEREDGESTATE 由面元压力估算边界层外缘状态。
% 本模块是工程闭合近似，不是激波/外流场求解器。
% 算法：沿用修正牛顿法的壁面压力，并假设外缘与自由来流具有相同
% 总温、总压；按定比热完全气体等熵关系反求 M_e、T_e、rho_e、U_e。
% 激波熵增会降低真实总压，因此本近似尤其不适用于鼻尖及强压缩区。

%% 1. 验证压力与当前飞行状态一致，避免用旧攻角结果计算新热流
if nargin<3, cfg=struct; end %#ok<NASGU>
% 优先使用经过网格、工况和状态方程校验的真实外流结果。
if isfield(mesh,'externalFlow')
    imported=mesh.externalFlow;
    if ~isfield(mesh,'localGeometry') || ...
            ~isequal(imported.V,mesh.V) || ~isequal(imported.F,mesh.F) || ...
            ~isequal(imported.faceCenter,mesh.localGeometry.center) || ...
            ~isequal(imported.faceNormal,mesh.localGeometry.normal)
        error('RocketAeroHeat:Edge:StaleExternalFlow','外流结果对应的网格已改变。');
    end
    for key={'Mach','alpha','beta','altitude','gamma','R', ...
            'pressure','temperature','density','speed','totalPressure','totalTemperature'}
        name=key{1};
        if ~isequal(imported.freestream.(name),freestream.(name))
            error('RocketAeroHeat:Edge:StaleExternalFlow','外流结果对应的工况已改变：%s。',name);
        end
    end
    mesh.boundaryLayerEdge=imported.edge;
    if isfield(mesh,'heating'), mesh=rmfield(mesh,'heating'); end
    if isfield(mesh,'heatingCoverage'), mesh=rmfield(mesh,'heatingCoverage'); end
    if isfield(mesh,'boundaryLayer'), mesh=rmfield(mesh,'boundaryLayer'); end
    if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
    if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
    if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
    return;
end
if ~isfield(mesh,'pressure') || ~isfield(mesh,'aerothermalRegion')
    error('RocketAeroHeat:Edge:Input','请先执行 modifiedNewtonianPressure。');
end
need={'Mach','alpha','beta','pressure','temperature','gamma','R','speed', ...
    'totalPressure','totalTemperature'};
for j=1:numel(need)
    key=need{j};
    if ~isfield(freestream,key) || ~isfield(mesh.pressure.freestream,key) || ...
            ~isequal(freestream.(key),mesh.pressure.freestream.(key))
        error('RocketAeroHeat:Edge:StalePressure','压力与当前 freestream 不一致：%s。',key);
    end
end
p=mesh.pressure.staticPressure(:); nf=size(mesh.F,1);
validateattributes(p,{'numeric'},{'real','finite','positive','size',[nf 1]});
g=freestream.gamma; R=freestream.R; p0=freestream.totalPressure;
if g<=1 || p0<=0 || freestream.Mach<1
    error('RocketAeroHeat:Edge:Domain','要求 M>=1、gamma>1 和正总压。');
end

%% 2. 从压力反求外缘马赫数和热力状态
% p>p0 时不存在等熵解；这些面元不伪造外缘结果。
valid=p<=p0*(1+1e-10);
M=nan(nf,1); T=M; rho=M; U=M; a=M; viscosity=M;
ratio=max(1,p0./p(valid));
M(valid)=sqrt(2/(g-1)*expm1((g-1)/g*log(ratio)));
T(valid)=freestream.totalTemperature./(1+0.5*(g-1)*M(valid).^2);
rho(valid)=p(valid)./(R*T(valid));
a(valid)=sqrt(g*R*T(valid)); U(valid)=M(valid).*a(valid);
% 标准空气 Sutherland 公式；与 standardAtmosphere 的输运模型一致。
viscosity(valid)=1.458e-6*T(valid).^(3/2)./(T(valid)+110.4);

%% 3. 保存算法来源、适用性和压力快照
edge=struct;
edge.method='修正牛顿壁压 + 自由来流总状态等熵反算（工程近似）';
edge.pressureSource=mesh.pressure.model;
edge.assumption='假设外缘总压等于自由来流总压；未计局部激波总压损失';
edge.pressure=p; edge.totalPressure=p0; edge.totalTemperature=freestream.totalTemperature;
edge.Mach=M; edge.temperature=T; edge.density=rho; edge.speed=U;
edge.speedOfSound=a; edge.dynamicViscosity=viscosity;
edge.algebraicallyValid=valid;
edge.physicallyUnresolved=false(nf,1);
if isfield(mesh,'stagnation')
    edge.physicallyUnresolved=mesh.stagnation.candidateFace;
end
edge.freestream=freestream;
edge.status='仅供初步边界层估算；驻点/强压缩区须由外流求解替换';
mesh.boundaryLayerEdge=edge;
if isfield(mesh,'heating'), mesh=rmfield(mesh,'heating'); end
if isfield(mesh,'heatingCoverage'), mesh=rmfield(mesh,'heatingCoverage'); end
if isfield(mesh,'boundaryLayer'), mesh=rmfield(mesh,'boundaryLayer'); end
if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
end
