function mesh = calcSurfaceHeating(mesh,freestream,cfg)
%CALCSURFACEHEATING 计算附着流侧壁的工程对流热流密度。
% 算法：Eckert 参考温度法估算输运性质；有三维准流线时，层流
% 使用沿程 Thwaites 型动量积分及局部速度梯度；已转捩湍流使用
% von Karman 动量积分 + Ludwieg-Tillmann Cf。没有三维准流线时
% 湍流使用局部平板 Cf_t=0.0592/Re_s*^0.2，层流也退回局部平板
% Cf_l=0.664/sqrt(Re_s*)，并明确标记为低阶近似。
% Reynolds 比拟 St=Cf/(2 Pr^(2/3))；恢复温度使用 r=sqrt(Pr)
% （层流）或 r=Pr^(1/3)（湍流）。默认干空气以焓方程求 Taw，
% q=St*rho_e*U_e*[h(Taw)-h(Tw)]；也可显式选旧版定比热。
% 以 Re_s=5e5 为默认经验转捩阈值；它不是稳定性预测。
% 附着流侧壁以外的球形鼻端可用 calcBluntNoseHeating 低阶延拓补齐；
% 其余封盖、坏拓扑及无有效状态的面元仍返回 NaN。
% 内置膨胀模型给出附着背风外流时，该区也参与工程热流计算。
% q>0 表示气流加热壁面；本模块不计算辐射、烧蚀或激波干扰。
% 连续介质 Knudsen 筛选在热流之前执行；滑移/过渡流面元及同一
% 准流线下游保持无解，不将无滑移关联式用于其适用域之外。
% 可选 cfg.heating.wallTemperatureByFace 按面元给出非均匀壁温；
% 沿程状态使用该场重新积分，未知壁温保持无解。

%% 1. 检查输入，并按当前工况更新外缘状态
if nargin<3, cfg=struct; end
if ~isfield(cfg,'flight') || ~isfield(cfg.flight,'Twall')
    error('RocketAeroHeat:Heating:WallTemperature','需要 cfg.flight.Twall [K]。');
end
Tw=surfaceWallTemperature(size(mesh.F,1),cfg);
Pr=0.72; criticalRe=5e5;
if isfield(cfg,'heating')
    if isfield(cfg.heating,'Prandtl'), Pr=cfg.heating.Prandtl; end
    if isfield(cfg.heating,'transitionRe'), criticalRe=cfg.heating.transitionRe; end
end
validateattributes(Pr,{'numeric'},{'real','finite','scalar','>',0,'<',1});
validateattributes(criticalRe,{'numeric'},{'real','finite','scalar','positive'});
if ~isfield(mesh,'stagnation'), mesh=findStagnationRegion(mesh,freestream,cfg); end
mesh=calcStagnationHeating(mesh,freestream,cfg);
mesh=calcBoundaryLayerEdgeState(mesh,freestream,cfg);
mesh=assessContinuumRegime(mesh,freestream,cfg);
e=mesh.boundaryLayerEdge; reg=mesh.aerothermalRegion;
continuum=mesh.continuumRegime;
if ~isfield(mesh,'surfaceCoord') || ~isfield(mesh.surfaceCoord,'faceS')
    error('RocketAeroHeat:Heating:SurfaceDistance','缺少鼻端起算的面元子午弧长。');
end
nf=size(mesh.F,1);
if isfield(mesh,'surfaceStreamlines')
    path=mesh.surfaceStreamlines;
    if ~strcmp(path.velocitySource,e.method) || numel(path.distance)~=nf
        error('RocketAeroHeat:Heating:StaleStreamlines', ...
            '三维流线与当前外缘状态不一致，请重新追踪。');
    end
    s=path.distance(:);
    distanceMethod=path.method;
else
    s=mesh.surfaceCoord.faceS(:);
    distanceMethod='使用鼻端起算子午弧长代替真实三维流线长度';
end
if numel(s)~=nf || any(s(isfinite(s))<0)
    error('RocketAeroHeat:Heating:SurfaceDistance','面元流动距离无效。');
end
[airMode,airScope]=selectAirThermalModel(freestream,cfg);
cpConstant=freestream.gamma*freestream.R/(freestream.gamma-1);

%% 2. 逐面元标记适用性；不能计算的面元保持 NaN
% 近切向面元允许参与：圆柱零攻角时 mu 的正负会受三角面法向舍入影响。
allowed=e.algebraicallyValid & reg.isSidewall & ...
    ~reg.isNoseTipAdjacent & ~reg.isAftTipAdjacent & ...
    ~reg.featureFace & ~reg.topologyIssue & ...
    ~mesh.stagnation.candidateFace & s>0 & e.speed>0;
% 最前端 1% 长度的尖鼻区没有可确认的流线起点，故不套用平板公式。
L=max(mesh.V(:,1))-min(mesh.V(:,1));
nearTip=s<0.01*L;
allowed=allowed & ~nearTip;

%% 3. Eckert 参考温度与局部 Reynolds 数
Te=e.temperature; U=e.speed;
recoverLam=sqrt(Pr); recoverTurb=Pr^(1/3);
if strcmp(airMode,'thermallyPerfect')
    % 使用 h(Taw)=h(Te)+r*U²/2，再把焓差换成壁温求解器所需的等效 cp。
    TawLam=airRecoveryTemperature(Te,U,recoverLam);
    TawTurb=airRecoveryTemperature(Te,U,recoverTurb);
    cpHeatLam=airEnthalpySecantCp(TawLam,Tw);
    cpHeatTurb=airEnthalpySecantCp(TawTurb,Tw);
else
    TawLam=Te+recoverLam*U.^2/(2*cpConstant);
    TawTurb=Te+recoverTurb*U.^2/(2*cpConstant);
    cpHeatLam=cpConstant*ones(nf,1);
    cpHeatTurb=cpHeatLam;
end
thermalPropertyUnavailableFace=allowed & ...
    (~isfinite(Tw) | ~isfinite(TawLam) | ~isfinite(TawTurb) | ...
    ~isfinite(cpHeatLam) | ~isfinite(cpHeatTurb));
if strcmp(airMode,'thermallyPerfect') && ...
        isfield(e,'thermalExpansionOutOfRangeFace')
    % 低于 80 K 的定比热外缘只用于维持离散准流线连通；该面元
    % 没有 NASA7 热力解，不能因占位温度偶然落在范围内而计算热流。
    thermalPropertyUnavailableFace=thermalPropertyUnavailableFace | ...
        (allowed & logical(e.thermalExpansionOutOfRangeFace(:)));
end
% 热物性空白不应重置流线上的已转捩状态：转捩仍由外流 Reynolds 数
% 沿完整可用路径判定，最后才筛去热物性越界的热流面元。
% 连续介质失效已沿路径传播，其面元不能再进入转捩 Reynolds 判据。
flowAllowed=allowed & ~continuum.unresolvedFace;
% 参考温度 T*=Te+0.5(Tw-Te)+0.22(Taw-Te)。
TrefLam=Te+0.5*(Tw-Te)+0.22*(TawLam-Te);
TrefTurb=Te+0.5*(Tw-Te)+0.22*(TawTurb-Te);
rhoRefLam=e.pressure./(freestream.R*TrefLam);
rhoRefTurb=e.pressure./(freestream.R*TrefTurb);
muRefLam=1.458e-6*TrefLam.^(3/2)./(TrefLam+110.4);
muRefTurb=1.458e-6*TrefTurb.^(3/2)./(TrefTurb+110.4);
nuRefTurb=muRefTurb./rhoRefTurb;
ReLam=rhoRefLam.*U.*s./muRefLam;
ReTurb=rhoRefTurb.*U.*s./muRefTurb;
ReEdge=e.density.*U.*s./e.dynamicViscosity;
hasIntegratedLayer=isfield(mesh,'surfaceStreamlines');
separationBlocked=false(nf,1);
separationTrace=[];
if hasIntegratedLayer
    mesh=calcBoundaryLayerState(mesh,freestream,cfg);
    boundaryLayer=mesh.boundaryLayer;
    % 先沿流线判转捩：下游即使 Re_s 暂时降低，也不能自动复位层流。
    firstPass=calcTransition(ReEdge,flowAllowed,cfg,mesh.surfaceStreamlines);
    % 仅未转捩层流的 Thwaites 分离才构成起点；当前没有再附模型，
    % 故该路径的全部下游面元均不再套用附着层流或平板湍流热流公式。
    separationOnset=flowAllowed & ...
        boundaryLayer.possibleLaminarSeparationFace & ...
        ~firstPass.turbulentFace;
    separationTrace=propagateSurfaceSeparation( ...
        mesh.surfaceStreamlines,separationOnset);
    separationBlocked=separationTrace.unresolvedFace;
    % 层流积分无解的面元保持 NaN；已转捩面仍用湍流工程模型。
    allowed=flowAllowed & ~thermalPropertyUnavailableFace & ...
        ~continuum.unresolvedFace & ...
        ~separationBlocked & ...
        (firstPass.turbulentFace | boundaryLayer.validFace);
    % 再筛选热流输出，但仍在完整有效外流路径上维持已转捩状态。
    % 热物性暂时无解的面元不能把下游误判回层流。
    transition=calcTransition(ReEdge,flowAllowed,cfg, ...
        mesh.surfaceStreamlines,allowed);
else
    allowed=flowAllowed & ~thermalPropertyUnavailableFace & ...
        ~continuum.unresolvedFace;
    transition=calcTransition(ReEdge,allowed,cfg);
end
if hasIntegratedLayer
    turbulentBoundaryLayer=calcTurbulentBoundaryLayerState( ...
        mesh,nuRefTurb,transition,cfg);
    % Head 疑似分离/积分失败已穿过所有下游面元传播；不能再由
    % 平板热流初始化或鼻端延拓补值，也不能重新判为附着层流。
    turbulentBlocked=turbulentBoundaryLayer.separationUnresolvedFace | ...
        turbulentBoundaryLayer.numericalUnresolvedFace;
    allowed=allowed & ~turbulentBlocked;
    separationBlocked=separationBlocked | turbulentBlocked;
    transition=calcTransition(ReEdge,flowAllowed,cfg, ...
        mesh.surfaceStreamlines,allowed);
else
    turbulentBoundaryLayer=[];
end
lam=nan(nf,1); turb=lam; selected=lam;
heatCoefficientLam=nan(nf,1); heatCoefficientTurb=heatCoefficientLam;
if hasIntegratedLayer
    laminarEligible=transition.laminarFace;
    heatCoefficientLam(laminarEligible)= ...
        boundaryLayer.skinFrictionCoefficient(laminarEligible) ./ ...
        (2*Pr^(2/3)).*e.density(laminarEligible).*U(laminarEligible).* ...
        cpHeatLam(laminarEligible);
else
    heatCoefficientLam(allowed)=0.664./sqrt(ReLam(allowed))./(2*Pr^(2/3)) .* ...
        e.density(allowed).*U(allowed).*cpHeatLam(allowed);
end
heatCoefficientTurb(allowed)=0.0592./ReTurb(allowed).^0.2./(2*Pr^(2/3)) .* ...
    e.density(allowed).*U(allowed).*cpHeatTurb(allowed);
if hasIntegratedLayer
    turbulentEligible=transition.turbulentFace & ...
        turbulentBoundaryLayer.validFace;
    heatCoefficientTurb(turbulentEligible)= ...
        turbulentBoundaryLayer.skinFrictionCoefficient(turbulentEligible) ./ ...
        (2*Pr^(2/3)).*e.density(turbulentEligible).* ...
        U(turbulentEligible).*cpHeatTurb(turbulentEligible);
end
if hasIntegratedLayer
    lam(laminarEligible)=heatCoefficientLam(laminarEligible).* ...
        (TawLam(laminarEligible)-Tw(laminarEligible));
else
    lam(allowed)=heatCoefficientLam(allowed).*(TawLam(allowed)-Tw(allowed));
end
turb(allowed)=heatCoefficientTurb(allowed).*(TawTurb(allowed)-Tw(allowed));
turbulent=transition.turbulentFace;
laminar=transition.laminarFace;
selected(laminar)=lam(laminar); selected(turbulent)=turb(turbulent);
selectedCoefficient=nan(nf,1); selectedRecovery=nan(nf,1);
selectedCoefficient(laminar)=heatCoefficientLam(laminar);
selectedCoefficient(turbulent)=heatCoefficientTurb(turbulent);
selectedRecovery(laminar)=TawLam(laminar);
selectedRecovery(turbulent)=TawTurb(turbulent);

%% 4. 用已验证球鼻的驻点锚定值填补侧壁模型未覆盖的鼻端
% 保留已有侧壁解；鼻端余弦延拓单独编号，不能误标成 Thwaites 层流。
nose=calcBluntNoseHeating(mesh,freestream,Tw);
noseFill=nose.faceMask & ~allowed & ~separationBlocked & ...
    ~continuum.unresolvedFace;
selected(noseFill)=nose.heatFlux(noseFill);
selectedCoefficient(noseFill)=nose.heatTransferCoefficient(noseFill);
selectedRecovery(noseFill)=nose.recoveryTemperature(noseFill);
faceMethodID=zeros(nf,1,'uint8');
faceMethodID(laminar)=1;
faceMethodID(turbulent)=2;
faceMethodID(noseFill)=3;
detachedHeat=false(nf,1);
if isfield(e,'faceMethodID')
    detachedHeat=(e.faceMethodID(:)==5) & (laminar|turbulent);
    faceMethodID(detachedHeat & laminar)=4;
    faceMethodID(detachedHeat & turbulent)=5;
end
calculated=allowed|noseFill;

%% 5. 保存分项结果与方法说明，供 GUI、报告和后续模块读取
h=struct;
if hasIntegratedLayer
    h.method=['Eckert 参考温度 + 沿流线 Thwaites 型层流积分/' ...
        '湍流积分 + Reynolds 比拟'];
    h.boundaryLayerMethod=boundaryLayer.method;
    h.boundaryLayerValidity=boundaryLayer.scope;
    h.turbulentMethod=turbulentBoundaryLayer.method;
    h.turbulentScope=turbulentBoundaryLayer.scope;
    h.turbulentBoundaryLayer=turbulentBoundaryLayer;
else
    h.method='Eckert 参考温度 + 局部平板 Cf + Reynolds 比拟';
    h.boundaryLayerMethod='无三维流线：局部平板层流近似';
    h.boundaryLayerValidity='忽略沿程压力梯度';
    h.turbulentMethod='无三维流线：局部平板 Cf_t=0.0592/Re_s^0.2';
    h.turbulentScope='忽略沿程压力梯度的附着湍流低阶估算';
    h.turbulentBoundaryLayer=[];
end
if strcmp(airMode,'thermallyPerfect')
    h.method=[h.method ' + 干空气 NASA7/低温振动冻结焓恢复与等效比热'];
else
    h.method=[h.method ' + 定比热恢复温度'];
end
h.edgeMethod=e.method;
h.thermalPropertyMethod=airMode;
h.thermalPropertyScope=airScope;
h.thermalPropertyUnavailableFace=thermalPropertyUnavailableFace;
h.wallTemperatureUnavailableFace=~isfinite(Tw);
if strcmp(airMode,'thermallyPerfect')
    h.wallTemperatureUnavailableFace=h.wallTemperatureUnavailableFace | Tw<80 | Tw>2500;
end
h.continuumMethod=continuum.method;
h.continuumScope=continuum.scope;
h.continuumUnresolvedFace=continuum.unresolvedFace;
h.transitionMethod=transition.method;
h.stagnationMethod=mesh.stagnationHeating.method;
h.stagnationStatus=mesh.stagnationHeating.status;
h.noseMethod=nose.method; h.noseScope=nose.scope;
h.noseDistribution=nose; h.noseFace=noseFill;
h.faceMethodID=faceMethodID;
h.faceMethodNames={'未计算','沿流线层流','沿流线湍流/显式平板回退', ...
    '球形鼻端偏置驻点余弦延拓', ...
    '脱体闭合外流 + 沿流线层流', ...
    '脱体闭合外流 + 沿流线湍流/平板回退'};
h.detachedClosureHeatFace=detachedHeat;
h.separationUnresolvedFace=false(nf,1);
if ~isempty(separationTrace)
    h.separationUnresolvedFace=separationTrace.unresolvedFace;
end
h.turbulentSeparationUnresolvedFace=false(nf,1);
h.turbulentNumericalUnresolvedFace=false(nf,1);
if hasIntegratedLayer
    h.turbulentSeparationUnresolvedFace=turbulentBoundaryLayer.separationUnresolvedFace;
    h.turbulentNumericalUnresolvedFace=turbulentBoundaryLayer.numericalUnresolvedFace;
end
h.separationTrace=separationTrace;
h.turbulentAssumptionMethod=['已转捩有效区采用 ' h.turbulentMethod ...
    '；未转捩侧壁仅给出局部平板湍流假设对照；已分离/积分失败路径保持 NaN'];
% 标量字段保留为配置的初始/单点壁温；逐面元实际值单独保存。
h.wallTemperature=cfg.flight.Twall;
h.wallTemperatureByFace=Tw;
h.wallTemperatureMethod='按逐面元壁温重算参考性质、沿程边界层及局部焓差';
h.Prandtl=Pr; h.transitionRe=criticalRe;
h.laminarHeatFlux=lam; h.turbulentHeatFlux=turb; h.heatFlux=selected;
h.heatTransferCoefficient=selectedCoefficient;
h.recoveryTemperature=selectedRecovery;
h.edgeReynolds=ReEdge; h.referenceReynoldsLaminar=ReLam;
h.referenceReynoldsTurbulent=ReTurb;
h.recoveryTemperatureLaminar=TawLam; h.recoveryTemperatureTurbulent=TawTurb;
h.effectiveCpLaminar=cpHeatLam; h.effectiveCpTurbulent=cpHeatTurb;
h.eligibleFace=calculated; h.laminarFace=laminar; h.turbulentFace=turbulent;
h.unresolvedFace=~calculated; h.nearTipFace=nearTip;
h.transition=transition;
h.streamlineDistanceAssumption=distanceMethod;
h.validity=['附着侧壁工程估算；已验证球形鼻端 <=30° 偏置驻点' ...
    '余弦延拓仅为低阶空间估算；前体脱体闭合属于更低阶外流估算；' ...
    '疑似层流/Head 湍流分离及其下游保持无解，因无再附模型；' ...
    '滑移/过渡流及其准流线下游不套用连续介质公式；' ...
    '湍流模型与回退范围见 turbulentScope，不预测真实三维分离线；' ...
    '激波干扰及真实气体效应未验证；' airScope];
h.summary=struct('calculatedFaces',nnz(calculated), ...
    'unresolvedFaces',nnz(~calculated), ...
    'laminarFaces',nnz(laminar),'turbulentFaces',nnz(turbulent), ...
    'noseFaces',nnz(noseFill), ...
    'detachedClosureHeatFaces',nnz(detachedHeat), ...
    'continuumUnresolvedFaces',nnz(continuum.unresolvedFace), ...
    'separationUnresolvedFaces',nnz(h.separationUnresolvedFace | ...
        h.turbulentSeparationUnresolvedFace), ...
    'turbulentNumericalUnresolvedFaces',nnz(h.turbulentNumericalUnresolvedFace));
if any(calculated)
    h.summary.minimumHeatFlux=min(selected(calculated));
    h.summary.maximumHeatFlux=max(selected(calculated));
else
    h.summary.minimumHeatFlux=NaN; h.summary.maximumHeatFlux=NaN;
end
mesh.heating=h;
if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
mesh=analyzeHeatingCoverage(mesh);
end
