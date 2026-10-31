function mesh = calcStagnationHeating(mesh,freestream,cfg)
%CALCSTAGNATIONHEATING 使用 Sutton-Graves 经验关联式估算钝鼻驻点热流。
% q_s = k * sqrt(rho_inf / R_n) * V_inf^3 [W/m^2]，
% 地球空气默认 k=1.83e-4（上述 SI 单位组合）。
% 鼻尖半径来自用户确认值，或通过多尺度检查的 STL 球形钝鼻拟合。
% 单一网格曲率、三角形尺寸或未经验证的局部球拟合不能代替 R_n。
% 冷壁经验值取配置的参考壁温（默认 300 K）。默认使用温度相关
% 干空气焓差比例修正，并输出驻点换热系数供单点壁温计算；不含辐射、
% 催化和真实气体细节，不能把单点结果填充为整块鼻部热流分布。
% 对多尺度验证通过的球形鼻端，驻点可沿逆来流方向偏移到球面上；
% 最大合成入射角为 30°。这只是旋转对称球的局部单点工程近似，
% 不是大攻角下真实三维激波层或球面热流分布的求解。
% 未经球形验证的用户半径仍只允许 <=5° 小攻角。
% 连续介质 Knudsen 数达到筛选阈值时不套用 Sutton-Graves 公式；
% 本模块不提供滑移流或自由分子流驻点热流。

%% 1. 取得几何候选和经过验证或用户确认的物理半径
if nargin<3, cfg=struct; end
if ~isfield(mesh,'stagnation')
    mesh=findStagnationRegion(mesh,freestream,cfg);
end
s=mesh.stagnation;
result=struct;
result.method='Sutton-Graves 冷壁经验式 + 总温焓差比例修正';
result.formula='q_ref=k*sqrt(rho_inf/R_n)*V_inf^3; q=q_ref*(T0-Tw)/(T0-Tref)';
result.available=false; result.heatFlux=NaN;
result.referenceHeatFlux=NaN;
result.heatTransferCoefficient=NaN;
result.adiabaticWallTemperature=NaN;
result.wallTemperature=NaN;
result.shockStateMethod='未提供；使用自由来流总温';
result.referenceWallTemperature=300;
result.radius=NaN; result.coefficient=1.83e-4;
result.radiusSource='未提供物理鼻尖半径';
result.faceID=s.representativeFaceID;
result.position=s.representativePosition;
result.positionMethod='几何驻点候选代表面面心';
result.incidenceAngleDeg=NaN;
result.freestreamMeanFreePath=NaN;
result.freestreamKnudsen=NaN;
result.maxContinuumKnudsen=0.01;
result.scope=['小攻角已知物理半径，或 <=30° 已验证球形钝鼻的' ...
    '偏置驻点单点估算；焓差修正仅为工程闭合，不是鼻部热流分布'];
[airMode,airScope]=selectAirThermalModel(freestream,cfg);
result.thermalPropertyMode=airMode;
result.thermalPropertyScope=airScope;
if strcmp(airMode,'thermallyPerfect')
    result.method='Sutton-Graves 冷壁经验式 + 干空气 NASA7 焓差比例修正';
    result.formula='q=q_ref*[h(T0)-h(Tw)]/[h(T0)-h(Tref)]';
else
    result.method='Sutton-Graves 冷壁经验式 + 定比热温差比例修正';
end
result.scope=[result.scope '；' airScope];
if isfield(cfg,'heating') && isfield(cfg.heating,'stagnationCoefficient')
    result.coefficient=cfg.heating.stagnationCoefficient;
end
if isfield(cfg,'heating') && ...
        isfield(cfg.heating,'stagnationReferenceWallTemperature')
    result.referenceWallTemperature= ...
        cfg.heating.stagnationReferenceWallTemperature;
end
validateattributes(result.coefficient,{'numeric'}, ...
    {'real','finite','scalar','positive'});
validateattributes(result.referenceWallTemperature,{'numeric'}, ...
    {'real','finite','scalar','positive'});
if isfield(s,'hasPhysicalRadius') && s.hasPhysicalRadius
    result.radius=s.physicalRadius;
    result.radiusSource=s.radiusSource;
end
if ~isfinite(result.radius) || result.radius<=0
    result.status='未计算：STL 未通过球形钝鼻检查，且未提供确认的物理鼻尖半径';
    mesh.stagnationHeating=result; return;
end

%% 2. 使用物理鼻端半径检查连续介质经验式的尺度适用性
if isfield(cfg,'heating') && isfield(cfg.heating,'maxContinuumKnudsen')
    result.maxContinuumKnudsen=cfg.heating.maxContinuumKnudsen;
end
validateattributes(result.maxContinuumKnudsen,{'numeric'}, ...
    {'real','finite','scalar','>',0,'<=',0.01});
result.scope=[result.scope sprintf( ...
    '；仅 Kn_{∞,Rn}<%.3g 的连续介质适用', ...
    result.maxContinuumKnudsen)];
for name={'pressure','temperature','R'}
    if ~isfield(freestream,name{1})
        error('RocketAeroHeat:StagnationHeating:Freestream', ...
            '连续介质筛选缺少 freestream.%s。',name{1});
    end
end
muInf=NaN;
if isfield(freestream,'dynamicViscosity')
    muInf=freestream.dynamicViscosity;
end
if ~isfinite(muInf) && isfinite(freestream.temperature) && ...
        freestream.temperature>0
    T=freestream.temperature;
    muInf=1.458e-6*T^(3/2)/(T+110.4);
end
result.freestreamMeanFreePath=airMeanFreePath( ...
    freestream.pressure,freestream.temperature,muInf,freestream.R);
result.freestreamKnudsen=result.freestreamMeanFreePath/result.radius;
if ~isfinite(result.freestreamKnudsen) || ...
        result.freestreamKnudsen>=result.maxContinuumKnudsen
    result.status=sprintf(['未计算：Kn_{∞,Rn}=%.3g 达到连续介质' ...
        '筛选阈值 %.3g；当前没有稀薄流驻点热流模型'], ...
        result.freestreamKnudsen,result.maxContinuumKnudsen);
    mesh.stagnationHeating=result; return;
end

%% 3. 明确经验式适用域
for name={'density','speed','Mach','alpha','beta','totalTemperature'}
    key=name{1};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:StagnationHeating:Freestream', ...
            '缺少 freestream.%s。',key);
    end
    validateattributes(freestream.(key),{'numeric'}, ...
        {'real','finite','scalar'},mfilename,['freestream.' key]);
end
if strcmp(airMode,'thermallyPerfect')
    stagnationTotalTemperature=airRecoveryTemperature( ...
        freestream.temperature,freestream.speed,1);
    result.shockStateMethod='未提供法向激波；使用自由来流总焓反解总温';
else
    stagnationTotalTemperature=freestream.totalTemperature;
end
if isfield(mesh,'stagnationShock') && mesh.stagnationShock.available
    shock=mesh.stagnationShock;
    shockMode='caloricallyPerfect';
    if isfield(shock,'thermalPropertyMode')
        shockMode=shock.thermalPropertyMode;
    end
    if ~strcmp(shockMode,airMode) || ...
            ~isfinite(stagnationTotalTemperature) || ...
            abs(shock.stagnationTemperature/stagnationTotalTemperature-1)>1e-6
        error('RocketAeroHeat:StagnationHeating:StaleShock', ...
            '驻点法向激波状态与当前空气热物性模式或来流总焓不一致。');
    end
    stagnationTotalTemperature=shock.stagnationTemperature;
    result.shockStateMethod=shock.method;
end
incidence=hypot(freestream.alpha,freestream.beta);
result.incidenceAngleDeg=incidence;
verifiedSphere=isfield(s,'noseRadiusEstimate') && ...
    s.noseRadiusEstimate.available && ...
    abs(result.radius/s.noseRadiusEstimate.radius-1)<=0.05 && ...
    isfield(mesh,'V') && isfield(freestream,'unitFlowDirection') && ...
    s.noseRadiusEstimate.outerDepthRadiusFraction>= ...
        1+freestream.unitFlowDirection(1)-1e-6;
if verifiedSphere
    % 经验证的球鼻：用球心到逆来流方向的半径定位偏置驻点。
    % 候选代表面仅用于检查局部法向和连通性，不能代替连续驻点坐标。
    V=double(mesh.V); R=s.noseRadiusEstimate.radius;
    xTip=max(V(:,1));
    tipIDs=abs(V(:,1)-xTip)<=max(1e-9,1e-8*R);
    tipYZ=mean(V(tipIDs,2:3),1);
    upstream=-double(freestream.unitFlowDirection(:)');
    result.position=[xTip-R,tipYZ]+R*upstream;
    result.positionMethod='已验证球形鼻端沿逆来流方向的解析球面交点';
end
if freestream.density<=0 || freestream.speed<=0 || freestream.Mach<=1
    result.status='未计算：经验式要求超声速正密度来流';
elseif incidence>30
    result.status='未计算：合成入射角超过 30°，球鼻偏置驻点近似不适用';
elseif incidence>5 && ~verifiedSphere
    result.status=['未计算：超过 5° 时必须验证球形鼻端、偏置驻点深度' ...
        '和物理半径一致性'];
elseif s.numberComponents~=1
    result.status='未计算：存在多个驻点候选区';
elseif isfield(s,'flags') && s.flags.flatCapCandidate
    result.status='未计算：驻点候选涉及平封盖';
elseif isfield(s,'peakMu') && s.peakMu<cosd(15)
    result.status='未计算：候选鼻端法向未与来流充分对齐';
elseif ~isfinite(stagnationTotalTemperature)
    result.status='未计算：来流总焓对应温度超出干空气热物性适用域';
elseif stagnationTotalTemperature<=result.referenceWallTemperature
    result.status='未计算：来流总温不高于冷壁参考温度';
elseif strcmp(airMode,'thermallyPerfect') && ...
        (stagnationTotalTemperature<80 || stagnationTotalTemperature>2500)
    result.status='未计算：驻点总温超出干空气焓模型 80~2500 K 适用域';
else
    %% 4. 冷壁参考值及当前壁温下的单点对流边界条件
    if ~isfield(cfg,'flight') || ~isfield(cfg.flight,'Twall')
        error('RocketAeroHeat:StagnationHeating:WallTemperature', ...
            '驻点焓差修正需要 cfg.flight.Twall。');
    end
    Tw=cfg.flight.Twall;
    validateattributes(Tw,{'numeric'}, ...
        {'real','finite','scalar','positive'});
    result.referenceHeatFlux=result.coefficient* ...
        sqrt(freestream.density/result.radius)*freestream.speed^3;
    if strcmp(airMode,'thermallyPerfect')
        cpReference=airEnthalpySecantCp( ...
            stagnationTotalTemperature,result.referenceWallTemperature);
        cpWall=airEnthalpySecantCp(stagnationTotalTemperature,Tw);
        if ~isfinite(cpReference) || ~isfinite(cpWall)
            result.status='未计算：壁温或参考壁温超出干空气焓模型适用域';
            mesh.stagnationHeating=result; return;
        end
        denominator=cpReference*(stagnationTotalTemperature- ...
            result.referenceWallTemperature);
        result.heatTransferCoefficient= ...
            result.referenceHeatFlux*cpWall/denominator;
    else
        result.heatTransferCoefficient=result.referenceHeatFlux/ ...
            (stagnationTotalTemperature-result.referenceWallTemperature);
    end
    result.adiabaticWallTemperature=stagnationTotalTemperature;
    result.wallTemperature=Tw;
    result.heatFlux=result.heatTransferCoefficient* ...
        (stagnationTotalTemperature-Tw);
    result.available=true;
    result.status='钝鼻驻点单点经验估算；已按壁温作一阶焓差修正';
end
mesh.stagnationHeating=result;
end
