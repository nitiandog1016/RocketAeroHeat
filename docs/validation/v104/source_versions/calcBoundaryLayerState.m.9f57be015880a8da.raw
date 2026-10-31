function mesh = calcBoundaryLayerState(mesh,freestream,cfg)
%CALCBOUNDARYLAYERSTATE 沿内置三维准流线积分层流动量厚度。
% 模块作用：把外流速度和表面路径转成 theta、Re_theta、局部压力梯度参数
% 及层流摩擦系数，供热流模块使用。采用变黏度 Thwaites 型工程积分：
% theta^2=0.45/U_e^6 * integral(nu_ref*U_e^5 ds)，
% lambda=theta^2/nu_ref*dU_e/ds，Cf=2*l(lambda)/Re_theta。
% l(lambda)=0.2227+1.57*lambda-1.8*lambda^2。
% lambda<=-0.09 只标记疑似层流分离；超出相关式范围时不外推。
% 速度梯度默认在上游 4 个准流线步长上求割线，抑制面元离散及
% 外流角度节点造成的局部数值跳变；仍不能替代连续流场梯度。
% 默认干空气以 h(Taw)=h(Te)+sqrt(Pr)*Ue²/2 求恢复温度，
% 与侧壁热流使用同一焓模型；显式定比热模式保留旧公式。
% 此法是附着层流的准二维近似，不是三维边界层或湍流分离求解器。

%% 1. 检查外流与路径，计算 Eckert 参考温度下的运动黏度
if nargin<3, cfg=struct; end
if ~isfield(mesh,'boundaryLayerEdge') || ~isfield(mesh,'surfaceStreamlines')
    error('RocketAeroHeat:BoundaryLayer:Input', ...
        '需要先计算内置外流并追踪三维表面准流线。');
end
if ~isfield(cfg,'flight') || ~isfield(cfg.flight,'Twall')
    error('RocketAeroHeat:BoundaryLayer:WallTemperature','需要壁温 Twall [K]。');
end
Tw=surfaceWallTemperature(size(mesh.F,1),cfg);
Pr=0.72;
if isfield(cfg,'heating') && isfield(cfg.heating,'Prandtl')
    Pr=cfg.heating.Prandtl;
end
validateattributes(Pr,{'numeric'},{'real','finite','scalar','>',0,'<',1});
gradientSteps=4;
if isfield(cfg,'heating') && isfield(cfg.heating,'gradientPathSteps')
    gradientSteps=cfg.heating.gradientPathSteps;
end
validateattributes(gradientSteps,{'numeric'}, ...
    {'real','finite','scalar','integer','>=',1,'<=',20});
e=mesh.boundaryLayerEdge; path=mesh.surfaceStreamlines;
nf=size(mesh.F,1);
if ~isfield(path,'predecessorFace') || ~isfield(path,'stepLength') || ...
        ~isfield(path,'pathCost') || ~strcmp(path.velocitySource,e.method)
    error('RocketAeroHeat:BoundaryLayer:Path', ...
        '流线路径缺少前驱/几何步长，或与当前外流状态不一致。');
end
U=e.speed(:); Te=e.temperature(:); pe=e.pressure(:);
if numel(U)~=nf || numel(path.distance)~=nf
    error('RocketAeroHeat:BoundaryLayer:Size','外流与网格面元数不一致。');
end
[airMode,airScope]=selectAirThermalModel(freestream,cfg);
if strcmp(airMode,'thermallyPerfect')
    Taw=airRecoveryTemperature(Te,U,sqrt(Pr));
    recoveryMethod='固定组分干空气 NASA7/低温振动冻结焓恢复';
else
    cp=freestream.gamma*freestream.R/(freestream.gamma-1);
    Taw=Te+sqrt(Pr)*U.^2/(2*cp);
    recoveryMethod='定比热温差恢复';
end
Tref=Te+0.5*(Tw-Te)+0.22*(Taw-Te);
rhoRef=pe./(freestream.R*Tref);
muRef=1.458e-6*Tref.^(3/2)./(Tref+110.4);
nuRef=muRef./rhoRef;

%% 2. 按 Dijkstra 选出的路径顺序，沿真实几何步长积分
theta=nan(nf,1); ReTheta=theta; lambda=theta; Cf=theta;
gradient=theta; integralTerm=theta;
eligible=path.connectedFace(:) & e.algebraicallyValid(:) & ...
    isfinite(U) & U>0 & isfinite(nuRef) & nuRef>0;
seed=eligible & path.predecessorFace(:)==0;
integralTerm(seed)=0;
[~,order]=sort(path.pathCost(eligible));
ids=find(eligible); ids=ids(order);
for k=1:numel(ids)
    id=ids(k); parent=path.predecessorFace(id);
    if parent==0, continue; end
    ds=path.stepLength(id);
    if ~eligible(parent) || ~isfinite(integralTerm(parent)) || ...
            ~isfinite(ds) || ds<=0
        continue;
    end
    integralTerm(id)=integralTerm(parent)+0.5*ds*( ...
        nuRef(parent)*U(parent)^5+nuRef(id)*U(id)^5);
    theta(id)=sqrt(0.45*integralTerm(id)/U(id)^6);
    ReTheta(id)=U(id)*theta(id)/nuRef(id);
    % 沿同一条前驱链取多步割线，避免局部数值跳变被单个很短的
    % 网格步长放大。靠近种子点时使用已有的较短路径。
    upstream=parent; span=ds;
    for step=2:gradientSteps
        older=path.predecessorFace(upstream);
        if older==0 || ~eligible(older) || ...
                ~isfinite(path.stepLength(upstream)) || ...
                path.stepLength(upstream)<=0
            break;
        end
        span=span+path.stepLength(upstream);
        upstream=older;
    end
    gradient(id)=(U(id)-U(upstream))/span;
    lambda(id)=theta(id)^2/nuRef(id)*gradient(id);
    if isfinite(ReTheta(id)) && ReTheta(id)>0 && ...
            isfinite(lambda(id)) && lambda(id)>-0.09 && lambda(id)<=0.1
        l=0.2227+1.57*lambda(id)-1.8*lambda(id)^2;
        Cf(id)=2*l/ReTheta(id);
    end
end

%% 3. 区分可用解、疑似层流分离及相关式范围外面元
separation=eligible & isfinite(lambda) & lambda<=-0.09;
outside=eligible & isfinite(lambda) & lambda>0.1;
valid=isfinite(Cf) & Cf>0;
b=struct;
b.method=['沿三维准流线的变黏度 Thwaites 型层流动量积分 + ' ...
    recoveryMethod];
b.scope=sprintf(['附着层流工程近似；速度梯度取上游最多 %d 步' ...
    '准流线割线；%s；不能解析三维边界层或分离再附'], ...
    gradientSteps,airScope);
b.gradientPathSteps=gradientSteps;
b.thermalPropertyMode=airMode;
b.recoveryTemperatureMethod=recoveryMethod;
b.recoveryTemperature=Taw;
b.momentumThickness=theta;
b.momentumThicknessReynolds=ReTheta;
b.pressureGradientParameter=lambda;
b.edgeVelocityGradient=gradient;
b.skinFrictionCoefficient=Cf;
b.referenceTemperature=Tref;
b.referenceKinematicViscosity=nuRef;
b.wallTemperatureByFace=Tw;
b.thermalPropertyUnavailableFace=e.algebraicallyValid(:) & ...
    (~isfinite(Taw) | ~isfinite(Tref) | ~isfinite(nuRef));
b.pathIntegral=integralTerm;
b.validFace=valid;
b.possibleLaminarSeparationFace=separation;
b.correlationOutOfRangeFace=outside;
b.unresolvedFace=~valid;
b.summary=struct('validFaces',nnz(valid), ...
    'possibleLaminarSeparationFaces',nnz(separation), ...
    'correlationOutOfRangeFaces',nnz(outside), ...
    'thermalPropertyUnavailableFaces',nnz(b.thermalPropertyUnavailableFace));
mesh.boundaryLayer=b;
if isfield(mesh,'heating'), mesh=rmfield(mesh,'heating'); end
if isfield(mesh,'heatingCoverage'), mesh=rmfield(mesh,'heatingCoverage'); end
if isfield(mesh,'wallThermal'), mesh=rmfield(mesh,'wallThermal'); end
if isfield(mesh,'aerothermalExposure'), mesh=rmfield(mesh,'aerothermalExposure'); end
if isfield(mesh,'flightHeatingHistory'), mesh=rmfield(mesh,'flightHeatingHistory'); end
end
