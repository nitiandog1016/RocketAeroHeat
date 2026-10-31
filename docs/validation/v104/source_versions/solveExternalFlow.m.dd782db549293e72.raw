function mesh = solveExternalFlow(mesh,freestream,cfg)
%SOLVEEXTERNALFLOW 从当前 STL 面元与飞行条件直接计算工程外流状态。
% 零攻角迎风曲面：局部切锥 Taylor-Maccoll 附体激波近似；
% 非零攻角迎风面：局部切楔二维附体斜激波近似；
% 背风凸转角：局部 Prandtl-Meyer 等熵膨胀近似；近切向面取自由来流。
% 面元彼此独立，因此无法解析脱体激波、肩部波系相互作用、分离与尾迹。
% 无附体解、驻点候选、封盖和特征折角返回 NaN，并标记对应面元。
% 非零攻角光滑前体若局部切楔超过附体极限，可选用修正牛顿壁压的
% 角度变化趋势，在附体极限状态与正激波滞止状态之间作受限闭合。
% 不调用外部 CFD；结果是有明确适用边界的内置工程外流模型。

%% 1. 验证几何与工况，取得当前入射角
if nargin<3, cfg=struct; end
if ~isfield(mesh,'localGeometry')
    error('RocketAeroHeat:ExternalSolver:Geometry','需要最终网格 localGeometry。');
end
mesh=classifyAerothermalRegion(mesh,freestream,cfg);
g=mesh.localGeometry; region=mesh.aerothermalRegion;
nf=size(mesh.F,1); N=g.normal;
flow=region.unitFlowDirection;
normalFrame=reconstructSectionNormals(mesh,cfg);
mu=-normalFrame.normal*flow';
angle=asind(max(-1,min(1,mu)));
angleStep=0.5;
if isfield(cfg,'externalFlow') && isfield(cfg.externalFlow,'angleStepDeg')
    angleStep=cfg.externalFlow.angleStepDeg;
end
validateattributes(angleStep,{'numeric'}, ...
    {'real','finite','scalar','>',0,'<=',2});
detachedClosureEnabled=true;
if isfield(cfg,'externalFlow') && ...
        isfield(cfg.externalFlow,'detachedClosureEnabled')
    detachedClosureEnabled=cfg.externalFlow.detachedClosureEnabled;
end
validateattributes(detachedClosureEnabled,{'logical','numeric'}, ...
    {'binary','scalar'});
if freestream.Mach<=1
    error('RocketAeroHeat:ExternalSolver:Mach','内置激波/膨胀模型要求 M>1。');
end
if ~isfield(mesh,'stagnation')
    mesh=findStagnationRegion(mesh,freestream,cfg);
end

%% 2. 球形钝鼻驻点的局部法向激波；不把它填成逐面元外缘场
stag=mesh.stagnation;
[airMode,~]=selectAirThermalModel(freestream,cfg);
% 侧壁切楔优先采用与鼻端相同的 NASA7 干空气热力模型；当总温
% 超出适用域时整组切楔退回有明确标注的定比热工程解。
thermalWedge=false;
wedgeThermalStatus='配置使用定比热切楔';
thermalCone=false;
coneThermalStatus='配置使用定比热切锥';
thermalExpansion=false;
expansionThermalStatus='零攻角轴向推进或配置使用定比热膨胀';
if strcmp(airMode,'thermallyPerfect')
    thermalConeReference=solveThermallyPerfectConicalShock(freestream,0);
    thermalCone=thermalConeReference.available;
    if thermalCone
        coneThermalStatus='切锥采用 NASA7 激波守恒与变比热 Taylor-Maccoll 积分';
    else
        coneThermalStatus=['切锥退回定比热：' thermalConeReference.status];
    end
    thermalWedgeReference=solveThermallyPerfectObliqueShock(freestream,0);
    thermalWedge=thermalWedgeReference.available;
    if thermalWedge
        wedgeThermalStatus='切楔采用 NASA7 焓/熵守恒附体弱斜激波';
    else
        wedgeThermalStatus=['切楔退回定比热：' thermalWedgeReference.status];
    end
    if hypot(freestream.alpha,freestream.beta)>1e-8
        thermalExpansionReference= ...
            solveThermallyPerfectPrandtlMeyerExpansion(freestream,0);
        thermalExpansion=thermalExpansionReference.available;
        if thermalExpansion
            expansionThermalStatus='背风局部膨胀采用 NASA7/低温振动冻结焓熵简单波';
        else
            expansionThermalStatus=['背风局部膨胀退回定比热：' ...
                thermalExpansionReference.status];
        end
    end
end
thermalGrazing=(hypot(freestream.alpha,freestream.beta)<1e-8 && ...
    thermalCone) || (hypot(freestream.alpha,freestream.beta)>1e-8 && ...
    thermalWedge && thermalExpansion);
if strcmp(airMode,'thermallyPerfect')
    stagnationMethod='固定组分干空气 NASA7 焓/熵法向激波';
else
    stagnationMethod='完全气体法向激波 Rankine-Hugoniot 跃迁';
end
stagnationShock=struct('available',false, ...
    'method',stagnationMethod, ...
    'thermalPropertyMode',airMode, ...
    'scope',['已验证球形鼻端、合成入射角 <=30° 的偏置驻点局部法向激波；' ...
        '仅单点热力状态，不是非对称三维激波面'], ...
    'status','鼻端几何或攻角不满足驻点法向激波近似');
incidence=hypot(freestream.alpha,freestream.beta);
verifiedSphere=isfield(stag,'noseRadiusEstimate') && ...
    stag.noseRadiusEstimate.available && ...
    abs(stag.physicalRadius/stag.noseRadiusEstimate.radius-1)<=0.05 && ...
    stag.noseRadiusEstimate.outerDepthRadiusFraction>= ...
        1+freestream.unitFlowDirection(1)-1e-6;
if stag.hasPhysicalRadius && stag.numberComponents==1 && ...
        ~stag.flags.flatCapCandidate && stag.peakMu>=cosd(15) && ...
        (incidence<=5 || (incidence<=30 && verifiedSphere))
    if strcmp(airMode,'thermallyPerfect')
        stagnationShock=solveThermallyPerfectNormalShock(freestream);
    else
        stagnationShock=solveNormalShock(freestream);
        stagnationShock.thermalPropertyMode='caloricallyPerfect';
    end
    stagnationShock.radius=stag.physicalRadius;
    stagnationShock.radiusSource=stag.radiusSource;
    stagnationShock.candidatePosition=stag.representativePosition;
    stagnationShock.incidenceAngleDeg=incidence;
    stagnationShock.scope=[stagnationShock.scope ...
        '；已验证球形鼻端可至 30° 偏置驻点；不代表非对称三维激波几何'];
end
mesh.stagnationShock=stagnationShock;
mesh.bowShockPreview=estimateBowShockSurface(mesh,freestream);
if mesh.bowShockPreview.available
    mesh.stagnationShock.estimatedStandOffDistance= ...
        mesh.bowShockPreview.standOffDistance;
    mesh.stagnationShock.standOffMethod=mesh.bowShockPreview.method;
else
    mesh.stagnationShock.estimatedStandOffDistance=NaN;
    mesh.stagnationShock.standOffMethod=mesh.bowShockPreview.status;
end

%% 3. 逐面元分区：在模型不适用的地方不制造外缘状态
baseValid=g.isSidewall & ~g.featureFace & ~g.topologyIssue & ...
    ~mesh.stagnation.candidateFace;
grazing=baseValid & abs(angle)<=region.grazingAngleDeg;
compress=baseValid & angle>region.grazingAngleDeg;
expand=baseValid & angle<-region.grazingAngleDeg;
zeroAngle=hypot(freestream.alpha,freestream.beta)<1e-8;
cone=compress & zeroAngle;
wedge=compress & ~zeroAngle;

%% 4. 初始化结果；0 表示无解，不用自由来流伪装解
p=nan(nf,1); T=p; rho=p; U=p; M=p; p0=p; T0=p;
shockAngle=p;
methodID=zeros(nf,1,'uint8');
% 近切向状态的表面压力/热力参数与自由来流相同。
p(grazing)=freestream.pressure;
T(grazing)=freestream.temperature;
rho(grazing)=freestream.pressure/(freestream.R*freestream.temperature);
U(grazing)=freestream.speed; M(grazing)=freestream.Mach;
p0(grazing)=freestream.totalPressure;
T0(grazing)=freestream.totalTemperature;
methodID(grazing)=1;
if thermalGrazing
    % 非零攻角的切楔/膨胀已使用同一总焓与标准熵；近切向自由流
    % 必须给出相同的 NASA7 总状态，避免沿圆柱周向出现人为总压跳变。
    if hypot(freestream.alpha,freestream.beta)<1e-8
        thermalReference=thermalConeReference;
    else
        thermalReference=thermalWedgeReference;
    end
    p0(grazing)=thermalReference.upstreamTotalPressureThermal;
    T0(grazing)=thermalReference.totalTemperature;
    M(grazing)=thermalReference.upstreamMachThermal;
end

%% 5. 按角度节点缓存解析解，在相邻节点间连续插值
% 逐面元直接求根过慢；最近邻分箱又会制造台阶速度，使离散
% Thwaites 梯度误判分离。这里仅对角度节点求解，再在节点间插值
% 静压与总压，用完全气体关系恢复其余热力状态。
if any(cone)
    indices=find(cone);
    if thermalCone
        coneSolver=@(turn)solveThermallyPerfectConicalShock(freestream,turn);
    else
        coneSolver=@(turn)solveConicalShock(freestream,turn);
    end
    interpolated=interpolateLocalState(angle(cone),angleStep, ...
        freestream,coneSolver,'attached',inf,true,thermalCone);
    [p,T,rho,U,M,p0,T0,shockAngle,methodID]=assignInterpolated( ...
        indices,interpolated,2,p,T,rho,U,M,p0,T0,shockAngle,methodID);
end
detachedWedge=false(nf,1);
if any(wedge)
    if thermalWedge
        wedgeSolver=@(turn)solveThermallyPerfectObliqueShock(freestream,turn);
        limit=thermalWedgeReference;
        wedgeValidField='available';
    else
        wedgeSolver=@(turn)solveObliqueShock(freestream,turn);
        limit=solveObliqueShock(freestream,0);
        wedgeValidField='attached';
    end
    maximumAttached=limit.maximumAttachedDeflectionDeg;
    detachedWedge=wedge & angle>maximumAttached;
    attachedWedge=wedge & ~detachedWedge;
    if any(attachedWedge)
        indices=find(attachedWedge);
        interpolated=interpolateLocalState(angle(attachedWedge), ...
            angleStep,freestream,wedgeSolver, ...
            wedgeValidField,maximumAttached,true,thermalWedge);
        [p,T,rho,U,M,p0,T0,shockAngle,methodID]=assignInterpolated( ...
            indices,interpolated,3,p,T,rho,U,M,p0,T0,shockAngle,methodID);
    end
end
thermalExpansionFace=false(nf,1);
expansionFallbackFace=false(nf,1);
thermalExpansionOutOfRangeFace=false(nf,1);
if any(expand)
    if thermalExpansion
        maximumExpansion=thermalExpansionReference.maximumAdditionalTurnDeg;
        inRange=expand & -angle<=maximumExpansion;
        thermalExpansionOutOfRangeFace=expand & ~inRange;
        expansionFallbackFace=thermalExpansionOutOfRangeFace;
        if any(inRange)
            indices=find(inRange);
            interpolated=interpolateLocalState(-angle(inRange),angleStep, ...
                freestream,@(turn) ...
                solveThermallyPerfectPrandtlMeyerExpansion(freestream,turn), ...
                'available',maximumExpansion,false,true);
            [p,T,rho,U,M,p0,T0,shockAngle,methodID]=assignInterpolated( ...
                indices,interpolated,4,p,T,rho,U,M,p0,T0,shockAngle,methodID);
            thermalExpansionFace(indices)=methodID(indices)==4;
        end
    else
        expansionFallbackFace=expand;
    end
    if any(expansionFallbackFace)
        indices=find(expansionFallbackFace);
        interpolated=interpolateLocalState( ...
            -angle(expansionFallbackFace),angleStep,freestream, ...
            @(turn)solvePrandtlMeyerExpansion(freestream,turn), ...
            'available',inf,false);
        [p,T,rho,U,M,p0,T0,shockAngle,methodID]=assignInterpolated( ...
            indices,interpolated,4,p,T,rho,U,M,p0,T0,shockAngle,methodID);
    end
    if thermalExpansion
        expansionThermalStatus=sprintf([ ...
            '背风局部膨胀采用 NASA7/低温振动冻结焓熵简单波；' ...
            '%d 面元低于 80 K 下限，定比热状态只作流向占位，热流留空'], ...
            nnz(expansionFallbackFace));
    end
end

%% 6. 局部切楔脱体时的受限工程闭合；不推断激波角或激波面
% 仅处理光滑前体 10% 轴向长度、局部转角 <=60° 的面元。
% 附体极限处锚定切楔弱解的静压/总压；将修正牛顿压力从附体极限
% 至正对来流的变化归一化，作为朝正激波滞止状态的插值权重。
% 同一个权重插值静压和总压，保证两者连续且静压始终小于总压。
% 定比热等熵关系只用于恢复局部外缘状态；这不是脱体波层求解。
closureStatus='关闭或没有超出附体极限的局部切楔面元';
if detachedClosureEnabled && any(detachedWedge)
    closureStatus='缺少与当前工况一致的修正牛顿壁压';
    pressureCurrent=isfield(mesh,'pressure') && ...
        isfield(mesh.pressure,'freestream') && ...
        isfield(mesh.pressure,'staticPressure') && ...
        numel(mesh.pressure.staticPressure)==nf;
    if pressureCurrent
        saved=mesh.pressure.freestream;
        for key={'Mach','alpha','beta','altitude','gamma','R', ...
                'pressure','temperature','density','speed', ...
                'totalPressure','totalTemperature'}
            name=key{1};
            pressureCurrent=pressureCurrent && isfield(saved,name) && ...
                isequal(saved.(name),freestream.(name));
        end
    end
    if pressureCurrent
        referenceNormal=solveNormalShock(freestream);
        if thermalWedge
            normalShock=solveThermallyPerfectNormalShock(freestream);
            angleLimit=thermalWedgeReference;
            attachedLimit=wedgeSolver(angleLimit.maximumAttachedDeflectionDeg);
        else
            normalShock=referenceNormal;
            angleLimit=solveObliqueShock(freestream,0);
            attachedLimit=solveObliqueShock( ...
                freestream,angleLimit.maximumAttachedDeflectionDeg);
        end
        maximumAttached=angleLimit.maximumAttachedDeflectionDeg;
        wallPressure=double(mesh.pressure.staticPressure(:));
        xTip=max(mesh.V(:,1));
        bodyLength=xTip-min(mesh.V(:,1));
        fallback=detachedWedge & normalShock.available & ...
            attachedLimit.attached & isfinite(attachedLimit.pressure) & ...
            (xTip-g.center(:,1))<=0.1*bodyLength & ...
            angle<=60 & ~region.isNoseTipAdjacent & ...
            ~region.isAftTipAdjacent & ...
            isfinite(wallPressure) & wallPressure>=freestream.pressure & ...
            wallPressure<referenceNormal.stagnationPressure;
        if any(fallback)
            wallPressureAtLimit=freestream.pressure+ ...
                (referenceNormal.stagnationPressure-freestream.pressure)* ...
                sind(maximumAttached)^2;
            fraction=max(0,min(1,(wallPressure(fallback)- ...
                wallPressureAtLimit)/ ...
                (referenceNormal.stagnationPressure-wallPressureAtLimit)));
            localPressure=(1-fraction)*attachedLimit.pressure+ ...
                fraction*normalShock.stagnationPressure;
            localTotalPressure=(1-fraction)*attachedLimit.totalPressure+ ...
                fraction*normalShock.stagnationPressure;
            ids=find(fallback);
            if thermalWedge
                recovered=recoverThermallyPerfectEdgeState( ...
                    localPressure,localTotalPressure, ...
                    attachedLimit.totalTemperature,freestream.R);
                keep=recovered.valid;
                ids=ids(keep);
                p(ids)=localPressure(keep);
                T(ids)=recovered.temperature(keep);
                rho(ids)=recovered.density(keep);
                M(ids)=recovered.Mach(keep);
                U(ids)=recovered.speed(keep);
                p0(ids)=localTotalPressure(keep);
                T0(ids)=attachedLimit.totalTemperature;
            else
                localMach=sqrt(2/(freestream.gamma-1)* ...
                    expm1((freestream.gamma-1)/freestream.gamma* ...
                    log(localTotalPressure./localPressure)));
                localTemperature=freestream.totalTemperature./ ...
                    (1+0.5*(freestream.gamma-1)*localMach.^2);
                p(ids)=localPressure;
                T(ids)=localTemperature;
                rho(ids)=localPressure./(freestream.R*localTemperature);
                M(ids)=localMach;
                U(ids)=localMach.*sqrt( ...
                    freestream.gamma*freestream.R*localTemperature);
                p0(ids)=localTotalPressure;
                T0(ids)=freestream.totalTemperature;
            end
            methodID(ids)=5;
        end
        closureStatus=sprintf('附体极限锚定 + 修正牛顿角度趋势/法向激波过渡：%d 面元', ...
            nnz(methodID==5));
    end
end

%% 7. 流向沿来流在局部切平面的投影；若投影为零则无外缘流向解
% mu=-dot(N,flow)，故投影 flow-(flow·N)N = flow+mu*N。
tangent=repmat(flow,nf,1)+region.mu.*N;
tangentLength=vecnorm(tangent,2,2);
valid=methodID>0 & tangentLength>1e-10;
methodID(~valid)=0;
p(~valid)=NaN; T(~valid)=NaN; rho(~valid)=NaN;
U(~valid)=NaN; M(~valid)=NaN; p0(~valid)=NaN; T0(~valid)=NaN;
shockAngle(~valid)=NaN;
velocity=nan(nf,3);
velocity(valid,:)=U(valid).*tangent(valid,:)./tangentLength(valid);
sound=sqrt(freestream.gamma*freestream.R*T);
thermalFace=false(nf,1);
thermalFace=(thermalCone & methodID==2) | ...
    (thermalWedge & (methodID==3 | methodID==5)) | ...
    thermalExpansionFace | (thermalGrazing & methodID==1);
if any(thermalFace)
    [localCp,~]=airThermalProperties(T(thermalFace));
    localGamma=localCp./(localCp-freestream.R);
    sound(thermalFace)=sqrt(localGamma.*freestream.R.*T(thermalFace));
end
viscosity=1.458e-6*T.^(3/2)./(T+110.4);

%% 8. 输出逐面元算法来源和适用性，并供气动热模块自动读取
edge=struct;
edge.method=['内置工程外流：角度节点连续插值的局部切锥/切楔附体激波 + ' ...
    'Prandtl-Meyer 膨胀 + 受限前体脱体闭合'];
edge.method=[edge.method '；' coneThermalStatus '；' ...
    wedgeThermalStatus '；' expansionThermalStatus '；入射角法向：' normalFrame.mode];
edge.incidenceNormal=normalFrame.normal;
edge.normalReconstruction=normalFrame;
edge.incidenceAngleDeg=angle;
edge.scope=['面元独立局部近似；脱体闭合仅光滑前体 10% 长度内' ...
    '的切楔超限压缩面元、转角 <=60°，采用附体极限到法向激波的' ...
    '归一化修正牛顿角度趋势插值，不求波面与波层流动；' ...
    '曲面、肩部、大攻角波系相互作用未解析'];
edge.pressureSource=edge.method;
edge.assumption=['切锥采用所选干空气热物性；' ...
    '切楔/背风膨胀/前体闭合的热物性模式逐面元标记；附体无粘流；' ...
    '脱体闭合采用附体极限状态、修正牛顿壁压角度趋势与法向激波' ...
    '滞止状态插值；局部流向沿面元切平面'];
edge.caseID=sprintf('built-in-M%.5g-a%.5g-b%.5g', ...
    freestream.Mach,freestream.alpha,freestream.beta);
edge.faceMethodID=methodID;
edge.thermalPropertyModeByFace=zeros(nf,1,'uint8');
edge.thermalPropertyModeByFace(methodID>0)=1;
edge.thermalPropertyModeByFace(thermalFace)=2;
edge.thermalPropertyModeByFace(thermalFace & T>=80 & T<150)=3;
edge.thermalPropertyModeNames={'无解','定比热', ...
    'NASA7 干空气变比热','80~150 K 振动冻结干空气近似'};
edge.wedgeThermalStatus=wedgeThermalStatus;
edge.coneThermalStatus=coneThermalStatus;
edge.expansionThermalStatus=expansionThermalStatus;
edge.thermalExpansionOutOfRangeFace=thermalExpansionOutOfRangeFace;
% 模块八补充：域外定比热状态只维持准流线诊断，不是已求物理外缘。
% 保留代数可追踪状态，另给字段有效性；显示和流态判断不能混用两者。
edge.diagnosticPlaceholderFace=valid & thermalExpansionOutOfRangeFace;
edge.thermalStateValidFace=valid & ~edge.diagnosticPlaceholderFace;
edge.stateValidityMethod=['代数可追踪与物理场有效性分开；低于所选热物性下限的' ...
    '定比热膨胀仅作流向诊断占位，不作为外缘p/T或真实稀薄流证据'];
edge.shockAngleDeg=shockAngle;
edge.shockAngleMethod=['切锥 ' coneThermalStatus '；切楔 ' ...
    wedgeThermalStatus ' 的局部附体弱激波角'];
edge.shockAngleScope='逐面元局部角度；不是三维激波面的位置或脱体激波距离';
edge.faceMethodNames={'无解/不适用','近切向自由来流', ...
    '局部切锥 Taylor-Maccoll','局部切楔斜激波', ...
    '局部 Prandtl-Meyer 膨胀','前体脱体工程闭合'};
if thermalCone
    edge.faceMethodNames{3}='局部切锥 NASA7 变比热 Taylor-Maccoll';
end
if thermalWedge
    edge.faceMethodNames{4}='局部切楔 NASA7 焓/熵附体斜激波';
    edge.faceMethodNames{6}='前体 NASA7 附体极限/法向激波工程闭合';
end
if thermalExpansion
    edge.faceMethodNames{5}='局部 NASA7/低温振动冻结焓熵膨胀';
end
if thermalGrazing
    edge.faceMethodNames{2}='近切向 NASA7 自由来流';
end
edge.detachedClosureStatus=closureStatus;
edge.pressure=p; edge.temperature=T; edge.density=rho;
edge.velocityVector=velocity; edge.speed=U; edge.Mach=M;
edge.speedOfSound=sound; edge.dynamicViscosity=viscosity;
edge.totalPressure=p0; edge.totalTemperature=T0;
edge.algebraicallyValid=valid;
edge.numericallyTraceableFace=valid;
edge.physicallyUnresolved=~edge.thermalStateValidFace;
edge.freestream=freestream;
edge.status='内置外流工程估算；无解面元为 NaN';
edge.summary=struct('validFaces',nnz(valid),'unresolvedFaces',nnz(~valid), ...
    'grazingFaces',nnz(methodID==1),'coneFaces',nnz(methodID==2), ...
    'wedgeFaces',nnz(methodID==3),'expansionFaces',nnz(methodID==4), ...
    'detachedClosureFaces',nnz(methodID==5), ...
    'thermalConeFaces',nnz(methodID==2 & thermalFace), ...
    'thermalWedgeFaces',nnz(methodID==3 & thermalFace), ...
    'thermalClosureFaces',nnz(methodID==5 & thermalFace), ...
    'thermalExpansionFaces',nnz(methodID==4 & thermalFace), ...
    'thermalGrazingFaces',nnz(methodID==1 & thermalFace), ...
    'lowTemperatureThermalFaces',nnz(edge.thermalPropertyModeByFace==3), ...
    'thermalExpansionFallbackFaces',nnz(expansionFallbackFace & methodID==4));
mesh.externalFlow=struct('edge',edge,'faceCenter',g.center, ...
    'faceNormal',g.normal,'V',mesh.V,'F',mesh.F, ...
    'freestream',freestream,'method',edge.method,'caseID',edge.caseID, ...
    'stagnationShock',mesh.stagnationShock, ...
    'bowShockPreview',mesh.bowShockPreview);
mesh.boundaryLayerEdge=edge;
if isfield(mesh,'axisymmetricMarch'), mesh=rmfield(mesh,'axisymmetricMarch'); end
for name={'surfaceStreamlines','boundaryLayer','heating', ...
        'heatingCoverage','wallThermal', ...
        'aerothermalExposure','flightHeatingHistory'}
    if isfield(mesh,name{1}), mesh=rmfield(mesh,name{1}); end
end
end

function [p,T,rho,U,M,p0,T0,shockAngle,methodID]=assignInterpolated( ...
    indices,s,code,p,T,rho,U,M,p0,T0,shockAngle,methodID)
% 只写入两端角度节点都有解析解的面元。
local=find(s.valid);
ids=indices(local);
p(ids)=s.pressure(local); T(ids)=s.temperature(local);
rho(ids)=s.density(local); U(ids)=s.speed(local); M(ids)=s.Mach(local);
p0(ids)=s.totalPressure(local); T0(ids)=s.totalTemperature(local);
shockAngle(ids)=s.shockAngleDeg(local);
methodID(ids)=code;
end

function result=interpolateLocalState(angles,step,freestream,solver, ...
    validField,maximumAngle,hasShockAngle,thermalMode)
% 只在离散角度节点调用切锥/切楔/膨胀求解器，不对面元重复求解。
if nargin<8, thermalMode=false; end
a=double(angles(:)); n=numel(a);
blank=nan(n,1);
result=struct('valid',false(n,1),'pressure',blank, ...
    'temperature',blank,'density',blank,'speed',blank,'Mach',blank, ...
    'totalPressure',blank,'totalTemperature',blank, ...
    'shockAngleDeg',blank);
if n==0, return; end
lower=floor(a/step)*step;
upper=lower+step;
% 在附体极限以极限弱解收尾，不能跨入没有附体解的角度节点。
upper=min(upper,maximumAngle);
upper=min(upper,89.999999);
exact=abs(a-lower)<=1e-10;
upper(exact)=lower(exact);
grid=unique([lower;upper]);
states=cell(numel(grid),1);
for k=1:numel(grid)
    states{k}=solver(grid(k));
end
[~,loIndex]=ismember(lower,grid);
[~,hiIndex]=ismember(upper,grid);
[pairs,~,group]=unique([loIndex hiIndex],'rows');
g=freestream.gamma; R=freestream.R;
for k=1:size(pairs,1)
    first=states{pairs(k,1)};
    last=states{pairs(k,2)};
    if ~first.(validField) || ~last.(validField), continue; end
    ids=find(group==k);
    span=upper(ids)-lower(ids);
    weight=zeros(numel(ids),1);
    moving=span>0;
    weight(moving)=(a(ids(moving))-lower(ids(moving)))./span(moving);
    weight=max(0,min(1,weight));
    localP=(1-weight)*first.pressure+weight*last.pressure;
    localP0=(1-weight)*first.totalPressure+weight*last.totalPressure;
    if any(~isfinite(localP) | ~isfinite(localP0) | ...
            localP<=0 | localP0<=localP)
        continue;
    end
    if thermalMode
        recovered=recoverThermallyPerfectEdgeState( ...
            localP,localP0,first.totalTemperature,R);
        result.valid(ids)=recovered.valid;
        result.temperature(ids)=recovered.temperature;
        result.density(ids)=recovered.density;
        result.speed(ids)=recovered.speed;
        result.Mach(ids)=recovered.Mach;
        result.totalTemperature(ids)=first.totalTemperature;
    else
        localM=sqrt(2/(g-1)*expm1((g-1)/g*log(localP0./localP)));
        localT=freestream.totalTemperature./(1+0.5*(g-1)*localM.^2);
        result.valid(ids)=true;
        result.Mach(ids)=localM;
        result.temperature(ids)=localT;
        result.density(ids)=localP./(R*localT);
        result.speed(ids)=localM.*sqrt(g*R*localT);
        result.totalTemperature(ids)=freestream.totalTemperature;
    end
    result.pressure(ids)=localP;
    result.totalPressure(ids)=localP0;
    if hasShockAngle
        result.shockAngleDeg(ids)= ...
            (1-weight)*first.shockAngleDeg+weight*last.shockAngleDeg;
    end
end
end
