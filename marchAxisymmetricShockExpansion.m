function mesh = marchAxisymmetricShockExpansion(mesh,freestream,cfg)
%MARCHAXISYMMETRICSHOCKEXPANSION 零攻角近轴对称箭体的逐段波系估算。
% 模块作用：在已求得的鼻锥切锥状态后，沿轴向截面从鼻到尾推进。
% 壁面向外流转开时用 Prandtl-Meyer 膨胀，向外流转回时用二维
% 附体弱斜激波；直段携带上游状态，不再恢复自由来流总压。
% 默认 NASA7 模式下，锥流锚点、膨胀与压缩均使用同一固定组分
% 干空气焓/熵模型；仍是代表流管近似，并非轴对称 Euler 全场。
% 这是假设各截面共享同一条代表流管的低阶 shock-expansion 模型，
% 不是轴对称 Euler/CFD 解，也不解析边界层、脱体波或波系干扰。

%% 1. 校验适用性；不满足条件时保留原逐面元外流
if nargin<3, cfg=struct; end
status=struct('available',false,'applied',false, ...
    'method','直鼻锥锚定 + 轴向逐段 Prandtl-Meyer/斜激波推进', ...
    'scope',['仅零攻角、近轴对称、带可识别直鼻锥的超声速箭体；' ...
        '二维局部波关系沿轴向携带总压，不能描述三维波系相互作用或分离'], ...
    'status','尚未检查几何与工况');
enabled=true;
if isfield(cfg,'externalFlow') && ...
        isfield(cfg.externalFlow,'axisymmetricMarchEnabled')
    enabled=cfg.externalFlow.axisymmetricMarchEnabled;
end
validateattributes(enabled,{'logical','numeric'},{'binary','scalar'});
if ~enabled
    status.status='配置已关闭轴向激波-膨胀推进';
    mesh.axisymmetricMarch=status; return;
end
if ~isfield(mesh,'boundaryLayerEdge') || ~isfield(mesh,'externalFlow') || ...
        ~isfield(mesh,'localGeometry') || ~isfield(mesh,'param') || ...
        ~isfield(mesh.param,'sectionX') || ...
        ~isfield(mesh.param,'sectionEquivalentRadius') || ...
        ~isfield(mesh.param,'faceBandID')
    status.status='缺少内置外流或原生轴向截面';
    mesh.axisymmetricMarch=status; return;
end
if hypot(freestream.alpha,freestream.beta)>1e-8 || freestream.Mach<=1
    status.status='轴向模型只用于零攻角、零侧滑的超声速工况';
    mesh.axisymmetricMarch=status; return;
end
if isfield(mesh,'axisymmetricMarch') && mesh.axisymmetricMarch.applied
    return;
end
edge=mesh.boundaryLayerEdge;
if ~isfield(edge,'faceMethodID') || ...
        ~contains(edge.method,'内置工程外流') || ...
        numel(edge.faceMethodID)~=size(mesh.F,1)
    status.status='需要与当前网格一致的内置局部外流';
    mesh.axisymmetricMarch=status; return;
end
x=double(mesh.param.sectionX(:));
r=double(mesh.param.sectionEquivalentRadius(:));
band=double(mesh.param.faceBandID(:));
nb=numel(x)-1; nf=size(mesh.F,1);
if numel(r)~=numel(x) || numel(band)~=nf || nb<6 || ...
        any(~isfinite(x) | ~isfinite(r)) || any(diff(x)<=0) || ...
        any(r<=0) || any(band~=floor(band))
    status.status='轴向截面半径或面元环带编号无效';
    mesh.axisymmetricMarch=status; return;
end
if isfield(mesh.param,'sectionCenterYZ')
    centers=double(mesh.param.sectionCenterYZ);
    if size(centers,1)~=numel(x) || ...
            max(vecnorm(centers-mean(centers,1),2,2))>0.05*max(r)
        status.status='截面中心线偏离直轴，不适用近轴对称推进';
        mesh.axisymmetricMarch=status; return;
    end
end
if isfield(mesh.param,'sectionPerimeter')
    perimeter=double(mesh.param.sectionPerimeter(:));
    circularity=perimeter./(2*pi*r);
    central=x>max(x)-0.2*(max(x)-min(x));
    if numel(perimeter)~=numel(x) || ...
            median(circularity(central),'omitnan')>1.08
        status.status='前体截面偏离近圆形，不适用轴对称推进';
        mesh.axisymmetricMarch=status; return;
    end
end

%% 2. 识别直鼻锥平台，取其最后一段作为上游锚点
phi=atan2d(r(1:end-1)-r(2:end),diff(x));
bodyLength=max(x)-min(x);
midX=0.5*(x(1:end-1)+x(2:end));
candidate=midX>max(x)-0.2*bodyLength & phi>=5 & phi<=60;
if nnz(candidate)<3
    status.status='前体没有足够的附体直鼻锥截面';
    mesh.axisymmetricMarch=status; return;
end
rounded=round(phi(candidate)*2)/2;
modeAngle=mode(rounded);
plateau=candidate & abs(phi-modeAngle)<=0.6;
edges=diff([false;plateau;false]);
starts=find(edges==1); ends=find(edges==-1)-1;
lengths=x(ends+1)-x(starts);
[longest,run]=max(lengths);
if isempty(run) || (ends(run)-starts(run)+1)<3 || ...
        longest<0.02*bodyLength
    status.status='没有达到长度阈值的直鼻锥平台';
    mesh.axisymmetricMarch=status; return;
end
anchorBand=starts(run);
anchorAngle=median(phi(anchorBand:ends(run)));
[airMode,~]=selectAirThermalModel(freestream,cfg);
thermalMarch=false;
if strcmp(airMode,'thermallyPerfect')
    thermalReference=solveThermallyPerfectConicalShock(freestream,0);
    thermalMarch=thermalReference.available;
end
if thermalMarch
    state=solveThermallyPerfectConicalShock(freestream,anchorAngle);
    status.method=['NASA7 变比热 Taylor-Maccoll 切锥锚点 + ' ...
        'NASA7 二维 Prandtl-Meyer/附体斜激波轴向推进'];
    status.thermalPropertyMode='thermallyPerfect';
else
    state=solveConicalShock(freestream,anchorAngle);
    status.thermalPropertyMode='caloricallyPerfect';
    if strcmp(airMode,'thermallyPerfect')
        status.method=[status.method '（热物性域外，定比热回退）'];
    end
end
if ~state.attached || ~isfinite(state.totalPressure)
    status.status='识别出的直鼻锥没有附体切锥解';
    mesh.axisymmetricMarch=status; return;
end
anchorTotalPressure=state.totalPressure;

%% 3. 从直鼻锥向尾部逐段推进；波后总压不能无故回升
waveType=zeros(nb,1,'uint8');
bandPressure=nan(nb,1); bandTotalPressure=bandPressure;
bandMach=bandPressure; bandShockAngle=bandPressure;
faceWaveType=zeros(nf,1,'uint8');
faceWaveShockAngle=nan(nf,1);
marchedFace=false(nf,1);
unresolvedBand=false(nb,1);
referenceAngle=anchorAngle;
failed=false;
thermalDomainFailure=false;
thermalFailureStatus='';
for j=anchorBand-1:-1:1
    if ~failed
        turn=referenceAngle-phi(j);
        wave=uint8(1); % 1=沿直段携带上游状态。
        next=struct('status','');
        if abs(turn)>=0.5
            % 变比热求解器的接口要求输入 Mach 与输入 gamma 对应；
            % 内部再用 NASA7 局部声速重新计算真实热力 Mach。
            upstream=struct('Mach',state.speed/sqrt( ...
                freestream.gamma*freestream.R*state.temperature), ...
                'gamma',freestream.gamma, ...
                'R',freestream.R,'pressure',state.pressure, ...
                'temperature',state.temperature,'density',state.density, ...
                'speed',state.speed,'totalPressure',state.totalPressure, ...
                'totalTemperature',state.totalTemperature);
            if turn>0 && turn<90 && upstream.Mach>1
                if thermalMarch
                    next=solveThermallyPerfectPrandtlMeyerExpansion( ...
                        upstream,turn);
                else
                    next=solvePrandtlMeyerExpansion(upstream,turn);
                end
                success=next.available;
                wave=uint8(2); % 2=凸转角膨胀。
            elseif turn<0 && -turn<90 && upstream.Mach>1
                if thermalMarch
                    next=solveThermallyPerfectObliqueShock(upstream,-turn);
                    success=next.available && next.attached;
                else
                    next=solveObliqueShock(upstream,-turn);
                    success=next.attached;
                end
                wave=uint8(3); % 3=压缩转角附体斜激波。
            else
                success=false;
            end
            if success && isfinite(next.totalPressure) && ...
                    next.totalPressure<=state.totalPressure*(1+1e-8) && ...
                    next.totalPressure>next.pressure
                if wave==2
                    % 等熵膨胀总压严格携带上游值；求根与插值产生的
                    % 微小舍入差不应在长距离推进中累积为总压回升。
                    next.totalPressure=state.totalPressure;
                    next.totalTemperature=state.totalTemperature;
                end
                state=next;
                referenceAngle=phi(j);
            else
                failed=true;
                status.status=sprintf('第 %d 轴向环带的局部转角无附体/等熵解',j);
                if thermalMarch && exist('next','var') && ...
                        isfield(next,'status') && ...
                        (contains(next.status,'80 K') || ...
                        contains(next.status,'热物性'))
                    thermalDomainFailure=true;
                    thermalFailureStatus=sprintf('第 %d 环带：%s',j,next.status);
                end
            end
        end
    end
    if failed
        unresolvedBand(j)=true;
    else
        waveType(j)=wave;
        bandPressure(j)=state.pressure;
        bandTotalPressure(j)=state.totalPressure;
        bandMach(j)=state.Mach;
        if wave==3, bandShockAngle(j)=state.shockAngleDeg; end
    end
    ids=find(band==j & mesh.localGeometry.isSidewall & ...
        ~mesh.localGeometry.featureFace & ...
        ~mesh.localGeometry.topologyIssue);
    if isempty(ids), continue; end
    if failed
        edge.faceMethodID(ids)=0;
        for name={'pressure','temperature','density','speed','Mach', ...
                'totalPressure','totalTemperature','shockAngleDeg', ...
                'speedOfSound','dynamicViscosity'}
            edge.(name{1})(ids)=NaN;
        end
        edge.velocityVector(ids,:)=NaN;
        edge.algebraicallyValid(ids)=false;
        edge.physicallyUnresolved(ids)=true;
        if isfield(edge,'thermalPropertyModeByFace')
            edge.thermalPropertyModeByFace(ids)=0;
        end
        continue;
    end
    edge.pressure(ids)=state.pressure;
    edge.temperature(ids)=state.temperature;
    edge.density(ids)=state.density;
    edge.speed(ids)=state.speed;
    edge.Mach(ids)=state.Mach;
    edge.totalPressure(ids)=state.totalPressure;
    edge.totalTemperature(ids)=state.totalTemperature;
    if thermalMarch
        edge.speedOfSound(ids)=state.speedOfSound;
    else
        edge.speedOfSound(ids)=sqrt( ...
            freestream.gamma*freestream.R*state.temperature);
    end
    edge.dynamicViscosity(ids)=1.458e-6*state.temperature^(3/2)/ ...
        (state.temperature+110.4);
    edge.shockAngleDeg(ids)=bandShockAngle(j);
    edge.faceMethodID(ids)=6;
    if isfield(edge,'thermalPropertyModeByFace')
        if thermalMarch
            code=uint8(2+double(state.temperature<150));
        else
            code=uint8(1);
        end
        edge.thermalPropertyModeByFace(ids)=code;
    end
    edge.algebraicallyValid(ids)=true;
    edge.physicallyUnresolved(ids)=false;
    faceWaveType(ids)=wave;
    faceWaveShockAngle(ids)=bandShockAngle(j);
    marchedFace(ids)=true;
end

% 不能在同一轴向流管上把域外定比热状态接到 NASA7 锥流状态后面。
% 若唯一障碍是热物性域，重算整条定比热外流，保留可用的工程热流；
% 结果明确标记为全链路定比热回退，不能冒充变比热解。
if thermalDomainFailure
    fallbackCfg=cfg;
    if ~isfield(fallbackCfg,'heating'), fallbackCfg.heating=struct; end
    fallbackCfg.heating.airThermalModel='caloricallyPerfect';
    fallback=solveExternalFlow(mesh,freestream,fallbackCfg);
    fallback=marchAxisymmetricShockExpansion(fallback,freestream,fallbackCfg);
    % 鼻端驻点法向激波独立于侧壁轴向流管，可继续保留 NASA7 解。
    if isfield(mesh,'stagnationShock')
        fallback.stagnationShock=mesh.stagnationShock;
        fallback.externalFlow.stagnationShock=mesh.stagnationShock;
    end
    fallback.axisymmetricMarch.requestedThermalPropertyMode='thermallyPerfect';
    fallback.axisymmetricMarch.thermalFallbackReason=thermalFailureStatus;
    fallback.axisymmetricMarch.status=[fallback.axisymmetricMarch.status ...
        '；干空气低温近似超出 80 K 下限，整条外流退回定比热'];
    fallback.axisymmetricMarch.method=[fallback.axisymmetricMarch.method ...
        '（NASA7 域外时整条链路回退）'];
    fallback.boundaryLayerEdge.method=[fallback.boundaryLayerEdge.method ...
        '；NASA7 域外，整条零攻角外流退回定比热'];
    fallback.externalFlow.edge=fallback.boundaryLayerEdge;
    fallback.externalFlow.method=fallback.boundaryLayerEdge.method;
    mesh=fallback;
    return;
end

%% 4. 更新三维切向速度与逐面元方法来源
if any(marchedFace)
    N=mesh.localGeometry.normal;
    flow=freestream.unitFlowDirection(:)';
    tangent=repmat(flow,nf,1)-sum(N.*flow,2).*N;
    tangentLength=vecnorm(tangent,2,2);
    usable=marchedFace & tangentLength>1e-10;
    edge.velocityVector(usable,:)= ...
        edge.speed(usable).*tangent(usable,:)./tangentLength(usable);
    unusable=marchedFace & ~usable;
    edge.velocityVector(unusable,:)=NaN;
    edge.faceMethodID(unusable)=0;
    if isfield(edge,'thermalPropertyModeByFace')
        edge.thermalPropertyModeByFace(unusable)=0;
    end
    edge.algebraicallyValid(unusable)=false;
    edge.physicallyUnresolved(unusable)=true;
    marchedFace(unusable)=false;
end
edge.method=[edge.method ' + 零攻角轴向激波-膨胀推进：' status.method];
edge.scope=[edge.scope '；零攻角近轴对称直鼻锥后沿轴向携带总压，' ...
    '以局部 Prandtl-Meyer 膨胀/二维附体斜激波处理截面转角；' ...
    '这是代表流管工程估算，不是连续激波层流场'];
edge.shockAngleMethod=[edge.shockAngleMethod ...
    '；轴向压缩转角的二维 theta-beta-M 附体弱激波角'];
edge.shockAngleScope=[edge.shockAngleScope ...
    '；轴向波系编号 6 的角度仅在附体压缩转角非 NaN，' ...
    '不代表三维激波面位置'];
edge.faceMethodNames{7}=status.method;
edge.axialWaveType=faceWaveType;
edge.axialWaveTypeNames={'未推进','携带上游状态', ...
    'Prandtl-Meyer 膨胀','二维附体斜激波'};
edge.axialWaveShockAngleDeg=faceWaveShockAngle;
edge.axialMarchFace=marchedFace;
edge.summary.validFaces=nnz(edge.algebraicallyValid);
edge.summary.unresolvedFaces=nf-edge.summary.validFaces;
edge.summary.grazingFaces=nnz(edge.faceMethodID==1);
edge.summary.coneFaces=nnz(edge.faceMethodID==2);
edge.summary.wedgeFaces=nnz(edge.faceMethodID==3);
edge.summary.expansionFaces=nnz(edge.faceMethodID==4);
edge.summary.detachedClosureFaces=nnz(edge.faceMethodID==5);
edge.summary.axialMarchFaces=nnz(marchedFace);
if isfield(edge,'thermalPropertyModeByFace')
    edge.summary.lowTemperatureThermalFaces= ...
        nnz(edge.thermalPropertyModeByFace==3);
end
mesh.boundaryLayerEdge=edge;
mesh.externalFlow.edge=edge;
mesh.externalFlow.method=edge.method;

status.available=true;
status.applied=any(marchedFace);
if ~failed
    status.status=sprintf('从直鼻锥后向尾部推进 %d 个面元',nnz(marchedFace));
end
status.anchorBand=anchorBand;
status.anchorAngleDeg=anchorAngle;
status.anchorTotalPressure=anchorTotalPressure;
status.anchorTotalTemperature=state.totalTemperature;
status.sectionX=x;
status.sectionRadius=r;
status.bandAngleDeg=phi;
status.bandWaveType=waveType;
status.bandPressure=bandPressure;
status.bandTotalPressure=bandTotalPressure;
status.bandMach=bandMach;
status.bandShockAngleDeg=bandShockAngle;
status.unresolvedBand=unresolvedBand;
status.faceMask=marchedFace;
status.summary=struct('marchedFaces',nnz(marchedFace), ...
    'expansionBands',nnz(waveType==2), ...
    'compressionBands',nnz(waveType==3), ...
    'unresolvedBands',nnz(unresolvedBand));
mesh.axisymmetricMarch=status;
for name={'surfaceStreamlines','boundaryLayer','heating','heatingCoverage', ...
        'wallThermal','aerothermalExposure','flightHeatingHistory'}
    if isfield(mesh,name{1}), mesh=rmfield(mesh,name{1}); end
end
end
