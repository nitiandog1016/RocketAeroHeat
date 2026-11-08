function mesh=calcPhaseAwareLaminarHeating(mesh,cfg)
%CALCPHASEAWARELAMINARHEATING 同模型外缘—参考焓—附着层流积分—热流。
% 模块用途：显式选择的新混合工程链；读取mesh.phaseAwareEdge，不重建外流。
% 保留原Thwaites型参考性质积分与Reynolds比拟，不宣称解出完整三维边界层。
% 旧NASA7/Sutherland默认及已有工程结果不会被静默用作相态模型输入。

%% 1. 绑定气体、原面、来源和严格上游路径；禁止标量掩码隐式广播
assert(isfield(mesh,'phaseAwareEdge'),'RAH:PhaseLayer:Input','缺少具名外缘包。');
e=mesh.phaseAwareEdge; n=size(mesh.F,1); identity=phaseAwareMeshIdentity(mesh);
required={'schema','model','transportModel','energyReference','meshSHA256','originalFaceID', ...
    'pressurePa','temperatureK','velocityMps','source','validFace','attachedFace','laminarFace', ...
    'quasiTwoDimensionalFace','continuumLengthM','freestreamMeanFreePathM','path'};
assert(all(isfield(e,required)),'RAH:PhaseLayer:Contract','外缘包缺少必要物理/来源字段。');
assert(strcmp(e.schema,'rah-phase-engineering-edge-v1') && ...
    strcmp(e.model,'rah-phase-aware-dry-air-of13-v1') && ...
    strcmp(e.transportModel,'rah-lemmon-dilute-icoTabulated-of13-v1') && ...
    strcmp(e.energyReference,'hf=0; hs(298.15K)=0; es=hs-RT'), ...
    'RAH:PhaseLayer:Identity','气体/输运/能量身份不相容，禁止混入内部+300000坐标。');
assert(strcmp(identity,e.meshSHA256) && isequal(e.originalFaceID(:),(1:n)'), ...
    'RAH:PhaseLayer:Geometry','原面编号、顺序或几何不匹配。');
assert(isstruct(e.source)&&all(isfield(e.source,{'id','sha256','role','method'})) && ...
    ischar(e.source.id)&&~isempty(e.source.id)&&ischar(e.source.method)&&~isempty(e.source.method) && ...
    ~isempty(regexp(e.source.sha256,'^[0-9a-f]{64}$','once')), ...
    'RAH:PhaseLayer:Source','来源必须明确且含SHA。');
assert(ismember(e.source.role,{'prescribed-control','accepted-edge','unqualified-research'}), ...
    'RAH:PhaseLayer:Source','未知来源角色。');
for name={'pressurePa','temperatureK','continuumLengthM','freestreamMeanFreePathM'}
    validateattributes(e.(name{1}),{'numeric'},{'real','column','numel',n});
end
validateattributes(e.velocityMps,{'numeric'},{'real','size',[n 3]});
for name={'validFace','attachedFace','laminarFace','quasiTwoDimensionalFace'}
    validateattributes(e.(name{1}),{'numeric','logical'},{'binary','column','numel',n});
end
p=e.path; required={'predecessorFace','stepLength','distance','pathCost','connectedFace','sourceSHA256','meshSHA256'};
assert(all(isfield(p,required)) && strcmp(p.sourceSHA256,e.source.sha256) && strcmp(p.meshSHA256,identity), ...
    'RAH:PhaseLayer:Path','路径没有绑定同一外缘与几何。');
for name=required(1:4),validateattributes(p.(name{1}),{'numeric'},{'real','column','numel',n});end
validateattributes(p.connectedFace,{'numeric','logical'},{'binary','column','numel',n});
parent=p.predecessorFace;connected=logical(p.connectedFace);
assert(all(isfinite(parent)&parent==fix(parent)&parent>=0&parent<=n), ...
    'RAH:PhaseLayer:Path','前驱编号无效。');
assert(all(isfinite(p.pathCost(connected))) && all(isfinite(p.distance(connected))&p.distance(connected)>=0), ...
    'RAH:PhaseLayer:Path','连通路径成本/距离无效。');
centers=(mesh.V(mesh.F(:,1),:)+mesh.V(mesh.F(:,2),:)+mesh.V(mesh.F(:,3),:))/3;
for id=find(connected)'
    j=parent(id);
    if j==0
        assert(p.distance(id)==0&&p.stepLength(id)==0,'RAH:PhaseLayer:Origin','没有非零入口厚度时，种子必须为真实零距离起点。');
    else
        assert(connected(j)&&p.pathCost(j)<p.pathCost(id)&&isfinite(p.stepLength(id))&&p.stepLength(id)>0 && ...
            abs(p.distance(id)-p.distance(j)-p.stepLength(id))<1e-10*max(1,p.distance(id)), ...
            'RAH:PhaseLayer:Path','前驱必须严格在上游，物理步长必须闭合。');
        delta=centers(id,:)-centers(j,:);
        assert(p.stepLength(id)>=norm(delta)-1e-10*max(1,p.stepLength(id)), ...
            'RAH:PhaseLayer:Path','弧长不能短于原面中心直线距离。');
        if e.validFace(id)&&e.validFace(j)&&all(isfinite(e.velocityMps([id j],:)),'all')
            assert(dot(e.velocityMps(id,:),delta)>0&&dot(e.velocityMps(j,:),delta)>0, ...
                'RAH:PhaseLayer:Direction','前驱链与实际三维流向相反，不能沿反向路径积分。');
        end
    end
end
sourceAccepted=true(n,1);
if strcmp(e.source.role,'unqualified-research'),sourceAccepted(:)=false;end
if strcmp(e.source.role,'accepted-edge')
    % 数值入口审核单独绑定；即使此记录通过，仍不自动授热流/实验资格。
    assert(isfield(e,'acceptance')&&all(isfield(e.acceptance,{'sourceSHA256','meshSHA256','acceptedFace','method','role'})), ...
        'RAH:PhaseLayer:Acceptance','缺少外缘接受记录。');
    a=e.acceptance;
    assert(strcmp(a.sourceSHA256,e.source.sha256)&&strcmp(a.meshSHA256,identity)&& ...
        strcmp(a.role,'engineering-edge-acceptance')&&ischar(a.method)&&~isempty(a.method), ...
        'RAH:PhaseLayer:Acceptance','接受记录与源不相符。');
    validateattributes(a.acceptedFace,{'numeric','logical'},{'binary','column','numel',n});
    sourceAccepted=logical(a.acceptedFace);
end
Tw=surfaceWallTemperature(n,cfg); U=sqrt(sum(e.velocityMps.^2,2));
area=sqrt(sum(cross(mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,1),:), ...
    mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,1),:),2).^2,2))/2;
assert(all(isfinite(area)&area>0),'RAH:PhaseLayer:Geometry','原面面积无效。');

%% 2. 分清直接拒绝原因；单个域外状态不能拖累其他独立流线
names={'层流工程候选','来源未接受','外缘/壁温缺失','非附着层流或强三维区域', ...
    '路径不连通','物性联合域外','连续介质筛选未过','上游受阻', ...
    '零距离起点热流未定义','Thwaites疑似分离','Thwaites压梯范围外','数值积分无效'};
code=zeros(n,1,'uint8');code(~sourceAccepted)=1;
stateOK=logical(e.validFace)&isfinite(U)&U>0 & isfinite(e.pressurePa)&e.pressurePa>0 & ...
    isfinite(e.temperatureK)&isfinite(Tw)&isfinite(e.continuumLengthM)&e.continuumLengthM>0 & ...
    isfinite(e.freestreamMeanFreePathM)&e.freestreamMeanFreePathM>0;
code(code==0&~stateOK)=2;
code(code==0&~(e.attachedFace&e.laminarFace&e.quasiTwoDimensionalFace))=3;
code(code==0&~connected)=4;
Taw=nan(n,1);Tref=Taw;href=Taw;drive=Taw;rho=Taw;mu=Taw;nu=Taw;Pr=Taw;PrEdge=Taw;kn=Taw;r=Taw;
referenceBatch(find(code==0));
code(code==0&kn>=.01)=6;

%% 3. 沿真实前驱链积分，任何拒绝/分离都传播，禁止下游重新从零起算
theta=nan(n,1);Cf=theta;lambda=theta;integralTerm=theta;gradient=theta;
[~,order]=sort(p.pathCost);
for id=order'
    if ~connected(id),continue;end
    j=parent(id);
    if j>0 && ~ismember(code(j),uint8([0 8]))
        if code(id)==0,code(id)=7;end
        continue;
    end
    if code(id)~=0,continue;end
    if j==0, integralTerm(id)=0;code(id)=8;continue;end
    ds=p.stepLength(id);
    integralTerm(id)=integralTerm(j)+.5*ds*(nu(j)*U(j)^5+nu(id)*U(id)^5);
    theta(id)=sqrt(.45*integralTerm(id)/U(id)^6);
    upstream=j;span=ds;
    for k=2:4
        older=parent(upstream);if older==0,break;end
        span=span+p.stepLength(upstream);upstream=older;
    end
    gradient(id)=(U(id)-U(upstream))/span;
    lambda(id)=theta(id)^2/nu(id)*gradient(id);
    if ~isfinite(theta(id))||theta(id)<=0||~isfinite(lambda(id)),code(id)=11;
    elseif lambda(id)<=-.09,code(id)=9;
    elseif lambda(id)>.1,code(id)=10;
    else
        Cf(id)=2*(.2227+1.57*lambda(id)-1.8*lambda(id)^2)/(U(id)*theta(id)/nu(id));
    end
end
valid=code==0&isfinite(Cf)&Cf>0;
q=nan(n,1);St=q;
St(valid)=Cf(valid)./(2*Pr(valid).^(2/3));
q(valid)=St(valid).*rho(valid).*U(valid).*drive(valid);
% q正负均保留；零热驱动产生物理零，而不是将未解值补零。
assert(all(isfinite(q(valid))),'RAH:PhaseLayer:Numerics','热流超出浮点数值范围。');

%% 4. 结果携带物理定义、面积及缺项；不从一个厚度伪造CFD入口剖面
method='相态标准Air/线性mu-k + Riley参考焓 + Thwaites型参考性质积分 + Reynolds比拟';
scope=['附着、层流、准二维工程近似；Kn<0.01为特征长度筛选而非梯度Kn；' ...
    '参考性质厚度不是压缩性质量加权积分；不覆盖强横流/分离再附/底面；未完成热实验资格'];
b=struct('method',method,'scope',scope,'model',e.model,'transportModel',e.transportModel, ...
    'source',e.source,'meshSHA256',identity,'originalFaceID',(1:n)', ...
    'pressurePa',e.pressurePa,'edgeTemperatureK',e.temperatureK,'velocityMps',e.velocityMps, ...
    'wallTemperatureK',Tw,'heatFluxWm2',q,'skinFrictionCoefficient',Cf,'StantonNumber',St, ...
    'referenceMomentumThicknessM',theta,'referenceKinematicViscosityM2s',nu, ...
    'referenceTemperatureK',Tref,'referenceEnthalpyJkg',href,'recoveryTemperatureK',Taw, ...
    'thermalDrivingEnthalpyJkg',drive,'referencePrandtl',Pr,'edgePrandtl',PrEdge, ...
    'recoveryFactor',r,'recoveryFactorMethod','sqrt(Pr_edge), laminar engineering approximation', ...
    'pressureGradientParameter',lambda,'velocityGradientPerSecond',gradient,'pathIntegral',integralTerm, ...
    'knudsen',kn,'reasonCode',code,'reasonNames',{names},'validFace',valid,'faceAreaM2',area, ...
    'compressibleMomentumThicknessM',nan(n,1),'displacementThicknessM',nan(n,1),'energyThicknessM',nan(n,1), ...
    'profile',[],'profileReady',false,'profileMissingReason','未求解或重建密度/速度/焓剖面，不能由Thwaites参考厚度唯一确定', ...
    'physicalQualification',false,'newQualifiedAreaM2',0, ...
    'inputAcceptanceEvidenceIndependentlyVerified',false, ...
    'summary',struct('totalAreaM2',sum(area),'candidateAreaM2',sum(area(valid)), ...
        'candidateFaces',nnz(valid),'heatPowerOverCandidateAreaW',sum(area(valid).*q(valid))));
mesh.phaseAwareLaminar=b;
mesh.heating=struct('sidewallModel','phaseAwareLaminarIntegral','method',method,'validity',scope, ...
    'heatFlux',q,'laminarHeatFlux',q,'turbulentHeatFlux',nan(n,1), ...
    'eligibleFace',valid,'unresolvedFace',~valid,'wallTemperatureByFace',Tw, ...
    'recoveryTemperature',Taw,'physicalQualification',false,'historyCompatible',false,'summary',b.summary);
mesh.heating.summary.calculatedFaces=nnz(valid);
mesh.heating.summary.unresolvedFaces=nnz(~valid);
for name={'heatingCoverage','wallThermal','aerothermalExposure','flightHeatingHistory'}
    if isfield(mesh,name{1}),mesh=rmfield(mesh,name{1});end
end

    function referenceBatch(ids)
        % 批量相态查询，只有已知物性域拒绝才二分隔离；程序错误不得吞掉。
        if isempty(ids),return;end
        try
            edge=phaseAwareAirState('PT',e.pressurePa(ids),e.temperatureK(ids),'linearTable');
            factor=sqrt(edge.Prandtl);
            ref=phaseAwareAirReferenceState(e.pressurePa(ids),e.temperatureK(ids),U(ids),Tw(ids),factor,'linearTable');
            Taw(ids)=ref.recoveryTemperatureK;Tref(ids)=ref.referenceTemperatureK;
            href(ids)=ref.referenceEnthalpyJkg;drive(ids)=ref.thermalDrivingEnthalpyJkg;
            rho(ids)=edge.densityKgm3;mu(ids)=ref.dynamicViscosityPas;
            nu(ids)=mu(ids)./ref.referenceDensityKgm3;Pr(ids)=ref.Prandtl;PrEdge(ids)=edge.Prandtl;r(ids)=factor;
            kn(ids)=max(edge.meanFreePathViscosityDefinitionM,e.freestreamMeanFreePathM(ids))./e.continuumLengthM(ids);
        catch err
            if ~ismember(err.identifier,{'RAH:PhaseAir:Domain','RAH:PhaseAir:Temperature', ...
                    'RAH:PhaseAir:Energy','RAH:PhaseAir:Inverse','RAH:PhaseAir:Numerics'}),rethrow(err);end
            if numel(ids)==1,code(ids)=5;
            else,k=floor(numel(ids)/2);referenceBatch(ids(1:k));referenceBatch(ids(k+1:end));end
        end
    end
end
