function h=solvePhaseAwareNoseHeating(mesh,e,Tw)
%SOLVEPHASEAWARENOSEHEATING 偏置驻点拟合球鼻的相态层流热候选。
% 模块用途：以已有内部球鼻外缘、完整拟合短弧及Rsin(theta)度量，
% 求固定2度节点的局部相似剖面，再映射到原面。真实三维/非相似项未解。
% 原零攻角NASA7模块保持；不把本旋转球面近似称为全箭轴对称解。
%% 1. 原面、气体、工况及来源绑定；不接受改标签的旧场
n=size(mesh.F,1);Tw=double(Tw(:));if isscalar(Tw),Tw=repmat(Tw,n,1);end
assert(numel(Tw)==n&&all(isfinite(Tw)),'RAH:PhaseNoseHeat:Wall','给定壁温尺寸或数值错误。');
assert(strcmp(e.meshSHA256,phaseAwareMeshIdentity(mesh)),'RAH:PhaseNoseHeat:Geometry','原面几何改变。');
digest=phaseAwarePayloadSHA(e.sourceText,e.pressurePa,e.temperatureK,e.speedMps,e.velocityMps, ...
    e.recoveryTemperatureK,e.validFace,e.reasonCode,e.arcLengthM,e.Mach,e.faceAreaM2);
assert(strcmp(digest,e.source.sha256),'RAH:PhaseNoseHeat:Source','外缘数组来源不符。');
bound=jsondecode(e.sourceText);
assert(strcmp(e.source.role,'engineering-candidate')&&strcmp(e.source.id,'internal-phase-aware-spherical-nose-v1') ...
    &&strcmp(e.model,'rah-phase-aware-dry-air-of13-v1') ...
    &&strcmp(e.transportModel,'rah-lemmon-dilute-icoTabulated-of13-v1') ...
    &&strcmp(e.energyReference,'hf=0; hs(298.15K)=0; es=hs-RT') ...
    &&isequaln(bound.condition,jsondecode(jsonencode(e.condition)))&&isequaln(bound.stagnation,e.stagnation) ...
    &&isequal(bound.sphereCenter(:),e.sphereCenter(:)),'RAH:PhaseNoseHeat:Identity','气体/工况/驻点来源不一致。');
names=fieldnames(bound.codeSHA256);
for j=1:numel(names)
    name=strrep(names{j},'_m','.m');
    assert(strcmp(fileSHA(name),bound.codeSHA256.(names{j})),'RAH:PhaseNoseHeat:Code','外缘生产代码已变化，须重算当前工况。');
end
R=max(mesh.V(:,1))-e.sphereCenter(1);valid=logical(e.validFace(:));
assert(any(valid),'RAH:PhaseNoseHeat:Empty','当前没有可供热剖面使用的球鼻外缘。');
C=(mesh.V(mesh.F(:,1),:)+mesh.V(mesh.F(:,2),:)+mesh.V(mesh.F(:,3),:))/3;
rad=C-e.sphereCenter(:)';rad=rad./vecnorm(rad,2,2);up=-e.condition.flowDirection(:);
theta=e.arcLengthM/R;
assert(all(isfinite(theta(valid)))&&max(abs(cos(theta(valid))-rad(valid,:)*up))<1e-10, ...
    'RAH:PhaseNoseHeat:Path','弧长与当前驻点/原面不一致。');
% 重新核对有效面完整球帽；不重建保存文件没有提供的无效面标志。
est=estimateNoseRadius(mesh);cap=max(0,1-est.outerDepthRadiusFraction);
v=mesh.V-e.sphereCenter(:)';vr=vecnorm(v,2,2);f=mesh.F(valid,:);
g=mesh.localGeometry;
assert(est.available&&abs(est.radius/R-1)<1e-10&&up(1)>=cap-1e-10 ...
    &&all(all(abs(vr(f)/R-1)<=.03))&&all(all(v(f,1)>=R*cap-1e-8)) ...
    &&~any(g.featureFace(valid)|g.topologyIssue(valid)|~g.isSidewall(valid)), ...
    'RAH:PhaseNoseHeat:Path','当前原面不满足完整拟合球帽/特征排除条件。');
assert(all(theta(valid)<=deg2rad(e.condition.maximumPolarAngleDeg)+1e-12));
wallValues=unique(Tw(valid));assert(numel(wallValues)<=16,'RAH:PhaseNoseHeat:Wall','本版最多16种给定壁温。');
%% 2. 固定角表/同熵外缘与Levy–Lees轴对称尺度
tmax=max(theta(valid));assert(tmax>0,'RAH:PhaseNoseHeat:Path','仅驻点单面缺少有限弧段。');
angles=unique([0:deg2rad(2):tmax tmax]);
states=edgeStates(angles);tables=cell(numel(wallValues),numel(angles));
% 32点Gauss-Legendre直接积分每一驻点弧；不是从首个面零厚度启动。
k=(1:31)';b=k./sqrt(4*k.^2-1);[Q,L]=eig(diag(b,1)+diag(b,-1));
[gx,order]=sort(diag(L));gw=2*Q(1,order).^2;gx=(gx+1)/2;gw=gw/2;
quadAngles=angles(:)*gx(:)';quad=edgeStates(quadAngles);
xi=R^3*angles(:).*sum(quad.densityKgm3.*quad.dynamicViscosityPas.*quad.U.*sin(quadAngles).^2.*gw,2);
stations=struct([]);edges=struct([]);
for j=1:numel(angles)
    th=angles(j);U=states.U(j);metric=R*sin(th);
    if th==0
        du=e.stagnation.velocityGradientPerSecond;dp=0;beta=.5;
    else
        dp=-(e.stagnation.pressurePa-e.condition.pressurePa)*sin(2*th)/R;
        du=-dp/(states.densityKgm3(j)*U);
        beta=2*xi(j)*du/(states.densityKgm3(j)*states.dynamicViscosityPas(j)*U^2*metric^2);
    end
    stations(j)=struct('thetaRad',th,'arcM',R*th,'metricM',metric,'xi',xi(j),'beta',beta, ...
        'pressureGradientPam',dp,'velocityGradientPerSecond',du, ...
        'pressureCurvaturePam2',e.stagnation.pressureCurvaturePam2);
    edges(j)=struct('pressurePa',states.pressurePa(j),'temperatureK',states.temperatureK(j), ...
        'speedMps',U,'totalEnthalpyJkg',e.stagnation.enthalpyJkg);
    for w=1:numel(wallValues)
        r=solvePhaseAwareLaminarProfile(edges(j),stations(j),wallValues(w));
        r.Re_s=states.densityKgm3(j)*U*R*th/states.dynamicViscosityPas(j);
        if r.available&&(r.Re_s>=5e5||r.delta/R>=.05)
            r.available=false;r.reason='Re_s或delta/R超出本版层流薄层筛选范围';
        end
        tables{w,j}=r;
    end
end
%% 3. 相邻已解节点映射；所有原面仍计入面积分母
blank=nan(n,1);h=struct('schema','rah-phase-nose-heat-v1','validFace',false(n,1), ...
    'heatFlux',blank,'wallShearStress',blank,'delta',blank,'deltaStar',blank, ...
    'thetaPhysical',blank,'massThickness',blank,'skinFriction',blank, ...
    'lowerStation',zeros(n,1),'upperStation',zeros(n,1),'upperWeight',blank, ...
    'reason',{repmat({'尚无适用球鼻外缘；后体/横流/底面仍待求解'},n,1)}, ...
    'wallTemperature',Tw,'wallValues',wallValues,'anglesRad',angles, ...
    'stations',stations,'edges',edges,'profiles',{tables},'meshSHA256',e.meshSHA256, ...
    'edgeSourceSHA256',e.source.sha256,'physicalQualification',false,'newQualifiedAreaM2',0, ...
    'method','相态变Pr Levy–Lees局部相似BVP；偏置拟合球鼻旋转近似；相邻2度节点插值', ...
    'scope','原完整球帽短弧/有利压梯附着层流；忽略非相似及真实三维耦合；Taw仍为sqrt(Pr)估计；未获实验/接口资格');
fields={'heatFlux','wallShearStress','delta','deltaStar','thetaPhysical','massThickness'};
for id=find(valid)'
    j=find(angles<=theta(id),1,'last');if j==numel(angles),j=j-1;end
    weight=(theta(id)-angles(j))/(angles(j+1)-angles(j));w=find(wallValues==Tw(id));
    a=tables{w,j};b=tables{w,j+1};h.lowerStation(id)=j;h.upperStation(id)=j+1;h.upperWeight(id)=weight;
    if ~a.available||~b.available
        h.reason{id}=['相邻剖面拒绝：' a.reason ' / ' b.reason];continue
    end
    for key=fields,h.(key{1})(id)=(1-weight)*a.(key{1})+weight*b.(key{1});end
    rho=e.pressurePa(id)/(states.gasConstantJkgK*e.temperatureK(id));
    if e.speedMps(id)>0,h.skinFriction(id)=h.wallShearStress(id)/(.5*rho*e.speedMps(id)^2);end
    h.validFace(id)=true;h.reason{id}='已算工程候选；相邻节点插值，未作角向/实验验证';
end
h.summary=struct('originalFaces',n,'totalAreaM2',sum(e.faceAreaM2), ...
    'heatCandidateFaces',nnz(h.validFace),'heatCandidateAreaM2',sum(e.faceAreaM2(h.validFace)), ...
    'heatCandidateAreaFraction',sum(e.faceAreaM2(h.validFace))/sum(e.faceAreaM2), ...
    'signedPowerW',sum(h.heatFlux(h.validFace).*e.faceAreaM2(h.validFace)), ...
    'profileCount',numel(tables),'acceptedProfiles',sum(cellfun(@(r)r.available,tables),'all'));
h.sourceText=jsonencode(struct('edgeSHA256',e.source.sha256,'method',h.method,'scope',h.scope, ...
    'profileCodeSHA256',fileSHA('solvePhaseAwareLaminarProfile.m'),'regionCodeSHA256',fileSHA('solvePhaseAwareNoseHeating.m'), ...
    'profilePayloadSHA256',phaseAwarePayloadSHA(jsonencode(tables)), ...
    'meshSHA256',e.meshSHA256,'stations',stations,'wallValues',wallValues));
h.sourceSHA256=phaseAwarePayloadSHA(h.sourceText,h.heatFlux,h.wallShearStress,h.delta,h.deltaStar, ...
    h.thetaPhysical,h.massThickness,h.skinFriction,h.validFace,h.wallTemperature, ...
    h.lowerStation,h.upperStation,h.upperWeight);
    function s=edgeStates(th)
        p=e.condition.pressurePa+(e.stagnation.pressurePa-e.condition.pressurePa)*cos(th).^2;
        c0=phaseAwareAirCaloric(e.stagnation.temperatureK);target=c0.s0+c0.R*log(p/e.stagnation.pressurePa);
        lo=60*ones(size(th));hi=e.stagnation.temperatureK*ones(size(th));
        for it=1:50
            mid=(lo+hi)/2;c=phaseAwareAirCaloric(mid);below=c.s0<target;lo(below)=mid(below);hi(~below)=mid(~below);
        end
        T=(lo+hi)/2;T(th==0)=e.stagnation.temperatureK;s=phaseAwareAirState('PT',p,T,'linearTable');
        dh=e.stagnation.enthalpyJkg-s.enthalpyJkg;dh(th==0)=0;
        assert(all(dh(:)>=0),'RAH:PhaseNoseHeat:Energy','沿弧动能为负，不裁剪。');s.U=sqrt(2*dh);
    end
end
function h=fileSHA(name)
fid=fopen(fullfile(fileparts(mfilename('fullpath')),name),'rb');assert(fid>=0);c=onCleanup(@()fclose(fid));
b=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(b,'int8'));
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
