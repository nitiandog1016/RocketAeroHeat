function e=solvePhaseAwareNoseEdge(mesh,fs,cfg)
%SOLVEPHASEAWARENOSEEDGE 同相态空气、原几何球鼻的内部工程外缘生产器。
% 模块用途：NASA SP-8062 §4.2.1驻点熵近似＋修正牛顿压力，求原面p/Te/U。
% 仅已核验球鼻前缀，不将驻点熵外推到后体，不计算三维激波面或壁面热流。
% hs(298.15K)=0；同Cp/线性mu-k；所有状态最终经过p/T联合域审核。
if nargin<3,cfg=struct;end
V=double(mesh.V);F=double(mesh.F);n=size(F,1);g=mesh.localGeometry;
C=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
raw=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
area=vecnorm(raw,2,2)/2;assert(all(area>0));N=raw./(2*area);
assert(isequal(size(g.normal),[n 3])&&max(abs(g.normal-N),[],'all')<1e-10 && ...
    max(abs(g.center-C),[],'all')<1e-10,'RAH:PhaseNose:Geometry','局部几何与当前原面不一致。');
flow=fs.unitFlowDirection(:)';assert(numel(flow)==3&&all(isfinite(flow))&&abs(norm(flow)-1)<1e-10);
validateattributes(fs.Mach,{'numeric'},{'real','finite','scalar','>',1});
air=phaseAwareAirState('PT',fs.pressure,fs.temperature,'linearTable');speed=fs.Mach*air.soundSpeedMps;
limit=60;if isfield(cfg,'externalFlow')&&isfield(cfg.externalFlow,'sphericalNoseMaxPolarAngleDeg'),limit=cfg.externalFlow.sphericalNoseMaxPolarAngleDeg;end
validateattributes(limit,{'numeric'},{'real','finite','scalar','>',0,'<=',80});
estimate=estimateNoseRadius(mesh);blank=nan(n,1);
e=struct('schema','rah-phase-nose-edge-v1','model',air.model,'transportModel',air.transportModel, ...
    'energyReference',air.energyReference,'meshSHA256',phaseAwareMeshIdentity(mesh),'originalFaceID',(1:n)', ...
    'pressurePa',blank,'temperatureK',blank,'speedMps',blank,'velocityMps',nan(n,3),'Mach',blank, ...
    'recoveryTemperatureK',blank,'densityKgm3',blank,'entropyResidual',blank,'energyResidual',blank, ...
    'polarAngleRad',blank,'arcLengthM',blank,'transverseMetricM',blank,'analyticTangent',nan(n,3), ...
    'validFace',false(n,1),'geometryFace',false(n,1),'completeFittedArc',false(n,1), ...
    'reasonCode',ones(n,1,'uint8'),'reasonNames',{{'鼻部外缘工程候选','球鼻几何/工况不适用', ...
    '原面不在完整拟合球鼻内','极角/完整弧线范围不符','热物性/总态或连续流筛选拒绝','切向方向退化'}}, ...
    'faceAreaM2',area,'geometryEstimate',estimate,'sphereCenter',[NaN NaN NaN], ...
    'stagnationPoint',[NaN NaN NaN],'stagnation',struct,'normalShock',struct, ...
    'physicalQualification',false,'newQualifiedAreaM2',0,'heatReady',false,'profileReady',false, ...
    'scope','拟合球鼻前缀；合成攻角<=30度；修正牛顿/驻点熵工程近似；不解后体熵层/三维横流或热剖面', ...
    'method','相态标准Air正激波 + 球鼻修正牛顿压力 + 同熵/总焓恢复', ...
    'sourceURL','https://ntrs.nasa.gov/api/citations/19710021703/downloads/19710021703.pdf', ...
    'sourceLocation','NASA SP-8062 §4.2.1, printed p30 / PDF p36', ...
    'condition',struct('Mach',fs.Mach,'pressurePa',fs.pressure,'temperatureK',fs.temperature, ...
        'speedMps',speed,'flowDirection',flow,'maximumPolarAngleDeg',limit));
incidence=acosd(max(-1,min(1,-flow(1))));
if estimate.available&&incidence<=30
    R=estimate.radius;xTip=max(V(:,1));tol=max(1e-10,1e-8*(xTip-min(V(:,1))));
    tip=mean(V(abs(V(:,1)-xTip)<=tol,:),1);origin=[xTip-R tip(2:3)];e.sphereCenter=origin;
    upstream=-flow;e.stagnationPoint=origin+R*upstream;
    radial=C-origin;radial=radial./vecnorm(radial,2,2);mu=radial*upstream';
    mu=min(1,max(-1,mu));mu(abs(mu-1)<=8*eps)=1;theta=acos(mu);
    onSphere=abs(vecnorm(V-origin,2,2)/R-1)<=.03;
    within=(xTip-V(:,1)<=estimate.outerDepthRadiusFraction*R+tol)&V(:,1)>=origin(1)-tol;
    geom=all(onSphere(F)&within(F),2)&sum(N.*radial,2)>=cosd(20) & ...
        ~g.featureFace&~g.topologyIssue&g.isSidewall;
    e.geometryFace=geom;e.reasonCode(:)=2;e.reasonCode(geom)=3;
    % 球面小圆帽是测地凸集：端点及驻点的X/R均在正帽内，短大圆弧不会离帽。
    % 这是拟合球面弧线，不是已经解出的三维流管或离散面元前驱树。
    cap=max(0,1-estimate.outerDepthRadiusFraction);
    complete=geom&upstream(1)>=cap-1e-10&radial(:,1)>=cap-1e-10&theta<=deg2rad(limit);
    e.completeFittedArc=complete;e.polarAngleRad(geom)=theta(geom);
    e.arcLengthM(complete)=R*theta(complete);e.transverseMetricM(complete)=R*sin(theta(complete));
    analytic=flow+mu.*radial;lengths=vecnorm(analytic,2,2);moving=complete&lengths>1e-12;
    e.analyticTangent(moving,:)=analytic(moving,:)./lengths(moving);
    e.reasonCode(complete)=4;
    if any(complete)
        shock=solvePhaseAwareWave('normal',fs.pressure,fs.temperature,speed);e.normalShock=shock;
        if shock.available
            try
                h0=air.enthalpyJkg+.5*speed^2;
                ref=phaseAwareAirState('Ph',shock.pressurePa,h0,'linearTable');T0=ref.temperatureK;
                st=phaseAwareAirCaloric(T0);post=phaseAwareAirCaloric(shock.temperatureK);
                p0=shock.pressurePa*exp((st.s0-post.s0)/air.gasConstantJkgK);
                total=phaseAwareAirState('PT',p0,T0,'linearTable');
                % dp/ds在驻点为零，速度由压力二阶曲率给有限梯度，不加任意epsilon。
                grad=sqrt(2*(p0-fs.pressure)/(total.densityKgm3*R^2));
                e.stagnation=struct('pressurePa',p0,'temperatureK',T0,'enthalpyJkg',h0, ...
                    'entropyJkgK',total.entropyJkgK,'velocityGradientPerSecond',grad, ...
                    'pressureCurvaturePam2',-2*(p0-fs.pressure)/R^2,'speedMps',0, ...
                    'gradientSource','modified-Newtonian spherical pressure curvature and tangential momentum', ...
                    'isMeshSeed',false,'boundaryLayerInitialProfileProvided',false);
                ids=find(complete);p=fs.pressure+(p0-fs.pressure)*mu(ids).^2;
                wanted=st.s0+air.gasConstantJkgK*log(p/p0);lowest=phaseAwareAirCaloric(60);
                good=wanted>=lowest.s0&wanted<=st.s0;ids=ids(good);p=p(good);wanted=wanted(good);
                if isempty(ids),error('RAH:PhaseAir:Energy','全部局部静压对应熵超出量热表。');end
                lo=60*ones(size(ids));hi=T0*ones(size(ids));
                for it=1:60
                    mid=(lo+hi)/2;c=phaseAwareAirCaloric(mid);below=c.s0<wanted;
                    lo(below)=mid(below);hi(~below)=mid(~below);
                end
                temp=(lo+hi)/2;stag=p==p0;temp(stag)=T0;
                states=phaseAwareAirState('PT',p,temp,'linearTable');dh=h0-states.enthalpyJkg;dh(stag)=0;
                assert(all(dh>=0),'RAH:PhaseNose:Energy','局部焓超过总焓，不裁负动能。');vel=sqrt(2*dh);
                recovery=phaseAwareAirState('Ph',p,states.enthalpyJkg+sqrt(states.Prandtl).*vel.^2/2,'linearTable');
                tangent=flow-sum(N(ids,:).*flow,2).*N(ids,:);tl=vecnorm(tangent,2,2);
                valid=(tl>1e-10|vel==0)&max(states.meanFreePathViscosityDefinitionM,air.meanFreePathViscosityDefinitionM)/R<.01;
                e.reasonCode(ids(tl<=1e-10&vel>0))=5;
                ent=states.entropyJkgK-total.entropyJkgK;energy=(states.enthalpyJkg+.5*vel.^2-h0)/max(1,abs(h0));
                assert(max(abs(ent))<1e-7&&max(abs(energy))<1e-9,'RAH:PhaseNose:Residual','同熵/总焓审核未过。');
                vector=zeros(numel(ids),3);mobile=tl>1e-10;vector(mobile,:)=vel(mobile).*tangent(mobile,:)./tl(mobile);
                keep=ids(valid);e.validFace(keep)=true;e.reasonCode(keep)=0;
                e.pressurePa(keep)=p(valid);e.temperatureK(keep)=temp(valid);e.speedMps(keep)=vel(valid);
                e.velocityMps(keep,:)=vector(valid,:);e.Mach(keep)=vel(valid)./states.soundSpeedMps(valid);
                e.densityKgm3(keep)=states.densityKgm3(valid);e.recoveryTemperatureK(keep)=recovery.temperatureK(valid);
                e.entropyResidual(keep)=ent(valid);e.energyResidual(keep)=energy(valid);
            catch err
                if ~ismember(err.identifier,{'RAH:PhaseAir:Domain','RAH:PhaseAir:Temperature','RAH:PhaseAir:Energy'}),rethrow(err);end
                % 本区域共用一个总态，域外时整区拒绝，不换物性或裁温。
                e.domainRejection=err.identifier;
            end
        end
    end
end

%% 来源绑定：实际原面、工况、算法代码与输出数组共同生成指纹
sources={'solvePhaseAwareNoseEdge.m','estimateNoseRadius.m','solvePhaseAwareWave.m', ...
    'phaseAwareAirState.m','phaseAwareAirCaloric.m','phaseAwarePayloadSHA.m'};hashes=struct;
for k=1:numel(sources)
    fid=fopen(fullfile(fileparts(mfilename('fullpath')),sources{k}),'rb');clean=onCleanup(@()fclose(fid));
    bytes=fread(fid,Inf,'*uint8');md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
    hashes.(matlab.lang.makeValidName(sources{k}))=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));clear clean;
end
e.sourceText=jsonencode(struct('model',e.model,'transportModel',e.transportModel,'energyReference',e.energyReference, ...
    'meshSHA256',e.meshSHA256,'condition',e.condition,'method',e.method,'scope',e.scope,'codeSHA256',hashes, ...
    'propertyTableSHA256',air.sourceSHA256,'sourceURL',e.sourceURL,'sphereCenter',e.sphereCenter,'stagnation',e.stagnation));
digest=phaseAwarePayloadSHA(e.sourceText,e.pressurePa,e.temperatureK,e.speedMps,e.velocityMps, ...
    e.recoveryTemperatureK,e.validFace,e.reasonCode,e.arcLengthM,e.Mach,e.faceAreaM2);
e.source=struct('id','internal-phase-aware-spherical-nose-v1','role','engineering-candidate', ...
    'sha256',digest,'method',e.method);
e.summary=struct('totalAreaM2',sum(area),'candidateAreaM2',sum(area(e.validFace)), ...
    'candidateFaces',nnz(e.validFace),'candidateAreaFraction',sum(area(e.validFace))/sum(area),'heatAreaM2',0);
end
