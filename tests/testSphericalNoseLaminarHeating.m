function summary=testSphericalNoseLaminarHeating(outputFolder)
%TESTSPHERICALNOSELAMINARHEATING 球鼻局部相似热边界的独立数值校核。
% 模块作用：独立复算正激波/球面压力、曲率及xi，核对同剖面壁面
% 导数、绝热根、角度表收敛和拒绝状态。不是实验或区域接口验证。
if nargin<1, outputFolder='outputs/v54_validation'; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
cfg.flight.Mach=4; cfg.flight.altitude=20000;
cfg.flight.alpha=0; cfg.flight.beta=0; cfg.flight.Twall=300;
cfg.heating.sphericalNoseAngleStepDeg=5;
fs=calcFreestream(cfg.flight);
[mesh,meta]=buildCanonicalBody('hemisphereCylinder',struct('resolution',24));
mesh=prepare(mesh,cfg); before=mesh;

%% 1. NASA7实际积分：全部可用球鼻热根/暖壁支/端点及逐场量来源
tic; nasa=solveSphericalNoseLaminarHeating(mesh,fs,300,cfg); nasaSeconds=toc;
fprintf('NASA7 observed sphere: %d/%d eligible faces, %d stations, %d BVP calls, %.3f s.\n', ...
    nnz(nasa.faceMask),nnz(nasa.edgeState.algebraicallyValid),nasa.uniqueStations,nasa.profileCalls,nasaSeconds);
assert(nasa.available,nasa.status);
assert(isequaln(mesh,before),'区域解不能改变网格/原几何结果');
assert(~nasa.historyCompatible && ~nasa.nonSimilarTermsSolved && ...
    ~nasa.entropySwallowingSolved && ~nasa.interfaceSolved && ~nasa.experimentalValidation);
assert(all(nasa.heatFlux(nasa.faceMask)>0));
assert(all(isnan(nasa.heatFlux(~nasa.faceMask))));
assert(~any(nasa.faceMask(strcmp(meta.faceRegion,'cylinder') | strcmp(meta.faceRegion,'base'))));
assert(strcmp(nasa.fieldMetadata.wallTemperatureInput.valueKind,'prescribed'));
highAngleStations=0; maximumHeatRootResidual=0;
for k=1:numel(nasa.stationResults)
    node=nasa.stationResults{k};
    assert(isfield(node,'adiabaticProfile'),'当前受支持球鼻站仍有未解热边界');
    checkProfile(node,nasa,fs);
    checkSources(node,nasa);
    t=node.adiabaticBracketTrials;
    maximumHeatRootResidual=max(maximumHeatRootResidual,abs(node.adiabaticProfile.heatFlux));
    if node.polarAngleDeg>0
        % 只接受实际有效且异号的端点；不强制制造历史错误门限的拒绝。
        ends=node.adiabaticBracket;
        endpointRows=zeros(1,2);
        for e=1:2
            row=find(abs(t(:,1)-ends(e))<1e-8,1,'last');
            assert(~isempty(row) && t(row,2)==1);
            endpointRows(e)=row;
        end
        assert(t(endpointRows(1),4)>=0 && t(endpointRows(2),4)<=0);
        assert(node.recoveryTemperature>=ends(1)-1e-7 && node.recoveryTemperature<=ends(2)+1e-7);
    end
    if node.polarAngleDeg>=45
        highAngleStations=highAngleStations+1;
        ad=node.adiabaticProfile;
        assert(ad.wallShearStress>0 && all(ad.profile.speed>=-1e-8));
        assert(ad.velocityAudit.maximumRatio>1+1e-6 || ...
            ad.velocityAudit.minimumEtaDerivative<-1e-6, ...
            '必须覆盖曾被全域单调门限拒绝的实际高角度暖壁支');
    end
end
assert(highAngleStations>=3 && max(nasa.angularProfileTable.polarAngleDeg)>54);
assert(all(nasa.angularProfileTable.statusCode==0) && ...
    all(nasa.faceMask(nasa.edgeState.algebraicallyValid)), ...
    '该冷壁M4/20km球鼻有效范围中仍存在未解节点/面元');
highAngleFaces=nasa.faceMask & nasa.polarAngleDeg>=45;
assert(nnz(highAngleFaces)>0);
checkFaceInterpolation(nasa,300);
fprintf('High-angle thermal roots: %d stations / %d faces >=45 deg; max |qad| %.9g W/m2.\n', ...
    highAngleStations,nnz(highAngleFaces),maximumHeatRootResidual);
result=nasa; save(fullfile(outputFolder,'sphere_laminar_case.mat'),'mesh','cfg','fs','meta','result');

%% 2. 定比热独立解析与均匀弧长积分；不调用被测模块的恢复/xi函数
caloric=cfg; caloric.heating.airThermalModel='caloricallyPerfect';
caloric.externalFlow.sphericalNoseMaxPolarAngleDeg=25;
c=solveSphericalNoseLaminarHeating(mesh,fs,300,caloric); assert(c.available,c.status);
g=fs.gamma; cp=g*fs.R/(g-1); M=fs.Mach; Rn=.1;
p2=fs.pressure*(2*g*M^2-(g-1))/(g+1);
M2sq=((g-1)*M^2+2)/(2*g*M^2-(g-1));
p0=p2*(1+(g-1)*M2sq/2)^(g/(g-1));
T0=fs.temperature+fs.speed^2/(2*cp); rho0=p0/(fs.R*T0);
a=sqrt(2*(p0-fs.pressure)/(rho0*Rn^2));
assert(abs(c.noseRadius/Rn-1)<1e-10);
assert(abs(c.stagnationVelocityGradient/a-1)<1e-10);
assert(c.stagnationMomentumCurvatureResidual<1e-12);
assert(abs(c.stagnationPressureCurvature/(-rho0*a^2)-1)<1e-10);
xiError=0; gradientError=0;
for k=1:numel(c.stationResults)
    node=c.stationResults{k}; phi=deg2rad(node.polarAngleDeg);
    e=node.edge; s=node.station;
    p=fs.pressure+(p0-fs.pressure)*cos(phi)^2;
    T=T0*(p/p0)^((g-1)/g); rho=p/(fs.R*T);
    U=sqrt(2*cp*(T0-T));
    assert(abs(e.pressure/p-1)<1e-11 && abs(e.temperature/T-1)<1e-11);
    assert(abs(e.density/rho-1)<1e-11);
    if phi==0
        assert(e.speed==0 && s.xi==0 && s.beta==.5);
        assert(abs(s.velocityGradient/a-1)<1e-10);
    else
        assert(abs(e.speed/U-1)<1e-10);
        t=linspace(0,phi,20001); pp=fs.pressure+(p0-fs.pressure)*cos(t).^2;
        tt=T0*(pp/p0).^((g-1)/g); rr=pp./(fs.R*tt);
        uu=sqrt(max(0,2*cp*(T0-tt)));
        independentXi=trapz(t,Rn*rr.*sutherland(tt).*uu.*(Rn*sin(t)).^2);
        xiError=max(xiError,abs(s.xi/independentXi-1));
        assert(abs(s.xi/independentXi-1)<5e-8);
        dP=-(p0-fs.pressure)*sin(2*phi)/Rn;
        gradientError=max(gradientError,abs((s.pressureGradient+rho*U*s.velocityGradient)/dP));
        assert(gradientError<1e-10);
        independentBeta=2*independentXi*(-dP/(rho*U))/(rho*sutherland(T)*U^2*(Rn*sin(phi))^2);
        assert(abs(s.beta/independentBeta-1)<5e-8);
    end
    if isfield(node,'adiabaticProfile'), checkProfile(node,c,fs); end
end
fprintf('Independent perfect-gas xi/curvature/momentum: max xi relative error %.3g, gradient %.3g.\n',xiError,gradientError);

%% 3. 角度表收敛、逐面元插值来源以及热/冷壁的实际符号
fineCfg=caloric; fineCfg.heating.sphericalNoseAngleStepDeg=2.5;
fine=solveSphericalNoseLaminarHeating(mesh,fs,300,fineCfg);
common=c.faceMask & fine.faceMask; assert(nnz(common)>=10);
qDifference=max(abs(c.heatFlux(common)./fine.heatFlux(common)-1));
tDifference=max(abs(c.recoveryTemperature(common)-fine.recoveryTemperature(common)));
assert(qDifference<.003 && tDifference<.15,'5/2.5度热边界表未达到本算例收敛检查');
ids=find(c.faceMask);
qcol=find(strcmp(c.angularProfileTable.fieldNames,'heatFlux'));
for id=ids(:)'
    pair=c.profileNodeIDs(id,:); weight=c.profileNodeWeight(id);
    assert(all(c.angularProfileTable.statusCode(pair)==0));
    expected=(1-weight)*c.angularProfileTable.values(pair(1),qcol)+ ...
        weight*c.angularProfileTable.values(pair(2),qcol);
    assert(abs(expected-c.heatFlux(id))<1e-10*max(abs(expected),1));
    assert(abs(c.heatTransferCoefficient(id)*(c.recoveryTemperature(id)-300)-c.heatFlux(id)) ...
        <1e-10*max(abs(c.heatFlux(id)),1));
end
warmCfg=caloric; warmCfg.externalFlow.sphericalNoseMaxPolarAngleDeg=25;
warm=solveSphericalNoseLaminarHeating(mesh,fs,T0-2,warmCfg);
assert(warm.available && any(warm.heatFlux(warm.faceMask)<0) && ...
    all(c.heatFlux(c.faceMask)>0), ...
    '冷壁和高于当地恢复温度的暖壁须分别保留真实加热/冷却符号');
assert(all(warm.heatTransferCoefficient(warm.faceMask)>0));
fprintf('Angle5/2.5 convergence: q %.3g relative, Taw %.3g K; heating/cooling checked.\n',qDifference,tDifference);

%% 4. 入射角/几何/壁温/Kn/Re/薄层/有限边缘拒绝，缺失不能补零
oblique=fs; oblique.unitFlowDirection=[-cosd(1) sind(1) 0];
bad=solveSphericalNoseLaminarHeating(mesh,oblique,300,cfg);
assert(~bad.available && any(bad.statusCode==7) && bad.profileCalls==0); checkBlank(bad);
wrong=mesh; wrong.V(:,2)=1.2*wrong.V(:,2); wrong=prepare(wrong,cfg);
bad=solveSphericalNoseLaminarHeating(wrong,fs,300,cfg);
assert(~bad.available && bad.profileCalls==0); checkBlank(bad);
bad=solveSphericalNoseLaminarHeating(mesh,fs,NaN,cfg);
assert(~bad.available && bad.profileCalls==0); checkBlank(bad);
bad=solveSphericalNoseLaminarHeating(mesh,fs,70,cfg);
assert(~bad.available && any(bad.statusCode==3)); checkBlank(bad);
rare=withPressure(fs,.01);
bad=solveSphericalNoseLaminarHeating(mesh,rare,300,cfg);
assert(~bad.available && bad.freestreamKnudsen>=.01 && any(bad.statusCode==2)); checkBlank(bad);
shortCfg=caloric; shortCfg.externalFlow.sphericalNoseMaxPolarAngleDeg=10;
shortCfg.massProfile=struct('etaMax',4,'edgeBuffer',1);
bad=solveSphericalNoseLaminarHeating(mesh,fs,300,shortCfg);
assert(~bad.available && any(bad.statusCode==5)); checkBlank(bad);
reCfg=caloric; reCfg.externalFlow.sphericalNoseMaxPolarAngleDeg=10;
reCfg.heating.transitionRe=1;
bad=solveSphericalNoseLaminarHeating(mesh,fs,300,reCfg);
assert(~bad.available && any(bad.statusCode==4)); checkBlank(bad);
thinCfg=caloric; thinCfg.externalFlow.sphericalNoseMaxPolarAngleDeg=10;
lowPressure=withPressure(fs,5);
bad=solveSphericalNoseLaminarHeating(mesh,lowPressure,300,thinCfg);
assert(~bad.available && bad.freestreamKnudsen<.01 && any(bad.statusCode==6), ...
    '应在连续介质筛选内显式拒绝delta/R过大的剖面'); checkBlank(bad);
fprintf('Incidence/geometry/unknown and out-of-domain wall/Kn/Re/delta-R/finite-edge rejection passed.\n');
summary=struct('tested',true,'experimentalValidation',false,'interfaceSolved',false, ...
    'nasaEligibleFaces',nnz(nasa.edgeState.algebraicallyValid),'nasaComputedFaces',nnz(nasa.faceMask), ...
    'nasaSeconds',nasaSeconds,'xiRelativeError',xiError,'momentumRelativeError',gradientError, ...
    'nasaProfileCalls',nasa.profileCalls,'nasaUniqueStations',nasa.uniqueStations, ...
    'highAngleComputedFaces',nnz(highAngleFaces),'highAngleStations',highAngleStations, ...
    'maximumAdiabaticHeatResidualWm2',maximumHeatRootResidual, ...
    'eligibleAreaM2',sum(mesh.area(nasa.edgeState.algebraicallyValid)), ...
    'computedAreaM2',sum(mesh.area(nasa.faceMask)), ...
    'totalMeshAreaM2',sum(mesh.area),'computedTotalAreaFraction',sum(mesh.area(nasa.faceMask))/sum(mesh.area), ...
    'angleHeatRelativeDifference',qDifference,'angleRecoveryDifferenceK',tDifference, ...
    'incidenceScope','zero combined incidence only');
save(fullfile(outputFolder,'sphere_laminar_results.mat'),'summary','cfg','fs','mesh','meta','nasa','c','fine','warm');
fid=fopen(fullfile(outputFolder,'sphere_laminar_summary.json'),'w');
guard=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(summary,PrettyPrint=true)); clear guard;
fprintf('Spherical nose local-similarity analytic/numerical tests passed; no experiment/interface claim.\n');
end
function checkSources(node,region)
% 模块作用：逐场量检查真实算法来源及给定/输入/求解状态。
for key={'heatFlux','recoveryTemperature','heatTransferCoefficient'}
    m=region.fieldMetadata.(key{1});
    assert(strcmp(m.valueKind,'computed') && ~isempty(m.method) && ~isempty(m.scope));
end
assert(strcmp(region.fieldMetadata.wallTemperatureInput.valueKind,'prescribed'));
assert(contains(region.heatFluxMethod,'壁面导数') && ...
    contains(region.recoveryTemperatureMethod,'q(Twall)=0') && ...
    contains(region.heatTransferCoefficientMethod,'q/(Taw-Twall)'));
assert(numel(region.sourceURLs)>=2 && contains(region.edgeState.pressureMethod,'修正牛顿') && ...
    contains(region.edgeState.stateMethod,'恢复'));
for key={'pressure','temperature','density','speed','totalEnthalpy'}
    assert(~isempty(node.edge.sources.(key{1})));
    assert(strcmp(node.wallProfile.fieldMetadata.(key{1}).valueKind,'input'));
end
for key={'streamlineMetric','xi','beta','pressureGradient','velocityGradient','wallTemperature'}
    assert(~isempty(node.station.sources.(key{1})));
end
assert(strcmp(node.wallProfile.fieldMetadata.heatFlux.valueKind,'computed') && ...
    node.wallProfile.fieldMetadata.heatFlux.valid && ...
    strcmp(node.wallProfile.fieldMetadata.wallTemperatureInput.valueKind,'prescribed'));
end
function checkFaceInterpolation(region,Twall)
% 模块作用：所有已解真面元对应有效BVP节点；不跨未解值/变更来源。
for id=find(region.faceMask)'
    pair=region.profileNodeIDs(id,:); w=region.profileNodeWeight(id);
    assert(all(isfinite(pair)) && all(region.angularProfileTable.statusCode(pair)==0));
    for key={'heatFlux','recoveryTemperature','wallShearStress','massThickness','delta','xi','beta'}
        col=find(strcmp(region.angularProfileTable.fieldNames,key{1}));
        value=(1-w)*region.angularProfileTable.values(pair(1),col)+ ...
            w*region.angularProfileTable.values(pair(2),col);
        assert(abs(value-region.(key{1})(id))<=1e-10*max(abs(value),1));
    end
    assert(region.statusCode(id)==0 && region.wallTemperatureInput(id)==Twall);
    assert(region.heatTransferCoefficient(id)>0 && ...
        abs(region.heatTransferCoefficient(id)*(region.recoveryTemperature(id)-Twall)- ...
        region.heatFlux(id))<=1e-10*max(abs(region.heatFlux(id)),1));
end
end
function checkProfile(node,region,fs)
% 物理壁坐标导数用相邻剖面重建，不调用被测模块的热流公式。
wp=node.wallProfile; ad=node.adiabaticProfile; assert(wp.available && ad.available);
assert(abs(ad.heatFlux)<=1e-5*max(abs(wp.heatFlux),1));
assert(wp.delta/region.noseRadius<=region.maximumDeltaRadiusFraction);
assert(wp.inputMomentumRelativeResidual<1e-10);
assert(node.recoveryTemperature>=node.edge.temperature-1e-7 && ...
    node.recoveryTemperature<=region.edgeState.normalShock.stagnationTemperature+1e-7);
Tw=wp.wallTemperatureInput; y=wp.profile.physicalY;
if strcmp(wp.thermalPropertyMode,'thermallyPerfect')
    [~,enthalpy]=airThermalProperties(wp.profile.temperature);
else
    enthalpy=fs.gamma*fs.R/(fs.gamma-1)*wp.profile.temperature;
end
% 五个近壁点的四阶多项式给独立物理dh/dy；单位k=mu*cp/Pr。
scale=y(5); polynomial=polyfit(y(1:5)/scale,enthalpy(1:5)-enthalpy(1),4);
physicalDerivative=polynomial(end-1)/scale;
expectedQ=sutherland(Tw)/wp.numericalOptions.Prandtl*physicalDerivative;
assert(abs(expectedQ/wp.heatFlux-1)<2e-4,'热流与同剖面物理焓壁导数不一致');
assert(ad.edgeChecks.heatGradientCriterionDefect<=1+1e-8);
end
function checkBlank(r)
assert(~any(r.faceMask) && all(isnan(r.heatFlux)) && all(isnan(r.recoveryTemperature)));
assert(all(isnan(r.heatTransferCoefficient)));
end
function mesh=prepare(mesh,cfg)
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
end
function value=sutherland(T)
value=1.458e-6*T.^(3/2)./(T+110.4);
end
function value=withPressure(fs,pressure)
% 同时变换理想气体密度与总压，避免不一致输入触发上游结构校验。
value=fs; ratio=pressure/fs.pressure;
value.pressure=pressure; value.density=fs.density*ratio;
value.totalPressure=fs.totalPressure*ratio;
value.dynamicPressure=fs.dynamicPressure*ratio; value.qinf=value.dynamicPressure;
end
