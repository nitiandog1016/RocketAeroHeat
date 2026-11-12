function report=testPhaseAwareNoseV226(folder)
%TESTPHASEAWARENOSEV226 真正原面生产/主工况/主工作台出口与相态球鼻控制。
d=load('outputs/v104_validation/baseline02/current_case.mat','current','cfg','fs');
original=d.current;cfg=d.cfg;cfg.analysisMode='phaseAwareNoseEdge';
[m,fs,caseReport]=runAeroHeatSingleCase(original,cfg);e=m.phaseAwareNoseEdge;
assert(e.summary.candidateFaces>0&&isequal(m.V,original.V)&&isequal(m.F,original.F));
assert(size(m.F,1)==19752&&abs(e.summary.totalAreaM2-1.9646292565715553)<2e-13);
assert(~isfield(m,'heating')&&~isfield(m,'pressure')&&~isfield(m,'externalFlow'));
assert(~e.heatReady&&~e.profileReady&&caseReport.heatComputedAreaFraction==0);
assert(abs(fs.speed-d.fs.speed)>.1&&strcmp(fs.model,e.model));
v=e.validFace;assert(all(e.completeFittedArc(v))&&all(e.geometryFace(v)));
assert(all(isnan(e.pressurePa(~v)))&&max(abs(sum(e.velocityMps(v,:).*original.normal(v,:),2)))<1e-8);
assert(max(abs(e.entropyResidual(v)))<1e-7&&max(abs(e.energyResidual(v)))<1e-9);
assert(all(e.recoveryTemperatureK(v)>=e.temperatureK(v)));
rho0=phaseAwareAirState('PT',e.stagnation.pressurePa,e.stagnation.temperatureK,'linearTable');
assert(abs(e.stagnation.pressureCurvaturePam2+rho0.densityKgm3*e.stagnation.velocityGradientPerSecond^2)<1e-6);
payload=buildAerothermalViewerData(m);assert(numel(payload.fields)==9&&strcmp(payload.defaultField,'pressure'));
assert(isequaln(payload.fields(1).values,e.pressurePa));
assert(all(isnan(payload.fields(7).values))&&all(payload.fields(8).status==2)&&all(isnan(payload.fields(9).values)));
write('rocket_workbench.json',payload);write('rocket_metadata.json',rmfield(e, ...
    {'pressurePa','temperatureK','speedMps','velocityMps','Mach','recoveryTemperatureK','densityKgm3', ...
    'entropyResidual','energyResidual','polarAngleRad','arcLengthM','transverseMetricM','analyticTangent', ...
    'validFace','geometryFace','completeFittedArc','reasonCode','faceAreaM2','originalFaceID'}));
rows=[(1:numel(v))' e.faceAreaM2 e.pressurePa e.temperatureK e.speedMps e.densityKgm3 ...
    e.velocityMps e.Mach e.recoveryTemperatureK e.polarAngleRad e.arcLengthM ...
    double(v) double(e.reasonCode)];
fid=fopen(fullfile(folder,'rocket_faces.csv'),'w');clean=onCleanup(@()fclose(fid));
fprintf(fid,'faceID,area,p,Te,U,rho,ux,uy,uz,Mach,Taw,theta,arc,valid,reason\n');
fprintf(fid,[repmat('%.17g,',1,14) '%.17g\n'],rows');clear clean;
% 显示入口必须拒绝原面或已绑定字段被换成旧结果。
for k=1:3
    bad=m;id=find(v,1);
    if k==1,bad.V(1,1)=bad.V(1,1)+1e-4;
    elseif k==2,bad.phaseAwareNoseEdge.pressurePa(id)=e.pressurePa(id)+1;
    else,bad.phaseAwareNoseEdge.Mach(id)=e.Mach(id)+1;end
    rejected=false;try,buildAerothermalViewerData(bad);catch err,rejected=startsWith(err.identifier,'RAH:PhaseNoseViewer:');end
    assert(rejected);
end
% 解析球鼻只查新恢复/零速度与几何域，复用成熟几何模块，不扫分辨率。
[control,~]=buildCanonicalBody('hemisphereCylinder',struct('resolution',32));
evalc('control=calcFaceGeometry(control);control=fixNormalDirection(control,cfg);control=buildSurfaceCoordinates(control,cfg);');
control=analyzeLocalGeometry(control,cfg);
ce=solvePhaseAwareNoseEdge(control,fs,cfg);assert(any(ce.validFace));
radial=control.center-ce.sphereCenter;radial=radial./vecnorm(radial,2,2);
id=find(ce.validFace&radial(:,1)>cosd(20),1);aligned=fs;aligned.unitFlowDirection=-radial(id,:);
zero=solvePhaseAwareNoseEdge(control,aligned,cfg);assert(zero.validFace(id)&&zero.speedMps(id)==0);
assert(zero.pressurePa(id)==zero.stagnation.pressurePa&&all(zero.velocityMps(id,:)==0));
assert(isfinite(zero.stagnation.velocityGradientPerSecond)&&zero.stagnation.velocityGradientPerSecond>0);
marked=control;marked.localGeometry.featureFace(id)=true;
excluded=solvePhaseAwareNoseEdge(marked,fs,cfg);assert(~excluded.validFace(id));
outside=fs;outside.unitFlowDirection=[-cosd(35) 0 sind(35)];
refused=solvePhaseAwareNoseEdge(control,outside,cfg);assert(~any(refused.validFace));
rare=fs;rare.pressure=1e-5;refused=solvePhaseAwareNoseEdge(control,rare,cfg);assert(~any(refused.validFace));
hot=fs;hot.Mach=12;refused=solvePhaseAwareNoseEdge(control,hot,cfg);assert(~any(refused.validFace));
limited=cfg;limited.externalFlow.sphericalNoseMaxPolarAngleDeg=25;
narrow=solvePhaseAwareNoseEdge(control,fs,limited);assert(nnz(narrow.validFace)<nnz(ce.validFace)&&~any(narrow.validFace&~ce.validFace));
bad=original;bad.V(1,1)=bad.V(1,1)+.01;rejected=false;
try,solvePhaseAwareNoseEdge(bad,fs,cfg);catch err,rejected=strcmp(err.identifier,'RAH:PhaseNose:Geometry');end
assert(rejected);
guard=cfg;guard.heating.sidewallModel='phaseAwareLaminarIntegral';rejected=false;
try,runAeroHeatSingleCase(original,guard);catch err,rejected=strcmp(err.identifier,'RAH:PhaseLayer:ProviderNotConnected');end
assert(rejected);
report=struct('passed',true,'realRocket',e.summary,'mainSingleCaseUsed',true,'mainViewerUsed',true, ...
    'fields',9,'controlFaces',size(control.F,1),'controlCandidateFaces',nnz(ce.validFace), ...
    'exactStagnationZero',true,'finiteStagnationGradient',true,'tamperRejections',3,'scopeRejections',6, ...
    'layerProviderGuardPreserved',true,'physicalQualification',false,'newQualifiedAreaM2',0, ...
    'heatComputedAreaM2',0,'geometrySHA256',e.meshSHA256,'sourceSHA256',e.source.sha256);
write('report.json',report);disp(report);
    function write(name,value)
        fid=fopen(fullfile(folder,name),'w','n','UTF-8');clean=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(value));
    end
end
