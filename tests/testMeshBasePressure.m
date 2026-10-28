function summary = testMeshBasePressure(outputDir)
%TESTMESHBASEPRESSURE 真实网格底压的区域解、面积、来源和拒绝校核。
% 模块用途：用独立解析圆锥区域、报告公式与物理剖面积分，核对
% 真实底面的计算值；检验平移、相似尺度、几何及肩部审核失败。
% 这些属于软件/解析/关联式校核，不是新的独立实验精度验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,outputDir=fullfile(root,'outputs','v54_validation','base_mesh');end
if ~exist(outputDir,'dir'),mkdir(outputDir);end
[m,meta]=buildCanonicalBody('cone',struct('radius',.2,'coneAngleDeg',9,'resolution',24));
f=makeState(180,250,3);
c=struct('heating',struct('airThermalModel','caloricallyPerfect'), ...
    'flight',struct('Twall',300), ...
    'basePressure',struct('hasRearSupport',false,'hasGasInjection',false, ...
    'operationSource','隔离解析尖圆锥算例：无尾撑、无底部喷流'));

%% 1. 实际底面值和独立报告公式/面面积；热流始终未映射
r=solveMeshBasePressure(m,f,c);
assert(r.pressureAvailable,r.status);
base=strcmp(meta.faceRegion,'base');
assert(isequal(r.baseFace,base) && isequal(r.computedFace,base));
assert(r.geometry.flatBaseVerified && r.geometry.verified);
independentArea=vecnorm(cross(m.V(m.F(:,2),:)-m.V(m.F(:,1),:), ...
    m.V(m.F(:,3),:)-m.V(m.F(:,1),:),2),2,2)/2;
e=r.shoulder;Re=e.density*e.speed*e.distance/e.dynamicViscosity;
C4=(1+5/e.Mach^2)^(-2);
pb=f.pressure*3.05*(e.Mach^2/sqrt(Re))^.6/C4;
cp=(pb-f.pressure)/(.5*f.density*f.speed^2);
assert(max(abs(r.pressure(base)/pb-1))<1e-12);
assert(max(abs(r.Cp(base)-cp))<1e-12 && cp<0);
assert(all(isnan(r.pressure(~base))) && all(isnan(r.Cp(~base))));
assert(abs(r.summary.computedPressureArea/sum(independentArea(base))-1)<1e-12);
assert(abs(r.summary.pressureAreaFraction-sum(independentArea(base))/sum(independentArea))<1e-12);
assert(r.summary.basePressureAreaFraction==1 && r.summary.computedHeatArea==0);
assert(r.pressureEligibleForSurfaceCoverage && ~r.heatEligibleForSurfaceCoverage);
assert(all(isnan(r.heatFlux)) && all(isnan(r.edgeTemperature)) && ...
    all(isnan(r.recoveryTemperature)) && all(isnan(r.solvedWallTemperature)));
assert(~r.radialHeatingDistributionSolved && ~r.regionRepresentative.eligibleForSurfaceCoverage);
assert(strcmp(r.fieldMetadata.pressure.valueKind,'computed'));
assert(strcmp(r.fieldMetadata.heatFlux.valueKind,'unresolved'));
assert(strcmp(r.fieldMetadata.wallTemperatureInput.valueKind,'prescribed'));
assert(strcmp(r.fieldMetadata.boundaryLayerRegime.valueKind,'computed') && ...
    r.fieldMetadata.boundaryLayerRegime.criticalReynolds==5e5);
assert(strcmp(r.operation.valueKind,'input') && ~isempty(r.operation.source));
assert(all(r.faceStatus(base)==0) && all(r.faceStatus(~base)==1));
assert(all(contains(r.faceMethod(base),'Eq.(8b)')));
assert(r.shoulderTransition.laminarFace && Re>=1e5 && Re<5e5);

%% 2. 同一物理剖面的厚度证据，不能以theta任意换算delta
b=r.shoulderBoundaryLayer;
assert(b.available && b.delta>0 && r.shoulderThicknessRatio<=.05);
assert(b.inputMomentumRelativeResidual<1e-12 && b.stationInput.beta==0);
assert(abs(b.stationInput.xi/(e.density*e.dynamicViscosity*e.speed*.2^2*e.distance/3)-1)<1e-12);
assert(abs(b.delta-b.profile.physicalY(end))/b.delta<2e-5);
mass=trapz(b.profile.physicalY,b.profile.density/e.density.*b.profile.speed/e.speed);
assert(abs(mass/b.massThickness-1)<2e-5);
assert(~b.nonSimilarTermsSolved && ~b.massMatchingSolved);

%% 3. 平移不变；维持Re和Kn的相似放大，压力尺度按1/k变化
translated=m;translated.V=m.V+[5 -3 7];
t=solveMeshBasePressure(translated,f,c);assert(t.available,t.status);
assert(max(abs(t.pressure(base)./r.pressure(base)-1))<1e-10);
assert(abs(t.geometry.generatorLength/r.geometry.generatorLength-1)<1e-12);
assert(abs(t.shoulderThicknessRatio/r.shoulderThicknessRatio-1)<1e-9);
factor=2; scaled=m;scaled.V=m.V*factor;
scaledF=makeState(f.temperature,f.pressure/factor,f.Mach);
ss=solveMeshBasePressure(scaled,scaledF,c);assert(ss.available,ss.status);
assert(max(abs(ss.pressure(base)./r.pressure(base)*factor-1))<1e-10);
assert(max(abs(ss.Cp(base)-r.Cp(base)))<1e-12);
assert(abs(ss.shoulderThicknessRatio/r.shoulderThicknessRatio-1)<1e-9);
assert(abs(ss.summary.computedPressureArea/r.summary.computedPressureArea/factor^2-1)<1e-12);
assert(abs(ss.correlation.continuumKnudsen/r.correlation.continuumKnudsen-1)<1e-10);

%% 4. 真实几何拒绝：名称/声明不能覆盖非圆底、钝鼻、开口或多组件
bad=m;bad.V(:,3)=bad.V(:,3)*1.05;
nonCircular=solveMeshBasePressure(bad,f,c);assertRejected(nonCircular,2);
assert(~nonCircular.geometry.flatBaseVerified);
bad=m;bad.V(:,1)=bad.V(:,1)+.01*bad.V(:,2);
tilted=solveMeshBasePressure(bad,f,c);assertRejected(tilted,2);
assert(~tilted.geometry.flatBaseVerified && ~any(tilted.baseFace));
bad=m;bad.V(1,1)=bad.V(1,1)+.002;
partial=solveMeshBasePressure(bad,f,c);assertRejected(partial,2);
assert(~partial.geometry.flatBaseVerified && nnz(partial.baseFace)==nnz(base)-2);
bad=m;bad.F(find(base,1),:)=bad.F(find(base,1),[1 3 2]);
wrongNormal=solveMeshBasePressure(bad,f,c);assertRejected(wrongNormal,2);
bad=m;bad.F(find(base,1),:)=[];
openBase=solveMeshBasePressure(bad,f,c);assertRejected(openBase,2);
[rounded,~]=buildCanonicalBody('bluntedCone',struct('radius',.2, ...
    'noseRadius',.02,'cylinderLength',0,'coneAngleDeg',9,'resolution',24));
blunt=solveMeshBasePressure(rounded,f,c);assertRejected(blunt,2);
[cylinder,~]=buildCanonicalBody('coneCylinder',struct('radius',.2, ...
    'coneAngleDeg',9,'resolution',24));
coneCylinder=solveMeshBasePressure(cylinder,f,c);assertRejected(coneCylinder,2);
assert(coneCylinder.geometry.flatBaseVerified && ~coneCylinder.geometry.verified);
bad=m;nv=size(bad.V,1);bad.V=[bad.V;2 2 2;2.1 2 2;2 2.1 2;2 2 2.1];
bad.F=[bad.F;nv+[1 3 2;1 2 4;2 3 4;3 1 4]];
multiple=solveMeshBasePressure(bad,f,c);assertRejected(multiple,2);

%% 5. 操作、攻角、同物性状态、转捩及物理薄层失败必须保留未解
ca=c;ca=rmfield(ca,'basePressure');missing=solveMeshBasePressure(m,f,ca);assertRejected(missing,3);
ca=c;ca.basePressure.operationSource='';missingSource=solveMeshBasePressure(m,f,ca);assertRejected(missingSource,3);
ca=c;ca.basePressure.hasRearSupport=true;support=solveMeshBasePressure(m,f,ca);assertRejected(support,3);
ca=c;ca.basePressure.hasGasInjection=true;jet=solveMeshBasePressure(m,f,ca);assertRejected(jet,3);
fa=f;fa.alpha=2;incidence=solveMeshBasePressure(m,fa,c);assertRejected(incidence,4);
fa=f;fa.beta=1;beta=solveMeshBasePressure(m,fa,c);assertRejected(beta,4);
ca=c;ca.heating.airThermalModel='thermallyPerfect';thermallyPerfect=solveMeshBasePressure(m,f,ca);assertRejected(thermallyPerfect,4);
fa=f;fa.totalPressure=fa.totalPressure*1.05;inconsistent=solveMeshBasePressure(m,fa,c);assertRejected(inconsistent,4);
fa=makeState(180,1000,3);turbulent=solveMeshBasePressure(m,fa,c);assertRejected(turbulent,6);
fa=makeState(180,200,3);thick=solveMeshBasePressure(m,fa,c);assertRejected(thick,6);
assert(thick.shoulderBoundaryLayer.available && thick.shoulderThicknessRatio>.05);
ca=c;ca.massProfile=struct('etaMax',4);unresolvedProfile=solveMeshBasePressure(m,f,ca);assertRejected(unresolvedProfile,6);
assert(isfield(unresolvedProfile,'shoulderBoundaryLayer') && ~unresolvedProfile.shoulderBoundaryLayer.available);
ca=c;ca.heating.wallTemperatureByFace=300*ones(size(m.F,1),1);ca.heating.wallTemperatureByFace(1)=310;
nonuniform=solveMeshBasePressure(m,f,ca);assertRejected(nonuniform,6);

%% 6. 保留定量证据与真实计算面积，绝不将区域审核称为全面积热解
summary=struct('baseFaces',nnz(base),'totalFaces',size(m.F,1), ...
    'basePressurePa',pb,'baseCp',cp,'shoulderReynolds',Re, ...
    'shoulderDeltaM',b.delta,'shoulderThicknessRatio',r.shoulderThicknessRatio, ...
    'computedPressureAreaM2',r.summary.computedPressureArea, ...
    'computedPressureAreaFraction',r.summary.pressureAreaFraction, ...
    'computedHeatAreaM2',0,'independentExperimentalValidation',false, ...
    'evidenceKind','真实网格/解析几何/公开关联式及同剖面数值校核', ...
    'unsupportedModelsRemainUnresolved',true);
fid=fopen(fullfile(outputDir,'mesh_base_summary.json'),'w');
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(summary));clear cleanup;
save(fullfile(outputDir,'mesh_base_results.mat'),'m','meta','f','c','r','summary', ...
    't','ss','thick','turbulent','unresolvedProfile');
fprintf('testMeshBasePressure 通过：%d底面，p=%.9g Pa，真实计算面积 %.9g m²，占总面积 %.3f%%；底热未解。\n', ...
    summary.baseFaces,pb,summary.computedPressureAreaM2,100*summary.computedPressureAreaFraction);
end
function assertRejected(r,code)
% 模块用途：拒绝时全场NaN、零计算面积、可追查底面拒绝状态。
assert(~r.available && ~r.pressureAvailable && ~any(r.computedFace));
assert(all(isnan(r.pressure)) && all(isnan(r.Cp)) && all(isnan(r.heatFlux)));
assert(r.summary.computedPressureArea==0 && r.summary.computedHeatArea==0);
assert(all(r.faceStatus(r.baseFace)==code),r.status);
end
function fs=makeState(T,p,M)
g=1.4;R=287.05287;
fs=struct('temperature',T,'pressure',p,'Mach',M,'gamma',g,'R',R, ...
    'density',p/(R*T),'speed',M*sqrt(g*R*T),'alpha',0,'beta',0, ...
    'totalTemperature',T*(1+.2*M^2),'totalPressure',p*(1+.2*M^2)^3.5);
end
