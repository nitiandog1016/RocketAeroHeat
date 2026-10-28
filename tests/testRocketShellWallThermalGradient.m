function report=testRocketShellWallThermalGradient(folder,domainFolder,initialFolder)
%TESTROCKETSHELLWALLTHERMALGRADIENT 实际20度新域上的默认全info复现及显式热边界/扩散导数检查。
% 模块一：解析三维斜壁温度场和缺失/矛盾切向输入拒绝，不制造隐含边界数据。
% 模块二：旧默认在同Q/同几何上完整逐字段复现；新策略只变墙热及对应能量残差。
% 模块三：完整空间算子的实际温度扰动与同一热扩散模板对照，不仅比较公式自身。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RocketShellThermalTest:Overwrite','不得覆盖实际测试。');mkdir(folder);

%% 模块一：独立给定的三维仿射场，切向数据来自其解析梯度
n=[1 2 -1]/sqrt(6);known=[20 -3 7];points=[.1 .2 .3;.2 -.1 .4;-.1 .3 .2];delta=[.01 .02 -.005;.02 .01 .005;.015 .01 -.02];
owners=points-delta;normals=repmat(n,3,1);wallT=200+points*known';cellT=200+owners*known';
bc=struct('type','wall','thermalGradient','dirichlet-tangent', ...
    'temperature',@(p,id)200+dot(p,known),'temperatureSurfaceGradient',@(p,id)known-dot(known,n)*n);
[g,audit]=cutCellWallThermalGradient(owners,points,normals,cellT,wallT,bc,(1:3)');
assert(max(abs(g-known),[],'all')<1e-10&&max(abs(sum(g.*delta,2)-(wallT-cellT)))<1e-12);
bad={rmfield(bc,'temperatureSurfaceGradient'), ...
    struct('type','wall','temperature',300,'thermalGradient','legacy','temperatureSurfaceGradient',[0 0 0]), ...
    struct('type','wall','temperature',300,'thermalGradient','wrong'), ...
    struct('type','wall','temperature',300,'thermalGradient','dirichlet-tangent','temperatureSurfaceGradient',n), ...
    struct('type','wall','temperature',300,'thermalGradient','dirichlet-tangent','temperatureSurfaceGradient',[2 -1 0])};
for j=1:numel(bad)
    rejected=false;try,cutCellWallThermalGradient(owners,points,normals,cellT,wallT,bad{j},(1:3)');catch ex,rejected=strcmp(ex.identifier,'CutNS:WallThermalBoundary');end
    assert(rejected,'RocketShellThermalTest:Boundary','缺失或矛盾墙边界未被拒绝。');
end
rejected=false;try,cutCellWallThermalGradient(points+delta,points,normals,cellT,wallT,bc,(1:3)');catch ex,rejected=strcmp(ex.identifier,'CutNS:WallThermalGeometry');end
assert(rejected);

%% 约束预分配只改变存储，独立保留原动态装配作为小几何完整矩阵参考
v=[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];f=convhull(v);
normal=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);flip=sum(normal.*v(f(:,1),:),2)<0;f(flip,:)=f(flip,[1 3 2]);
small=buildRocketTriangularShellDomain(v,f,struct('normalLayerHeightsM',[0 .01],'outerScaleFactors',2));
gasSmall=struct('gamma',1.4,'R',287,'mu',1.8e-5,'Pr',.72,'temperatureRange',[80 2000]);
p=[1 -100 0 25 300];q=[p(1),p(1)*p(2:4),p(1)*(287/.4*p(5)+.5*sum(p(2:4).^2))];q=repmat(q,numel(small.cells.volume),1);
given=struct('type','characteristic','mode','riemannFarfield','reservoirPrimitive',p);
bs=struct('xmin',given,'xmax',given,'ymin',given,'ymax',given,'zmin',given,'zmax',given, ...
    'solidWall',struct('type','wall','temperature',320,'velocity',[0 0 0]));
for mode=1:2
    if mode==2,bs.solidWall.thermalGradient='dirichlet-tangent';end
    tag=cutCellBoundaryStencilState(q,small,gasSmall,bs);fast=buildCutCellDiffusionStencil(small,tag);reference=rocketShellStencilAllocationReference(small,tag);
    assert(isequaln(fast,reference),'RocketShellThermalTest:Allocation','预分配改变了完整扩散/应力模板。');
end
fprintf('SHELL_THERMAL_ALLOCATION_REFERENCE_PASSED\n');

%% 模块二：真实原网格/原初态的完整默认空间残差、info与所有护栏保持
d=load(fullfile(domainFolder,'shell_domain_evidence.mat'),'domain');cut=d.domain;
old=load(fullfile(initialFolder,'initial_residual_evidence.mat'),'Q','rate','info','gas','boundary');Q=old.Q;gas=old.gas;
[legacyRate,legacyInfo]=cutCellNavierStokesResidual(Q,cut,gas,old.boundary);
assert(isequaln(legacyRate,old.rate)&&isequaln(legacyInfo,old.info),'RocketShellThermalTest:Legacy','默认完整空间结果改变。');
fprintf('SHELL_THERMAL_REAL_LEGACY_INFO_PASSED\n');
clear legacyRate legacyInfo;
boundary=old.boundary;boundary.solidWall.thermalGradient='dirichlet-tangent';
t=tic;[rateNew,infoNew]=cutCellNavierStokesResidual(Q,cut,gas,boundary);evaluationSeconds=toc(t);
assert(isequal(rateNew(:,1:4),old.rate(:,1:4))&&isequal(infoNew.faceFlux,old.info.faceFlux)&& ...
    isequal(infoNew.gradient,old.info.gradient)&&isequal(infoNew.wallPressure,old.info.wallPressure)&& ...
    isequal(infoNew.wallFlux(:,1:4),old.info.wallFlux(:,1:4)), ...
    'RocketShellThermalTest:Scope','新热边界意外改变质量/动量/内面通量/WLS/墙压。');
assert(all(infoNew.wallHeatFluxIntoWall<0)&&nnz(old.info.wallHeatFluxIntoWall>0)==184, ...
    'RocketShellThermalTest:Direction','冷气体/热等温壁的导热方向未正确恢复。');
mapped=mapCutCellWallFields(cut,infoNew,false);assert(~mapped.eligibleForSurfaceCoverage&&mapped.heatPowerMappingDifference==0);
balance=max(abs(infoNew.balanceResidual)./max(1,abs(infoNew.volumeIntegratedResidual)+abs(infoNew.netBoundaryFlux)));
assert(balance<1e-9,'RocketShellThermalTest:Balance','新热边界不满足原五量账。');
fprintf('SHELL_THERMAL_REAL_NEW_WALL_AND_BALANCE_PASSED\n');

%% 模块三：给定边界冻结时，对靠墙真实单元的温度扰动检查热模板导数
state=cutCellBoundaryStencilState(Q,cut,gas,boundary);t=tic;stencil=buildCutCellDiffusionStencil(cut,state);stencilSeconds=toc(t);
fprintf('SHELL_THERMAL_REAL_STENCIL_BUILT_SECONDS=%.9g\n',stencilSeconds);
perturbation=zeros(size(Q,1),1);ids=[1 6 100 14169 14170 14171 14260 15000 19000];perturbation(ids)=[.5 -.25 .7 .3 -.5 .8 .6 -.2 .4];
Qperturbed=Q;cv=gas.R/(gas.gamma-1);Qperturbed(:,5)=Qperturbed(:,5)+Q(:,1)*cv.*perturbation;
[~,infoPerturbed]=cutCellNavierStokesResidual(Qperturbed,cut,gas,boundary);
assert(isequal(infoPerturbed.primitive(:,2:4),infoNew.primitive(:,2:4)));
deltaFlux=infoPerturbed.viscousFaceFlux(:,5)-infoNew.viscousFaceFlux(:,5);faces=cut.openFaces;
actual=accumarray(faces.cellMinus,faces.area.*deltaFlux,[size(Q,1) 1]);active=faces.cellPlus>0;
actual=actual-accumarray(faces.cellPlus(active),faces.area(active).*deltaFlux(active),[size(Q,1) 1]);
actual=actual-accumarray(cut.walls.ownerCell,cut.walls.area.*(infoPerturbed.wallHeatFluxIntoWall-infoNew.wallHeatFluxIntoWall),[size(Q,1) 1]);
actual=actual./cut.cells.volume;
heatGeometryMatrix=stencil.heatGeometricMatrix;k=gas.mu*gas.gamma*gas.R/(gas.gamma-1)/gas.Pr;
expected=k*(heatGeometryMatrix*perturbation);derivativeRelative=max(abs(actual-expected))/max(1,max(abs(expected)));
assert(derivativeRelative<1e-8,'RocketShellThermalTest:Stencil','实际完整空间导热扰动与步长模板不一致。');
fprintf('SHELL_THERMAL_REAL_KERNEL_DERIVATIVE_PASSED=%.9g\n',derivativeRelative);
report=struct('status','passed-actual-20-degree-legacy-info-explicit-wall-temperature-and-actual-diffusion-derivative', ...
    'cellCount',size(Q,1),'sourceFaces',numel(cut.walls.area),'defaultCompleteInfoExactlyReplayed',true, ...
    'legacyInitialHeatWrongDirectionFaces',184,'newInitialHeatWrongDirectionFaces',0, ...
    'newInitialHeatRangeWm2',[min(infoNew.wallHeatFluxIntoWall),max(infoNew.wallHeatFluxIntoWall)], ...
    'fullMassMomentumAndInteriorFluxExactlyKept',true,'fiveQuantityBalanceRelative',balance, ...
    'actualHeatDerivativeRelativeDifference',derivativeRelative,'oneNewSpatialEvaluationSeconds',evaluationSeconds, ...
    'diffusionStencilBuildSeconds',stencilSeconds,'wallTemperatureValueKind','prescribed', ...
    'boundaryRejections',numel(bad)+1,'steadyConverged',false,'physicalSolution',false, ...
    'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'wall_thermal_evidence.mat'),'report','rateNew','infoNew','infoPerturbed','perturbation', ...
    'heatGeometryMatrix','boundary','gas','audit','-v7');
fid=fopen(fullfile(folder,'wall_thermal_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
