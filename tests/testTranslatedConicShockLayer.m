function report=testTranslatedConicShockLayer(folder,conicFile)
%TESTTRANSLATEDCONICSHOCKLAYER 释放激波顶点的两个横向位置，并用未参与拟合点审核。
% 模块一：七个全局物理几何参数共用一张激波，不对各面元单独拟合压力。
% 模块二：9个匹配列与36个独立检查列分开，后者不进入优化目标。
% 模块三：采样拟合通过不等于整个区域、驻点或边界层外缘已解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupEngineeringShockLayer;
assert(~exist(folder,'dir'),'ShockLayer:Overwrite','顶点位置研究须在新目录运行。');
d=load(conicFile);geometry=d.geometry;fs=d.fs;seed=d.seed;best=d.bodyFit.result.shape;
upstream=struct('pressure',double(fs.pressure),'temperature',double(fs.temperature), ...
    'density',double(fs.density),'velocityVector',double(fs.flowVelocityVector(:)), ...
    'method','原19752面箭体M5/30km/alpha20均匀自由流；未使用旧压力/外缘/NS状态');
seed.K=best.K;seed.profileB=best.profileB;
seed.standOffM=dot(seed.noseCenter-best.origin,seed.frame(:,1))-seed.referenceRadiusM;
gas=d.cases{2}.column.shock.gas;coordinates=d.coordinates;
timer=tic;bodyFit=fitTranslatedConicShockToBody(double(geometry.V),double(geometry.F),upstream,gas,seed,coordinates, ...
    struct('maximumIterations',10,'bodyMatchToleranceM',1e-5,'quadratureOrder',8));seconds=toc(timer);
mkdir(folder);save(fullfile(folder,'fit_before_checks.mat'),'bodyFit','seconds','-v7');
assert(max(abs(bodyFit.initial.residualM-d.bodyFit.result.residualM))<1e-10);
assert(norm(bodyFit.result.residualM)<norm(d.bodyFit.result.residualM));

%% 二、在三个新半径和十二个周向角上检查，绝不让这些点回流到拟合
shape=bodyFit.result.shape;heldCoordinates=zeros(36,2);k=0;
for r=[.016 .025 .030]
    for angle=(0:11)*pi/6,k=k+1;heldCoordinates(k,:)=[r*cos(angle) r*sin(angle)];end
end
columns=cell(36,1);points=zeros(36,3);normals=points;thickness=zeros(36,1);pressure=thickness;
for j=1:36
    columns{j}=solveEllipticConicShockLayerColumn(upstream,gas,shape,heldCoordinates(j,:)',[1;0],struct('quadratureOrder',16));
    points(j,:)=columns{j}.shock.position';normals(j,:)=columns{j}.shock.normal';
    thickness(j)=columns{j}.normalDistanceM(end);pressure(j)=columns{j}.state.pressure(end);
end
hits=intersectShockNormalsWithBody(double(geometry.V),double(geometry.F),points,normals);
assert(all(hits.hit));heldOut=struct('columns',{columns},'coordinates',heldCoordinates, ...
    'geometricMapping',hits,'residualM',thickness-hits.distanceM,'candidatePressurePa',pressure, ...
    'usedDuringFit',false,'flowValuesAssignedToProductionBody',false);
passed=max(abs(heldOut.residualM))<=bodyFit.bodyMatchToleranceM;
report=struct('status','tested-shock-apex-position-fitting-with-independent-real-body-check-points', ...
    'actualAngleDeg',20,'originalFaces',19752,'globalGeometricParameters',7,'fittedColumns',9,'heldOutColumns',36, ...
    'sampledBodyMatchingConverged',bodyFit.shapeMatchingConverged,'heldOutBodyMatchingPassed',passed, ...
    'bodyMatchToleranceM',bodyFit.bodyMatchToleranceM,'maximumFitResidualM',max(abs(bodyFit.result.residualM)), ...
    'fitRMSM',norm(bodyFit.result.residualM)/3,'maximumHeldOutResidualM',max(abs(heldOut.residualM)), ...
    'heldOutRMSM',norm(heldOut.residualM)/6,'bodyFitSeconds',seconds,'bodyFitEvaluations',bodyFit.evaluations, ...
    'stagnationGlobalMatchingSolved',false,'boundaryLayerEdgeAvailable',false,'heatFluxSolved',false, ...
    'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'translated_conic_layer_evidence.mat'), ...
    'bodyFit','heldOut','seed','coordinates','geometry','fs','gas','report','-v7');
fid=fopen(fullfile(folder,'translated_conic_layer_report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));disp(report);
end
