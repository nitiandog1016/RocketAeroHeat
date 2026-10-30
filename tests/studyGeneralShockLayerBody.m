function report=studyGeneralShockLayerBody(folder,sourceFile)
%STUDYGENERALSHOCKLAYERBODY 在同一原箭体上计算一般共同曲面质量/熵与真实几何残差。
% 模块一：沿用v98七参数端点作几何基线；比较固定整体非椭圆形变，不增加逐点拟合。
% 模块二：9原拟合位置与36独立检查位置均用新连续式映射计算，热量仍未求。
% 模块三：保留原面号/重心和候选端点；没有通过门槛不得赋值到生产壁面。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);setupEngineeringShockLayer;
assert(~exist(folder,'dir'),'ShockLayer:Overwrite','真实一般曲面研究必须在新目录运行。');
d=load(sourceFile);geometry=d.geometry;fs=d.fs;gas=d.gas;
up=struct('pressure',double(fs.pressure),'temperature',double(fs.temperature),'density',double(fs.density), ...
    'velocityVector',double(fs.flowVelocityVector(:)),'method','冻结原箭体M5/30km/alpha20均匀上游；未读取旧壁压或外部CFD');
shape=d.bodyFit.result.shape;shape.deformationDirections=[1 0 sqrt(.5);0 1 sqrt(.5)];
shape.deformationReferenceLengthM=d.seed.referenceRadiusM;shape.deformationWeights=[.12;.07;.05];
coordinateBanks={d.coordinates,d.heldOut.coordinates};groups=cell(2,1);timer=tic;
for k=1:2
    coordinates=coordinateBanks{k};n=size(coordinates,1);columns=cell(n,1);points=zeros(n,3);normals=points;thickness=zeros(n,1);
    for j=1:n
        column=solveGeneralShockLayerColumn(up,gas,shape,coordinates(j,:)',[1;.7;.3;0]);
        points(j,:)=column.shock.position';normals(j,:)=column.shock.normal';thickness(j)=column.normalDistanceM(end);
        % 大回溯网格只在单位核验中保留；真实列保存完整物理量及质量摘要，避免累积缓存。
        m=column.massCoordinateMap;column.massMapSummary=rmfield(m,{'traceParameter','sourceCoordinateGridM', ...
            'coordinateDerivative','logMassRatio','logMassDerivative','coefficientGridM','logProjectedAreaRatio'});
        columns{j}=rmfield(column,'massCoordinateMap');
        fprintf('GENERAL_BODY group=%d column=%d/%d twoPass=%.3g\n',k,j,n,m.twoPassMassRelativeDifference);
    end
    hits=intersectShockNormalsWithBody(double(geometry.V),double(geometry.F),points,normals);
    assert(all(hits.hit),'ShockLayer:BodyFit','固定共同曲面没有全部命中真实原面。');
    groups{k}=struct('coordinates',coordinates,'columns',{columns},'geometricMapping',hits, ...
        'residualM',thickness-hits.distanceM,'usedToFitThisShape',false,'shapeFittedThisStudy',false, ...
        'productionQualified',false,'newProductionSolvedAreaM2',0);
end
report=struct('status','tested-fixed-global-nonelliptic-shock-mass-map-on-real-body-points', ...
    'actualAngleDeg',20,'actualFaces',size(geometry.F,1),'originalFitLocations',9,'independentLocations',36, ...
    'firstGroupMaximumResidualM',max(abs(groups{1}.residualM)),'secondGroupMaximumResidualM',max(abs(groups{2}.residualM)), ...
    'bodyMatchToleranceM',1e-5,'newShapeFits',0,'bodyMatchingPassed', ...
    max(abs([groups{1}.residualM;groups{2}.residualM]))<=1e-5,'seconds',toc(timer), ...
    'productionQualified',false,'boundaryLayerEdgeAvailable',false,'heatFluxSolved',false,'newProductionSolvedAreaM2',0);
mkdir(folder);save(fullfile(folder,'general_body_evidence.mat'),'groups','geometry','fs','shape','gas','up','report','-v7');
fid=fopen(fullfile(folder,'general_body_report.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));disp(report);
end
