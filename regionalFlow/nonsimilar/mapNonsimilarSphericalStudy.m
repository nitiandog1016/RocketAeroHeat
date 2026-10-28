function mapped=mapNonsimilarSphericalStudy(mesh,study,freestream,cfg)
%MAPNONSIMILARSPHERICALSTUDY 将条件非相似站映回当前真实球鼻三角面。
% 模块用途：核验当前几何和工况后，仅在完整三角形落于已接受角度段时
% 线性插值站值。研究数值不能写入生产heating或计作新物理解面积。
% 不使用旧95/186计数，不跨未解角度，不补驻点/后体/恢复温度或壁温解。
if nargin<4,cfg=struct;end
if ~isstruct(study)||~isfield(study,'kind')|| ...
        ~strcmp(study.kind,'conditional-spherical-nonsimilar-study')|| ...
        ~study.currentSourceReplayPassed||study.physicalSolution||study.eligibleForSurfaceCoverage
    error('RocketAeroHeat:NonSimilarMap:Study','需要已重放核对且保持条件研究身份的数据。');
end
if ~isfield(study,'replayBindingSHA256')||~strcmp(study.replayBindingSHA256,nonsimilarStudyDigest(study))
    error('RocketAeroHeat:NonSimilarMap:Binding','条件/站表自重放后改变，不能跨缺站或使用旧缓存。');
end
keys={'Mach','pressure','temperature','speed','R','gamma','alpha','beta'};
for k=1:numel(keys)
    key=keys{k};a=study.freestream.(key);b=freestream.(key);
    if ~isnumeric(b)||~isreal(b)||~isscalar(b)||~isfinite(b)|| ...
            abs(a-b)>1e-10*max([abs(a),abs(b),1])
        error('RocketAeroHeat:NonSimilarMap:Condition','目标工况与已保存条件段不同：%s。',key);
    end
end
[mode,~]=selectAirThermalModel(freestream,cfg);
if ~strcmp(mode,'thermallyPerfect')||~strcmp(study.gas.thermalModel,'thermallyPerfect')|| ...
        study.gas.R~=freestream.R||study.gas.gamma~=freestream.gamma
    error('RocketAeroHeat:NonSimilarMap:ThermalModel','不能用定比热配置承接NASA7条件研究。');
end
flow=freestream.unitFlowDirection;
if norm(flow-[-1 0 0])>1e-12||freestream.alpha~=0||freestream.beta~=0
    error('RocketAeroHeat:NonSimilarMap:Incidence','该轴对称条件段不能映射到非零攻角/侧滑。');
end
Tw=surfaceWallTemperature(size(mesh.F,1),cfg);
if any(~isfinite(Tw))||any(abs(Tw-study.uniformWallTemperature)>1e-10)
    error('RocketAeroHeat:NonSimilarMap:WallHistory','必须匹配整段给定均匀壁温历史，不能只借当前局部壁温。');
end
% 当前V/F必须与几何缓存逐值相符，避免拿修改后的顶点配旧法向/区域。
V=mesh.V;F=mesh.F;nf=size(F,1);
actualCenter=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
crossing=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
actualNormal=crossing./vecnorm(crossing,2,2);g=mesh.localGeometry;
if ~isequal(size(g.center),size(actualCenter))||~isequal(size(g.normal),size(actualNormal))|| ...
        any(~isfinite(g.center(:)))||any(~isfinite(g.normal(:)))|| ...
        max(abs(g.center(:)-actualCenter(:)))>1e-10*max([abs(V(:));1])|| ...
        max(abs(g.normal(:)-actualNormal(:)))>1e-8
    error('RocketAeroHeat:NonSimilarMap:StaleGeometry','当前顶点/三角面与局部几何缓存不一致。');
end
edge=solveSphericalNoseEdge(mesh,freestream,cfg);
if ~edge.available||abs(edge.noseRadius/study.noseRadius-1)>1e-6
    error('RocketAeroHeat:NonSimilarMap:Geometry','当前几何未核验为同半径的球鼻。');
end
angles=[study.station.angleDeg]';
if numel(angles)<2||~isreal(angles)||any(~isfinite(angles))||any(diff(angles)<=0)|| ...
        any(diff(angles)>study.maximumAngularStepDeg+1e-9)|| ...
        ~all([study.station.thinLayerChecksPassed])||~all([study.station.positiveWallShear])
    error('RocketAeroHeat:NonSimilarMap:Stations','角度表非连续有效接受站。');
end
% 模块一：整面几何必须在有效角度段内，不能仅靠质心掩盖跨边界三角形。
angularBounds=sphericalTrianglePolarBounds(V,F,edge.sphereCenter);
minAngle=angularBounds(:,1);maxAngle=angularBounds(:,2);
inInterval=minAngle>=angles(1)-1e-9 & maxAngle<=angles(end)+1e-9;
theta=edge.polarAngleDeg;
keep=edge.geometryFace & edge.algebraicallyValid & inInterval & ...
    isfinite(theta) & theta>=angles(1) & theta<=angles(end);
reason=repmat({'几何/球鼻外缘不适用'},nf,1);
reason(edge.geometryFace & ~inInterval)={'完整三角形不在已接受非相似角度段；不外推或跨入口/边界'};
reason(keep)={''};
blank=nan(nf,1);mapped=struct('kind','conditional-spherical-face-map','schemaVersion',1, ...
    'available',any(keep),'faceMask',keep,'heatFlux',blank,'wallShearStress',blank, ...
    'edgePressure',blank,'edgeTemperature',blank,'edgeSpeed',blank, ...
    'wallTemperatureInput',Tw,'wallTemperatureSolved',blank,'recoveryTemperature',blank, ...
    'energyDefectOverWallPower',blank,'polarAngleDeg',theta, ...
    'profileNodeIDs',nan(nf,2),'profileNodeWeight',blank,'missingReason',{reason}, ...
    'method',study.method,'scope',study.scope, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'experimentalValidation',false, ...
    'newSolvedAreaM2',0,'noseRadius',edge.noseRadius,'sphereCenter',edge.sphereCenter, ...
    'acceptedIntervalDeg',[angles(1),angles(end)],'inletAngleDeg',study.inletAngleDeg, ...
    'source',study.source,'currentSourceHashes',study.currentSourceHashes, ...
    'condition',freestream,'geometryVertices',V,'geometryFaces',F, ...
    'sourceNodes',study.station, ...
    'interpolation','真实面质心极角线性插值，完整平面三角形均在有效段内；不外推', ...
    'angularBoundingMethod','正X投影三角形到原点的精确距离；包括边内部和轴线交点，不仅检查顶点', ...
    'summary',struct);
fields={'heatFlux','wallShearStress','edgePressure','edgeTemperature','edgeSpeed','energyDefectOverWallPower'};
for k=1:numel(fields)
    key=fields{k};values=[study.station.(key)]';
    if ~isreal(values)||any(~isfinite(values)),error('RocketAeroHeat:NonSimilarMap:Finite','接受站包含非有限/复数%s。',key);end
    mapped.(key)(keep)=interp1(angles,values,theta(keep),'linear');
end
for id=find(keep)'
    lower=find(angles<=theta(id),1,'last');lower=min(lower,numel(angles)-1);
    mapped.profileNodeIDs(id,:)=[lower,lower+1];
    mapped.profileNodeWeight(id)=(theta(id)-angles(lower))/(angles(lower+1)-angles(lower));
end
% 模块二：面积与热功率只描述条件研究映射，不是生产已解覆盖。
crossing=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
area=vecnorm(crossing,2,2)/2;
mapped.summary=struct('mappedFaces',nnz(keep),'totalFaces',nf, ...
    'mappedConditionalAreaM2',sum(area(keep)), ...
    'mappedConditionalAreaFraction',sum(area(keep))/sum(area), ...
    'conditionalHeatPowerW',sum(area(keep).*mapped.heatFlux(keep)), ...
    'newProductionSolvedAreaM2',0,'minimumFullFaceAngleDeg',NaN,'maximumFullFaceAngleDeg',NaN);
if any(keep)
    mapped.summary.minimumFullFaceAngleDeg=min(minAngle(keep));
    mapped.summary.maximumFullFaceAngleDeg=max(maxAngle(keep));
end
end
