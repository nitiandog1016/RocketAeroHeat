function report=check_resolved_export(folder,checkFolder)
%CHECK_RESOLVED_EXPORT 按独立环带 CSV 核对真正送给前端的值，并检查两条网格入口。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
evidence=fullfile(root,'docs','validation','v263','evidence','tmr_resolved');
wall=readtable(fullfile(evidence,'wall.csv'));
assert(~isfolder(checkFolder));mkdir(checkFolder);

% 模块一：浏览器离线数据与 JSON 逐字一致，所有显示三角形均能回到真实环带。
payload=fileread(fullfile(folder,'results.json'));
assert(strcmp(strtrim(fileread(fullfile(folder,'case-data.js'))),['window.ROCKET_CASE=' payload ';']));
d=jsondecode(payload);b=d.wallBinding;fields=d.fields;
[found,ring]=ismember(b.sourceFaceID,wall.faceID);n=numel(ring);
assert(n==4992 && all(found) && all(accumarray(ring,1)==96) && height(wall)==52);
assert(isequal(b.ownerCell,wall.nearWallCellID(ring)));
assert(all(abs(accumarray(ring,b.curvedAreaPerTriangleM2)-wall.ringArea_m2)<1e-12));
assert(d.surfaceIntegralAudit.endcapFaces==0 && ~d.researchIdentity.physicalSolution && ...
    ~d.researchIdentity.eligibleForSurfaceCoverage && ~d.researchIdentity.experimentalValidation);
assert(abs(sum(b.curvedAreaPerTriangleM2.*b.sourceHeatFluxIntoWallWm2)- ...
    sum(wall.ringArea_m2.*wall.heatIntoWall_W_m2))<1e-7);
source={wall.pressure_Pa(ring)/1000,wall.heatIntoWall_W_m2(ring)/1000,wall.givenTw_K(ring),wall.nearWallCellT_K(ring)};
keys={'pressure','heatFlux','wallInput','nearWallCellTemperature'};status=[1 1 2 1];
for k=1:numel(keys)
    f=fields(strcmp({fields.key},keys{k}));
    assert(isscalar(f) && all(f.status==status(k)) && all(abs(f.values-source{k})<=1e-10) && ...
        all(f.methodIndex==0) && ~isempty(f.methodNames) && numel(f.scopeNames)==height(wall));
end
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    f=fields(strcmp({fields.key},key{1}));
    assert(all(isnan(f.values)) && all(f.status==0));
end

% 模块二：尺寸相同的非官方内点须在写页面之前拒绝。
s=load(fullfile(folder,'tmr_surface_snapshot.mat'),'mesh','result','boundary');
bad=s.mesh;bad.X(3,4)=bad.X(3,4)+1e-7;
badFolder=fullfile(checkFolder,'invalid-grid');
try
    exportTmrAxisymmetricSAStudy(bad,s.result,s.boundary,badFolder);
    error('TmrStations:UnexpectedAdmission','改动后的节点未被拒绝。');
catch err
    assert(strcmp(err.identifier,'AxisFV:TmrSAViewGeometry') && ~isfolder(badFolder));
end

% 模块三：用真实旧稳态经过原入口，核对规则网格身份和原字段。
legacy=load(fullfile(root,'outputs','v66_validation','tmr_wall_thermal','tmr_wall_thermal_evidence.mat'), ...
    'mesh','implicit','boundary');
old=exportTmrAxisymmetricSAStudy(legacy.mesh,legacy.implicit,legacy.boundary, ...
    fullfile(checkFolder,'legacy'),struct('azimuthSegments',8));
assert(strcmp(old.gridIdentity.kind,'official-TMR-level3-nested-original-nodes-v1') && ...
    old.gridIdentity.streamwiseStride==2 && old.gridIdentity.normalStride==2 && old.audit.sourceWallFaces==40);
oldData=jsondecode(fileread(old.jsonFile));assert(numel(oldData.fields)==6);
report=struct('sourceRings',height(wall),'displayFaces',n,'matchedQuantities',{keys}, ...
    'givenTwK',unique(wall.givenTw_K),'unknownTemperatureFields',3,'curvedWallHeatW',sum(wall.ringArea_m2.*wall.heatIntoWall_W_m2), ...
    'nonofficialNodeRejected',true,'legacyGridReplayed',true,'offlineDataMatchesJSON',true, ...
    'visualInteractionTested',false,'newMainCaseAreaM2',0);
fid=fopen(fullfile(checkFolder,'report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
end
