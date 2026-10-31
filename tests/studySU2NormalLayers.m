function studySU2NormalLayers(out)
%STUDYSU2NORMALLAYERS 核对原面不变的近壁偏置层与原生方向，先只做几何研究。
% 模块用途：比较同层数/同首层参数的近壁方向，不把几何通过算作流场资格。
assert(~isfolder(out));mkdir(out);
[surface,truth]=buildRocketTestBody(struct('resolution',24,'protrusionEnabled',false));
options=struct('radialLayers',40,'firstRadialSpacing',2e-6,'maxWallAngleStepDeg',10,'baseRadialSegments',4, ...
    'nearWall',struct('normalLayers',24,'normalGrowth',1.25));
volume=buildSU2RocketVolume(surface,truth,options);
q=volume.meshConstruction;
assert(isequal(surface.V,volume.sourceVertices) && isequal(surface.F,volume.sourceFaces));
assert(volume.topology.declaredBoundaryExactlyMatched && all(volume.cellVolume>0));
assert(q.minimumWallAlignmentDeg>q.baselineMinimumWallAlignmentDeg+20);
projection=buildSU2WallProjection(volume,surface);
assert(projection.maximumTargetAreaRelativeError<2e-8 && projection.maximumSourceAreaRelativeError<2e-8);
info=writeSU2RocketCase(volume,fullfile(out,'native_input'),struct('Mach',5,'temperature',226.51,'pressure',1171.9, ...
    'alpha',20,'beta',0,'wallTemperature',300,'iterations',2));
save(fullfile(out,'normal_volume.mat'),'surface','truth','options','volume','projection','info','-v7');
report=struct('nodes',size(volume.vertices,1),'cells',numel(volume.cells),'originalFaces',size(surface.F,1), ...
    'wallArea',volume.bodyArea,'nearWall',q,'minimumSubTetVolume',volume.minimumSubTetVolume, ...
    'targetMappingClosure',projection.maximumTargetAreaRelativeError,'physicsQualified',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'geometry_report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Normal-layer complete rocket: %d cells, wall alignment min %.6g deg versus ray %.6g; first distance %.6g m; no physical qualification.\n', ...
    numel(volume.cells),q.minimumWallAlignmentDeg,q.baselineMinimumWallAlignmentDeg,q.firstNormalSpacing);
end
