function report=extractRocketFlowOrderingGeometry(folder)
%EXTRACTROCKETFLOWORDERINGGEOMETRY 从已验证真实域提取完整质心/给定来流供代数排序。
% 模块一：全部CV和原墙编号保持，原域来源按字节绑定，不重新生成几何。
% 模块二：给定U∞投影只用于排序，不冒称当地流线、外缘速度或物理解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FlowOrderGeometry:Overwrite','不得覆盖真实排序几何。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat');s=load(sourceFile,'domain','fs');cut=s.domain;fs=s.fs;clear s;
centroid=cut.cells.centroid;volume=cut.cells.volume;sourceFaces=numel(cut.walls.originalFaceID);layers=cut.quality.layerCount;
assert(isequal(cut.walls.originalFaceID,(1:sourceFaces)')&&size(centroid,1)==sourceFaces*layers&&all(volume>0));
unitFlowDirection=fs.unitFlowDirection(:);flowVelocityVector=fs.flowVelocityVector(:);
assert(norm(unitFlowDirection-flowVelocityVector/norm(flowVelocityVector))<1e-12);clear cut;
report=struct('status','extracted-source-bound-real-cell-centroids-and-prescribed-freestream-direction', ...
    'cellCount',size(centroid,1),'sourceFaces',sourceFaces,'layers',layers,'sourceFile',sourceFile, ...
    'sourceSHA256',rocketStreamingFileSHA256(sourceFile),'givenUnitFlowDirection',unitFlowDirection', ...
    'localSolvedFlowClaimed',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'flow_ordering_geometry.mat'),'report','centroid','volume','unitFlowDirection','flowVelocityVector','sourceFile','-v7');
fid=fopen(fullfile(folder,'flow_ordering_geometry_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
