function [mesh,pair]=attachShortSupportGridPair(mesh,report,geometry)
% 模块一：将已核验的五层基线与十层原面按同一编号配对，保留两份独立来源。
root=fileparts(mfilename('fullpath'));folder=fullfile(root,'cases','short_sphere');
pair=readShortSupportGridPair(folder);
assert(strcmp(report.frameSHA256,pair.fineFrame) && strcmp(report.thermal.meshSHA256,pair.fineMesh) && ...
    strcmp(report.bodySHA256,pair.bodySHA256) && report.faces==pair.faces);
archive=fullfile(folder,'supportedSurface65.zip');assert(strcmp(cm011FileSHA256(archive),pair.coarseArchiveSHA256));
temp=tempname;mkdir(temp);done=onCleanup(@()discard(temp)); %#ok<NASGU>
unzip(archive,temp);s=jsondecode(fileread(fullfile(temp,'summary.json')));
assert(strcmp(s.frameSHA256,pair.coarseFrame) && strcmp(s.thermal.meshSHA256,pair.coarseMesh) && ...
    strcmp(s.transferSHA256,report.transferSHA256) && strcmp(s.caseSHA256,report.caseSHA256));
wallFile=fullfile(temp,'wall.csv');assert(strcmp(cm011FileSHA256(wallFile),s.files.wall_csv));
w=readtable(wallFile);
assert(isequal(w.faceID,(1:pair.faces)') && isequal(w.regionID,geometry.regionID) && ...
    max(abs(w.areaM2-geometry.areaM2))<1e-17 && isequal(w.wallInputK,mesh.hybridWallTemperatureK));
mesh.shortSupportGridPair=struct('coarseFrame',pair.coarseFrame,'fineFrame',pair.fineFrame, ...
    'coarseMesh',pair.coarseMesh,'fineMesh',pair.fineMesh,'coarseHeat',w.fluxInWm2,'coarseGauge',w.gaugePa, ...
    'fluxRange',pair.sharedFluxRangeKWm2,'gaugeRange',pair.sharedGaugeRangePa);
mesh.hybridCondition.caseLabel='Short Test 65 · 双支管银球 · 十层近壁网主值';
mesh.evidenceIdentity.physicalCondition=mesh.hybridCondition;
mesh.evidenceIdentity.scope='十层近壁网研究解；五层与十层共用原面、外域、给定壁温，三维对照差值为十层减五层。';
end
function discard(folder)
% 仅释放本次展开的已登记文件。
for name={'summary.json','wall.csv','stations.csv'}
    path=fullfile(folder,name{1});if isfile(path),delete(path);end
end
rmdir(folder);
end
