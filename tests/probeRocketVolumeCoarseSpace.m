function report=probeRocketVolumeCoarseSpace(folder)
%PROBEROCKETVOLUMECOARSESPACE 只检查真实域聚合规模、连通及真实体积传递，不作NS调用。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RocketCoarseProbe:Overwrite','不得覆盖粗空间探针。');
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');
space=buildRocketVolumeCoarseSpace(d.domain,struct('cartesianBins',[2 2 2],'layerBreaks',4,'maximumCoarseUnknowns',200));
map=space.cellToGroup;R=space.restriction;P=space.prolongation;fineVolume=space.fineVolume;groupVolume=space.groupVolume;
report=struct('status','measured-real-connected-volume-coarse-space','fineCells',numel(map),'coarseGroups',space.groupCount, ...
    'coarseUnknowns',space.coarseUnknowns,'restrictionProlongationIdentityError',space.restrictionProlongationIdentityError, ...
    'totalVolumeRelativeDifference',space.totalVolumeRelativeDifference,'options',space.options, ...
    'newOriginalNSEvaluations',0,'newModelApplications',0,'physicalMeshChanged',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'volume_coarse_space_evidence.mat'),'report','map','R','P','fineVolume','groupVolume','-v7');
fid=fopen(fullfile(folder,'volume_coarse_space_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
