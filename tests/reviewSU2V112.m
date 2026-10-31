function reviewSU2V112(out)
%REVIEWSU2V112 核对多重网格原面结果，并检查旧成功快照仍可在主工作台复开。
% 不启动求解；拒绝把原生终态、映射通过或旧UI复核当成物理资格。
assert(~isfolder(out));mkdir(out);
s=load('outputs/v108_validation/normal01/normal_volume.mat','projection');
folder='outputs/v112_validation/mg01';data=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
result=mapSU2WallResult(data,s.projection);result.receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
audit=auditSU2ResultConvergence(folder,result);
assert(audit.checkpointIdentityConfirmed,audit.reason);
assert(audit.checkpointIteration==599 && ~audit.checkpointResidualPassed && ~result.productionQualified);
save(fullfile(out,'mg_mapped.mat'),'result','audit','-v7');
testSU2ResearchReplay(fullfile(out,'legacy'),'docs/validation/v107/replay_seed/mapped_research.mat',false);
report=struct('mgOriginalFaces',numel(result.heatFluxRaw),'mgAudit',audit,'legacySnapshotStillAccepted',true, ...
    'solverInvoked',false,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Multigrid same-face map/audit and legacy successful snapshot checks passed; no solver rerun.\n');
end
