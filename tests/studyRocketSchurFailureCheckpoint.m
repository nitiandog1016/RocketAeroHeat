function report=studyRocketSchurFailureCheckpoint(folder)
%STUDYROCKETSCHURFAILURECHECKPOINT 仅重建一次原候选以保存此前缺失的拒绝主元/组装检查点。
% 模块一：使用冻结C与端点，原算术/门槛保持；今后直接读取检查点，不再重复组装。
% 模块二：结构化失败摘要包含原CV/面/层，不依赖被宿主破坏的中文异常字符串。
% 模块三：本研究不执行右端迭代、新J/原NS或物理更新。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SchurCheckpointStudy:Overwrite','不得覆盖真实检查点研究。');mkdir(folder);
packedFile=fullfile(root,'outputs','v83_validation','blockilureal01','steady_graph_packed_evidence.mat');s=load(packedFile,'packed','sourceFile','sourceBinding');packed=s.packed;sourceFile=s.sourceFile;sourceBinding=s.sourceBinding;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceFile),sourceBinding.sourceFileSHA256));
sourceBinding.packedFile=packedFile;sourceBinding.packedFileSHA256=rocketStreamingFileSHA256(packedFile);sourceBinding.sourceFaces=packed.sourceFaces;sourceBinding.layers=packed.layers;
observation=struct('folder',fullfile(folder,'checkpoint'),'sourceBinding',sourceBinding);caught=false;identifier='';
try,pre=buildRocketFieldSchurPreconditioner(packed,observation);save(fullfile(folder,'unexpected_success.mat'),'pre','sourceBinding','-v7');
catch ex,identifier=ex.identifier;if ~strcmp(identifier,'FieldILU:FactorDefect'),rethrow(ex);end;caught=true;end
assert(caught,'SchurCheckpointStudy:Unexpected','原拒绝未复现，不能冒称相同数值路径。');
file=fullfile(observation.folder,'density_energy_failure.mat.json');failureSummary=jsondecode(fileread(file));
assert(failureSummary.algebraicRow==188887&&failureSummary.originalControlVolume==61925&&failureSummary.originalFace==2669&&failureSummary.layer==4&& ...
    failureSummary.originalInverseDefectThreshold==1e-8&&failureSummary.maximumInverseDefect>=1e-8&&~failureSummary.preconditionerUsable);
report=struct('status','reproduced-original-Schur-rejection-with-raw-assembly-and-structured-pivot-ledger', ...
    'errorIdentifier',identifier,'failureSummary',failureSummary,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
fid=fopen(fullfile(folder,'schur_failure_checkpoint_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
