function report=testRocketKrylovFailureEvidence(folder)
%TESTROCKETKRYLOVFAILUREEVIDENCE 原增强内核完整回放及失败观察器的非放行性质。
% 模块一：缺省与显式观察器成功路径都逐字段复现旧80维完整输出。
% 模块二：故意非线性作用仍由原Arnoldi门槛拒绝，额外保存证据而不返回方向。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'KrylovFailureTest:Overwrite','不得覆盖失败观察器测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','recycledkernel01','recycled_kernel_evidence.mat');s=load(sourceFile,'A','b','x','recycled');
[x,detail]=solveRocketRecycledGMRES(@(v)s.A*v,@(v)v,s.b,s.recycled.options);
assert(isequaln(x,s.x)&&isequaln(detail,s.recycled),'KrylovFailureTest:Default','缺省完整增强内核输出没有复现。');
options=s.recycled.options;options.failureEvidenceFile=fullfile(folder,'must_not_exist.mat');
[explicit,explicitReport]=solveRocketRecycledGMRES(@(v)s.A*v,@(v)v,s.b,options);
assert(isequaln(explicit,x)&&isequaln(explicitReport,detail)&&exist(options.failureEvidenceFile,'file')~=2);
bad=struct('restartLength',4,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',4, ...
    'absoluteTolerance',1e-12,'maximumBasisBytes',16*1024^2,'saveBasisEvidence',false, ...
    'failureEvidenceFile',fullfile(folder,'intentional_nonlinear_failure.mat'));
caught=false;try,solveRocketRecycledGMRES(@(v)diag(1:8)*v+.1*norm(v)^2*ones(8,1),@(v)v,ones(8,1),bad);
catch ex,caught=strcmp(ex.identifier,'RecycledGMRES:Arnoldi');end
assert(caught&&exist(bad.failureEvidenceFile,'file')==2);f=load(bad.failureEvidenceFile,'failure');
assert(f.failure.orthogonalityThreshold==1e-8&&f.failure.modelRelationThreshold==1e-8&& ...
    ~f.failure.usableDirectionClaimed&&f.failure.modelRelationRelativeDifference>=1e-8);
report=struct('status','passed-default-complete-kernel-replay-and-explicit-failure-evidence-with-original-rejection', ...
    'oldCompleteKernelOutputExactlyReplayed',true,'observerSuccessOutputExactlyReplayed',true,'intentionalNonlinearActionRejected',true, ...
    'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'failure_observer_unit_evidence.mat'),'report','sourceFile','x','detail','explicit','explicitReport','bad','-v7');
fid=fopen(fullfile(folder,'failure_observer_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
