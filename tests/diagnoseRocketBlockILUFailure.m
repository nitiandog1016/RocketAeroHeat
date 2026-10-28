function report=diagnoseRocketBlockILUFailure(folder)
%DIAGNOSEROCKETBLOCKILUFAILURE 复用已保存因子和已验证同Q剩余，只观察原Arnoldi拒绝。
% 模块一：不重构因子/粗像，也不调用新的完整J/NS；所有输入来源单列绑定。
% 模块二：观察器只保存原1e-8门槛下的失败基、实际剩余和投影剩余。
% 模块三：拒绝候选仍不可用于物理解或覆盖，不以放宽阈值消除停滞。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'BlockFailure:Overwrite','不得覆盖失败诊断。');mkdir(folder);
factorFile=fullfile(root,'outputs','v83_validation','blockilureal01','block_ilu_factor_evidence.mat');s=load(factorFile,'factor','Q','sourceFile','sourceBinding','sourceIdentity');
factor=s.factor;Q=s.Q;sourceFile=s.sourceFile;sourceBinding=s.sourceBinding;sourceIdentity=s.sourceIdentity;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceFile),sourceBinding.sourceFileSHA256));
for key={'packRocketSteadyGraph','factorRocketBlockILU0','solveRocketBlockILU0'}
    assert(strcmp(rocketStreamingFileSHA256(which(key{1})),sourceIdentity.(key{1})),'BlockFailure:Source','原因子数学来源已变化。');
end
anchorFile=fullfile(root,'outputs','v83_validation','innergraphreal01','inner_graph_real_evidence.mat');a=load(anchorFile,'Q','remaining','sourceFile');
assert(isequaln(Q,a.Q)&&strcmp(sourceFile,a.sourceFile));remaining=a.remaining;clear a;
s=load(sourceFile,'result');assert(isequaln(s.result.Q,Q));graph=s.result.firstOperator;clear s;
C=graph;for k=1:5,C.diagonal(k,k,:)=C.diagonal(k,k,:)-1/C.pseudoTimeStepSeconds;end
[single,singleReport]=solveRocketBlockILU0(factor,remaining);Csingle=applyRocketWlsNormalOperator(C,single);clear C;
save(fullfile(folder,'single_ilu_failure_context.mat'),'factorFile','anchorFile','sourceFile','sourceBinding','remaining','single','singleReport','Csingle','-v7');
options=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false, ...
    'failureEvidenceFile',fullfile(folder,'arnoldi_failure_evidence.mat'));
reproduced=false;identifier='';
try,[candidate,innerReport]=solveRocketBlockGraphKrylov(graph,factor,remaining,options);save(fullfile(folder,'unexpected_replay_result.mat'),'candidate','innerReport','-v7');
catch ex,identifier=ex.identifier;if ~strcmp(identifier,'RecycledGMRES:Arnoldi'),rethrow(ex);end;reproduced=true;end
assert(reproduced&&exist(options.failureEvidenceFile,'file')==2,'BlockFailure:NotReproduced','原拒绝未复现，需要另行比较输入/舍入，不能冒称一致。');
f=load(options.failureEvidenceFile,'failure');failure=f.failure;
report=struct('status','reproduced-original-real-block-ILU-Arnoldi-rejection-with-explicit-evidence', ...
    'errorIdentifier',identifier,'orthogonalityError',failure.orthogonalityError,'modelRelationRelativeDifference',failure.modelRelationRelativeDifference, ...
    'arnoldiSteps',failure.arnoldiSteps,'CApplications',failure.modelApplications,'ILUApplications',failure.preconditionerApplications, ...
    'singleScanRelativeCResidual',norm(remaining-Csingle,'fro')/norm(remaining,'fro'),'singleDirectionNormRelativeToRhs',norm(single,'fro')/norm(remaining,'fro'), ...
    'newFactorConstructions',0,'newCoarseConstructions',0,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'rejectedDirectionUsable',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
fid=fopen(fullfile(folder,'block_ilu_failure_diagnostic_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
