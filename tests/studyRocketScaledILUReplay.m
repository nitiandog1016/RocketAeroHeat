function report=studyRocketScaledILUReplay(folder)
%STUDYROCKETSCALEDILUREPLAY 冻结原因子与同Q剩余，仅复测缩放复正交的32方向。
% 模块一：绑定原ILU0、原NS端点及已验证当前模型剩余，不重构因子或粗像。
% 模块二：只更换显式内部Arnoldi数值模式，实际C剩余与原1e-8关系均保留。
% 模块三：不调用新的完整J/NS；内部失败和未达目标均不得变成物理解。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'ScaledILUReplay:Overwrite','不得覆盖真实复测。');mkdir(folder);
factorFile=fullfile(root,'outputs','v83_validation','blockilureal01','block_ilu_factor_evidence.mat');s=load(factorFile,'factor','Q','sourceFile','sourceBinding','sourceIdentity');
factor=s.factor;Q=s.Q;sourceFile=s.sourceFile;sourceBinding=s.sourceBinding;sourceIdentity=s.sourceIdentity;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceFile),sourceBinding.sourceFileSHA256));
for key={'packRocketSteadyGraph','factorRocketBlockILU0','solveRocketBlockILU0'}
    assert(strcmp(rocketStreamingFileSHA256(which(key{1})),sourceIdentity.(key{1})),'ScaledILUReplay:Source','原因子来源已变化。');
end
anchorFile=fullfile(root,'outputs','v83_validation','innergraphreal01','inner_graph_real_evidence.mat');s=load(anchorFile,'Q','remaining','sourceFile');
assert(isequaln(Q,s.Q)&&strcmp(sourceFile,s.sourceFile));rhs=s.remaining;clear s;
% 使用已经逐块独立核对的冻结C，避免再载入1.2GB压缩的完整主链结构。
packedFile=fullfile(root,'outputs','v83_validation','blockilureal01','steady_graph_packed_evidence.mat');
s=load(packedFile,'packed','sourceBinding');assert(isequaln(sourceBinding,s.sourceBinding));packed=s.packed;clear s;
for key={'rowPointer','column','diagonalPosition','order','inverseOrder'}
    assert(isequaln(packed.(key{1}),factor.(key{1})),'ScaledILUReplay:Pattern','冻结C与因子排序/模式不符。');
end
options=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1536*1024^2,'saveBasisEvidence',true, ...
    'robustArnoldi',true,'failureEvidenceFile',fullfile(folder,'scaled_arnoldi_failure_evidence.mat'));
save(fullfile(folder,'input_binding.mat'),'factorFile','packedFile','anchorFile','sourceFile','sourceBinding','rhs','options','-v7');
configuration=struct('restartLength',32,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',32, ...
    'absoluteTolerance',1e-3*norm(rhs,'fro'),'maximumBasisBytes',1536*1024^2,'saveBasisEvidence',true, ...
    'robustArnoldi',true,'failureEvidenceFile',options.failureEvidenceFile);
timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),configuration);value=reshape(x,size(rhs));clear x;
detail=struct('method','frozen-complete-C-block-ILU0-scaled-adaptive-one-cycle-GMRES', ...
    'relativeCResidual',krylov.relativeResidual,'krylov',krylov,'CApplications',krylov.modelApplications, ...
    'blockILUApplications',krylov.preconditionerApplications,'elapsedSeconds',toc(timer), ...
    'outerFixedLinearMapClaimed',false,'requiresFlexibleOuterMethod',true,'originalNSAcceptanceClaimed',false);
report=struct('status','tested-real-frozen-ILU0-explicit-scaled-adaptive-Arnoldi-no-physical-update', ...
    'relativeCResidual',detail.relativeCResidual,'innerTargetMet',detail.krylov.flag==0, ...
    'orthogonalityError',detail.krylov.maximumBasisOrthogonalityError,'CApplications',detail.CApplications,'ILUApplications',detail.blockILUApplications, ...
    'newFactorConstructions',0,'newCoarseConstructions',0,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'scaled_ilu_replay_evidence.mat'),'report','factorFile','packedFile','anchorFile','sourceFile','sourceBinding','rhs','options','value','detail','-v7');
fid=fopen(fullfile(folder,'scaled_ilu_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function out=action(v)
        product=applyRocketPackedGraph(packed,reshape(v,size(rhs)));out=product(:);
    end
    function out=precondition(v)
        product=solveRocketBlockILU0(factor,reshape(v,size(rhs)));out=product(:);
    end
end
