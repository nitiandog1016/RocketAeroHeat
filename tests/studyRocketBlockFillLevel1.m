function report=studyRocketBlockFillLevel1(folder)
%STUDYROCKETBLOCKFILLLEVEL1 冻结真实C/排序的一级块ILU及有界内部方向研究。
% 模块一：原C不重构，新符号模式不改变原矩阵；新增填充只服务于预条件。
% 模块二：预算、局部逆及完整Arnoldi门槛保留，拒绝时保留输入及失败证据。
% 模块三：内部模型研究不调用完整J/NS，不增加已解面积或接受端点。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'BlockFillStudy:Overwrite','不得覆盖真实一级填充研究。');mkdir(folder);
packedFile=fullfile(root,'outputs','v83_validation','blockilureal01','steady_graph_packed_evidence.mat');s=load(packedFile,'packed','sourceFile','sourceBinding');
packed=s.packed;sourceFile=s.sourceFile;sourceBinding=s.sourceBinding;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceFile),sourceBinding.sourceFileSHA256));
anchorFile=fullfile(root,'outputs','v83_validation','innergraphreal01','inner_graph_real_evidence.mat');s=load(anchorFile,'remaining','sourceFile');
assert(strcmp(sourceFile,s.sourceFile));rhs=s.remaining;clear s;
stage='enrichment';factorAvailable=false;innerRejected=false;errorIdentifier='';detail=struct;singleReport=struct;singleResidual=NaN;singleAmplification=NaN;
try
    enriched=enrichRocketBlockGraphLevel1(packed);enrichmentReport=enriched.enrichmentReport;
    pattern=struct('rowPointer',enriched.rowPointer,'column',enriched.column,'diagonalPosition',enriched.diagonalPosition, ...
        'order',enriched.order,'inverseOrder',enriched.inverseOrder,'enrichmentReport',enrichmentReport);
    save(fullfile(folder,'expanded_pattern_evidence.mat'),'packedFile','sourceFile','sourceBinding','pattern','-v7');
    probe=reshape(sin((1:numel(rhs))'*.37),size(rhs));reference=applyRocketPackedGraph(packed,probe);actual=applyRocketPackedGraph(enriched,probe);
    actionError=norm(actual-reference,'fro')/max(1,norm(reference,'fro'));assert(actionError<2e-10,'BlockFillStudy:C','新增零块改变了原C作用。');clear packed probe reference actual;
    stage='factorization';enriched=factorRocketBlockILU0(enriched);factor=enriched;clear enriched;factorAvailable=true;
    factorFile=fullfile(folder,'block_fill1_factor_evidence.mat');factorIdentity=struct;
    for key={'enrichRocketBlockGraphLevel1','factorRocketBlockILU0','solveRocketBlockILU0','solveRocketRecycledGMRES'}
        factorIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));
    end
    save(factorFile,'factor','factorIdentity','packedFile','anchorFile','sourceFile','sourceBinding','-v7');
    s=load(packedFile,'packed');packed=s.packed;clear s;
    [single,singleReport]=solveRocketBlockILU0(factor,rhs);singleImage=applyRocketPackedGraph(packed,single);
    singleResidual=norm(rhs-singleImage,'fro')/norm(rhs,'fro');singleAmplification=norm(single,'fro')/norm(rhs,'fro');clear single singleImage;
    options=struct('restartLength',32,'maximumCycles',1,'retainedDirections',0,'maximumArnoldiSteps',32,'absoluteTolerance',1e-3*norm(rhs,'fro'), ...
        'maximumBasisBytes',1536*1024^2,'saveBasisEvidence',true,'robustArnoldi',true,'failureEvidenceFile',fullfile(folder,'scaled_arnoldi_failure_evidence.mat'));
    save(fullfile(folder,'input_binding.mat'),'factorFile','packedFile','anchorFile','sourceFile','sourceBinding','rhs','options','-v7');stage='inner-Arnoldi';
    timer=tic;[x,krylov]=solveRocketRecycledGMRES(@action,@precondition,rhs(:),options);value=reshape(x,size(rhs));clear x;
    detail=struct('method','symbolic-level1-five-conservation-block-ILU-scaled-adaptive-one-cycle-GMRES','krylov',krylov, ...
        'relativeCResidual',krylov.relativeResidual,'CApplications',krylov.modelApplications,'blockILUApplications',krylov.preconditionerApplications, ...
        'elapsedSeconds',toc(timer),'outerFixedLinearMapClaimed',false,'requiresFlexibleOuterMethod',true,'originalNSAcceptanceClaimed',false);
    stage='bounded-inner-complete';
catch ex
    if ~(startsWith(ex.identifier,'BlockFill:')||startsWith(ex.identifier,'BlockILU:')||strcmp(ex.identifier,'RecycledGMRES:Arnoldi')),rethrow(ex);end
    errorIdentifier=ex.identifier;innerRejected=true;
end
report=struct('status','tested-real-level1-block-pattern-research-with-explicit-rejection-status', ...
    'stage',stage,'factorAvailable',factorAvailable,'candidateRejected',innerRejected,'errorIdentifier',errorIdentifier, ...
    'singleRelativeCResidual',singleResidual,'singleDirectionNormRelativeToRhs',singleAmplification, ...
    'newFactorConstructionAttempts',double(~strcmp(stage,'enrichment')),'newCoarseConstructions',0,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
if exist('enrichmentReport','var'),report.enrichmentReport=rmfield(enrichmentReport,'originalToEnrichedPosition');end
if factorAvailable,report.factorizationReport=factor.factorizationReport;end
if ~innerRejected
    report.relativeCResidual=detail.relativeCResidual;report.innerTargetMet=krylov.flag==0;
    save(fullfile(folder,'scaled_ilu_replay_evidence.mat'),'report','packedFile','factorFile','anchorFile','sourceFile','sourceBinding','rhs','options','value','detail','-v7');
elseif strcmp(stage,'inner-Arnoldi')
    f=load(options.failureEvidenceFile,'failure');report.failureOrthogonality=f.failure.orthogonalityError;report.failureModelRelation=f.failure.modelRelationRelativeDifference;
    report.arnoldiSteps=f.failure.arnoldiSteps;report.CApplications=f.failure.modelApplications;report.ILUApplications=f.failure.preconditionerApplications;
end
fid=fopen(fullfile(folder,'block_fill1_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function out=action(v)
        product=applyRocketPackedGraph(packed,reshape(v,size(rhs)));out=product(:);
    end
    function out=precondition(v)
        product=solveRocketBlockILU0(factor,reshape(v,size(rhs)));out=product(:);
    end
end
