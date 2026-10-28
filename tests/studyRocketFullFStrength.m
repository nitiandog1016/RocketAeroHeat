function report=studyRocketFullFStrength(folder)
%STUDYROCKETFULLFSTRENGTH 仅在已验完整F图内重新强连接分组，比较同右端32步。
% 模块一：不再继承原C点移除后残缺的F组，但不改变任何物理CV或F/C角色。
% 模块二：新局部F缓存独立认证，完整A_FF全部连接保持，旧S不冒称与新映射一致。
% 模块三：固定原C值回代仅作完整原矩阵诊断，不更新原NS流场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FullFStrengthStudy:Overwrite','不得覆盖F内强连接真实试算。');mkdir(folder);
graphFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_graph_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(graphFile),'b8ac74781499c0054ea49af77f3bad41f6ae6600dea9045264bba0230c5e4a04'));
s=load(graphFile,'graph','packedFile','sourceBinding');graph=s.graph;packedFile=s.packedFile;sourceBinding=s.sourceBinding;clear s;
inputFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_inputs.mat');
assert(strcmp(rocketStreamingFileSHA256(inputFile),'8d90e8f6a9416b189ecae56e029a500c889f47d413534d7bb31a73f53d06bbf6'));
s=load(inputFile,'rhsF','originalRHS','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));rhsF=s.rhsF;originalRHS=s.originalRHS;clear s;
baseFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_real_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(baseFile),'1bf2bfef8fd8f5641f0e2441ab358991b553cd3a5ab77c3a96d15a3bf4b576d5'));
s=load(baseFile,'value','diagnostic');baseValue=s.value;oldDiagnostic=s.diagnostic;clear s;
sourceIdentity=struct;for key={'buildRocketStrongPatches','factorRocketCoupledPatches','certifyRocketLocalInverse','rocketCompensatedInverseResidual','rocketFullFStrengthCoverage', ...
        'solveRocketCoupledPatches','solveRocketCoupledPatchKrylov','compactRocketCoupledKrylovReport','applyRocketPackedGraph','solveRocketRecycledGMRES'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
rejected=false;identifier='';message='';factorAvailable=false;singleComplete=false;innerComplete=false;parentComplete=false;stage='F-only-strength-local-factor';
try
    partition=buildRocketStrongPatches(graph,8);factor=factorRocketCoupledPatches(graph,partition, ...
        struct('maximumFactorBytes',384*1024^2,'maximumConstructionBytes',1536*1024^2,'failureEvidenceFile',fullfile(folder,'F_strength_inverse_failure.mat')));
    cache=whos('graph','factor');assert(sum([cache.bytes])<=768*1024^2,'FullFStrengthStudy:Budget','完整F图与新局部缓存共同超预算。');
    old=load(fullfile(root,'outputs','v83_validation','spatialreductionpreflight01','spatial_reduction_input.mat'),'originalGroupIDOneBased');oldGroup=old.originalGroupIDOneBased(double(graph.originalPhysicalCVs));clear old;
    coverage=rocketFullFStrengthCoverage(graph,factor,oldGroup);clear oldGroup;
    factorFile=fullfile(folder,'full_F_strength_factor_evidence.mat');save(factorFile,'factor','coverage','graphFile','inputFile','baseFile','sourceBinding','sourceIdentity','-v7');factorAvailable=true;
    stage='single-new-F-local-map';single=solveRocketCoupledPatches(factor,rhsF);singleImage=applyRocketPackedGraph(graph,single);singleResidual=norm(rhsF-singleImage,'fro')/norm(rhsF,'fro');singleComplete=true;
    stage='bounded-full-FF-GMRES-with-new-fixed-strength-map';options=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2, ...
        'saveBasisEvidence',false,'failureEvidenceFile',fullfile(folder,'F_strength_arnoldi_failure.mat'));
    [valueF,detail]=solveRocketCoupledPatchKrylov(graph,factor,rhsF,options);Fimage=applyRocketPackedGraph(graph,valueF);innerComplete=true;
    save(fullfile(folder,'full_F_strength_inner_evidence.mat'),'rhsF','valueF','Fimage','detail','single','singleImage','singleResidual','factorFile','graphFile','sourceBinding','sourceIdentity','-v7');
    stage='fixed-C-original-parent-diagnostic';s=load(packedFile,'packed');packed=s.packed;clear s;value=baseValue;value(graph.originalPhysicalCVs,:)=valueF;
    image=applyRocketPackedGraph(packed,value);remaining=originalRHS-image;isF=false(size(value,1),1);isF(graph.originalPhysicalCVs)=true;
    diagnostic=struct('relativeFullFResidual',norm(rhsF-Fimage,'fro')/norm(rhsF,'fro'),'relativeParentOriginalCResidual',norm(remaining,'fro')/norm(originalRHS,'fro'), ...
        'FRowsRelativeToOriginalRHS',norm(remaining(isF,:),'fro')/norm(originalRHS,'fro'),'CRowsRelativeToOriginalRHS',norm(remaining(~isF,:),'fro')/norm(originalRHS,'fro'), ...
        'FRowsSubmatrixConsistency',norm(remaining(graph.originalPhysicalCVs,:)-(rhsF-Fimage),'fro')/max(1,norm(originalRHS,'fro')), ...
        'parentOriginalCTargetMet',norm(remaining,'fro')/norm(originalRHS,'fro')<=1e-3,'oldApproximateSchurIdentityClaimedForNewMap',false);
    save(fullfile(folder,'full_F_strength_real_evidence.mat'),'rhsF','valueF','Fimage','detail','single','singleImage','singleResidual','value','image','remaining','originalRHS', ...
        'diagnostic','factorFile','graphFile','inputFile','baseFile','packedFile','sourceBinding','sourceIdentity','-v7');
    assert(diagnostic.FRowsSubmatrixConsistency<1e-8,'FullFStrengthStudy:ParentRows','新F值与完整原C的F行不一致。');parentComplete=true;stage='bounded-F-strength-study-complete';
catch ex
    if ~(startsWith(ex.identifier,'FullFStrengthStudy:')||startsWith(ex.identifier,'StrongPatches:')||startsWith(ex.identifier,'CoupledPatches:')|| ...
            startsWith(ex.identifier,'LocalInverse:')||startsWith(ex.identifier,'InverseResidual:')||startsWith(ex.identifier,'RecycledGMRES:')||startsWith(ex.identifier,'BlockGraph:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-real-F-only-strength-repartition-with-fixed-C-parent-diagnostic','stage',stage,'candidateRejected',rejected,'errorIdentifier',identifier,'message',message, ...
    'factorAvailable',factorAvailable,'singleComplete',singleComplete,'innerRunComplete',innerComplete,'parentDiagnosticComplete',parentComplete, ...
    'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'systemScope','same-complete-spatial-A_FF-only-not-original-full-C-solve', ...
    'oldThirtyTwoStepFullFResidual',oldDiagnostic.relativeFullFResidual,'newFullFGraphExtractions',0,'newApproximateSchurAssemblies',0,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if factorAvailable,report.factorizationReport=factor.factorizationReport;report.strengthCoverage=coverage;end
if singleComplete,report.singleFullFResidual=singleResidual;end
if innerComplete,report.innerFullFReport=compactRocketCoupledKrylovReport(detail);report.fullFInternalTargetMet=detail.krylov.flag==0;report.fullFArnoldiBasisSaved=false;report.extraFullFDiagnosticApplications=2;end
if parentComplete,report.parentDiagnostic=diagnostic;report.newOriginalCDiagnosticApplications=1;end
raw=jsonencode(report);assert(numel(raw)<128*1024,'FullFStrengthStudy:Summary','F强连接摘要超过128KiB。');fid=fopen(fullfile(folder,'full_F_strength_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
