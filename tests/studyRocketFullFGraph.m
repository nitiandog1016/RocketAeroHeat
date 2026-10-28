function report=studyRocketFullFGraph(folder)
%STUDYROCKETFULLFGRAPH 首次对真实完整A_FF作有界同右端比较，不重建原F逆或旧S。
% 模块一：精确提取所有原F-F系数，旧粗S向量只提供固定C值和对照。
% 模块二：完整F迭代期间释放原C内存；数值证据先存，再重读原C作诊断。
% 模块三：F自身目标、原C行残差和旧近似S身份分开，禁止直接更新物理场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FullFStudy:Overwrite','不得覆盖真实完整F研究。');mkdir(folder);
assemblyFile=fullfile(root,'outputs','v83_validation','spatialreductionreal01','spatial_reduction_assembly_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(assemblyFile),'b04c5c8d604c52c9289ce5d08127270a41c31322eef5ac44cd0fda4ebc0c76ca'));
a=load(assemblyFile,'reduction','packedFile','sourceBinding');split=a.reduction.split;Ffactor=a.reduction.Ffactor;B=a.reduction.B;packedFile=a.packedFile;sourceBinding=a.sourceBinding;clear a;
assert(strcmp(rocketStreamingFileSHA256(packedFile),'beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'));
sourceCandidateFile=fullfile(root,'outputs','v83_validation','spatialschurreal01','spatial_schur_krylov_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(sourceCandidateFile),'6a9f2798bc996cbe9e1cfdeab594bb7d6e7f79e53be7d05d188ca2ab90b52b3c'));
s=load(sourceCandidateFile,'value','rhs','diagnostic','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));originalRHS=s.rhs;fixedC=s.value(split.coarseOriginalCVs,:);previousFValue=s.value(split.fineOriginalCVs,:);previousParentDiagnostic=s.diagnostic;clear s;
sourceIdentity=struct;for key={'extractRocketFullFGraph','solveRocketFullFKrylov','compactRocketFullFReport','applyRocketPackedGraph','solveRocketRecycledGMRES', ...
        'solveRocketSpatialFBlocks','applyRocketRectangularCoupling'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
rejected=false;identifier='';message='';stage='exact-all-FF-extraction';graphAvailable=false;inputsAvailable=false;innerComplete=false;parentComplete=false;
try
    s=load(packedFile,'packed');packed=s.packed;clear s;[graph,extractionReport]=extractRocketFullFGraph(packed,split);
    assert(extractionReport.allOriginalFFBlocks==1285801&&extractionReport.crossPatchFFBlocks==682056,'FullFStudy:Source','完整F结构与已验源不符。');graphAvailable=true;
    graphFile=fullfile(folder,'full_F_graph_evidence.mat');save(graphFile,'graph','extractionReport','packedFile','assemblyFile','sourceBinding','sourceIdentity','-v7');
    rhsF=originalRHS(split.fineOriginalCVs,:)-applyRocketRectangularCoupling(B,fixedC);previousFImage=applyRocketPackedGraph(graph,previousFValue);
    previousFullFResidual=norm(rhsF-previousFImage,'fro')/norm(rhsF,'fro');
    save(fullfile(folder,'full_F_inputs.mat'),'rhsF','fixedC','previousFValue','previousFImage','previousFullFResidual','originalRHS', ...
        'previousParentDiagnostic','graphFile','assemblyFile','sourceCandidateFile','sourceBinding','-v7');inputsAvailable=true;clear packed;
    stage='first-bounded-full-FF-GMRES';options=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2, ...
        'saveBasisEvidence',false,'failureEvidenceFile',fullfile(folder,'full_F_arnoldi_failure.mat'));
    [valueF,detail]=solveRocketFullFKrylov(graph,Ffactor,rhsF,options);innerComplete=true;Fimage=applyRocketPackedGraph(graph,valueF);
    save(fullfile(folder,'full_F_inner_evidence.mat'),'rhsF','valueF','Fimage','detail','graphFile','assemblyFile','sourceBinding','sourceIdentity','-v7');
    stage='fixed-C-parent-original-C-diagnostic';s=load(packedFile,'packed');packed=s.packed;clear s;value=zeros(size(originalRHS));
    value(split.coarseOriginalCVs,:)=fixedC;value(split.fineOriginalCVs,:)=valueF;image=applyRocketPackedGraph(packed,value);remaining=originalRHS-image;
    diagnostic=struct('relativeFullFResidual',norm(rhsF-Fimage,'fro')/norm(rhsF,'fro'), ...
        'relativeParentOriginalCResidual',norm(remaining,'fro')/norm(originalRHS,'fro'), ...
        'FRowsRelativeToOriginalRHS',norm(remaining(split.fineOriginalCVs,:),'fro')/norm(originalRHS,'fro'), ...
        'CRowsRelativeToOriginalRHS',norm(remaining(split.coarseOriginalCVs,:),'fro')/norm(originalRHS,'fro'), ...
        'FRowsSubmatrixConsistency',norm(remaining(split.fineOriginalCVs,:)-(rhsF-Fimage),'fro')/max(1,norm(originalRHS,'fro')), ...
        'oldApproximateSchurIdentityClaimedForNewFValue',false,'parentOriginalCTargetMet',norm(remaining,'fro')/norm(originalRHS,'fro')<=1e-3);
    save(fullfile(folder,'full_F_real_evidence.mat'),'rhsF','valueF','Fimage','detail','value','image','remaining','originalRHS','fixedC','diagnostic', ...
        'graphFile','assemblyFile','packedFile','sourceCandidateFile','sourceBinding','sourceIdentity','-v7');
    assert(diagnostic.FRowsSubmatrixConsistency<1e-8,'FullFStudy:ParentRows','原C的F行与完整F子系统不一致。');
    parentComplete=true;stage='bounded-full-FF-and-fixed-C-parent-study-complete';
catch ex
    if ~(startsWith(ex.identifier,'FullF')||startsWith(ex.identifier,'RecycledGMRES:')||startsWith(ex.identifier,'SpatialF:')|| ...
            startsWith(ex.identifier,'SpatialCoupling:')||startsWith(ex.identifier,'BlockGraph:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-first-bounded-all-original-five-field-FF-subsystem-study','stage',stage,'candidateRejected',rejected,'errorIdentifier',identifier,'message',message, ...
    'graphAvailable',graphAvailable,'inputsAvailable',inputsAvailable,'innerRunComplete',innerComplete,'parentDiagnosticComplete',parentComplete, ...
    'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'systemScope','exact-all-five-field-A_FF-with-fixed-C-values-not-original-full-C-solve', ...
    'newLocalFInverseConstructions',0,'newApproximateSchurAssemblies',0,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if graphAvailable,report.extractionReport=extractionReport;end
if inputsAvailable,report.previousFullFResidual=previousFullFResidual;report.previousParentDiagnostic=previousParentDiagnostic;end
if innerComplete,report.innerFullFReport=compactRocketFullFReport(detail);report.fullFInternalTargetMet=detail.krylov.flag==0;report.fullFArnoldiBasisSaved=false;report.extraFullFDiagnosticApplications=2;end
if parentComplete,report.parentDiagnostic=diagnostic;report.newOriginalCDiagnosticApplications=1;end
raw=jsonencode(report);assert(numel(raw)<128*1024,'FullFStudy:Summary','完整F摘要超过128KiB。');fid=fopen(fullfile(folder,'full_F_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
