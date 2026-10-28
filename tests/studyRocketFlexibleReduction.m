function report=studyRocketFlexibleReduction(folder)
%STUDYROCKETFLEXIBLEREDUCTION 复用冻结全部原C/F对象，作8方向完整父C柔性校正。
% 模块一：原图/代表/局部逆/旧S/像列尺度按原文件身份复用，无新装配或物理删点。
% 模块二：每次原预条件向量和真实外层V/Z保存，供独立全方向及工作量核对。
% 模块三：即使代数目标通过，仍不代表选中J/原NS或区域气动热验收。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
baseFile=fullfile(root,'outputs','v83_validation','spatialreductionreal01','spatial_reduction_real_evidence.mat');
s=load(baseFile,'assemblyFile','coarseFactorFile','packedFile','sourceBinding','rhs','detail');
assemblyFile=s.assemblyFile;coarseFactorFile=s.coarseFactorFile;packedFile=s.packedFile;sourceBinding=s.sourceBinding;rhs=s.rhs;previousFixedSpatial32StepCResidual=s.detail.relativeCResidual;clear s;
assert(strcmp(rocketStreamingFileSHA256(assemblyFile),'b04c5c8d604c52c9289ce5d08127270a41c31322eef5ac44cd0fda4ebc0c76ca'));
assert(strcmp(rocketStreamingFileSHA256(coarseFactorFile),'fc7fcdf76174f97829eab02cc7b0405f42cdef327abeac2cb9ea691be076efaf'));
assert(strcmp(rocketStreamingFileSHA256(packedFile),'beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'));
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
graphFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_graph_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(graphFile),'b8ac74781499c0054ea49af77f3bad41f6ae6600dea9045264bba0230c5e4a04'));
scaleFile=fullfile(root,'outputs','v83_validation','fullFimagescalereal01','full_F_image_scale_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(scaleFile),'99b4e6d25b1255e29ad0d48e379dca2fcae2f442ac50144d144baf8a73b5128c'));
s=load(assemblyFile,'reduction');reduction=s.reduction;clear s;s=load(coarseFactorFile,'coarseFactor');coarseFactor=s.coarseFactor;clear s;
s=load(packedFile,'packed');packed=s.packed;clear s;s=load(graphFile,'graph');Fgraph=s.graph;clear s;s=load(scaleFile,'columnScale');columnScale=s.columnScale;clear s;
sourceIdentity=struct;for key={'applyRocketFlexibleReduction','solveRocketFlexibleReductionKrylov','compactRocketFlexibleReductionReport','solveRocketFullFLeastSquares', ...
        'applyRocketPackedGraph','applyRocketPackedGraphTranspose','solveRocketSpatialFBlocks','solveRocketSpatialFBlocksTranspose','applyRocketRectangularCoupling','solveRocketCoupledPatches','solveRocketRecycledGMRES'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
rejected=false;identifier='';message='';innerComplete=false;stage='input-cache-and-bounded-complete-C-FGMRES';
measured=whos('packed','reduction','coarseFactor','Fgraph','columnScale');cacheBytes=sum([measured.bytes]);
try
    assert(cacheBytes<=1280*1024^2,'FlexibleReduction:Budget','原C/消元/F图/粗对象/列尺度共同缓存超1280MiB。');
    options=struct('maximumIterations',8,'relativeTolerance',1e-3,'innerIterationSequence',32,'innerRelativeTolerance',1e-3, ...
        'maximumInnerOperatorApplications',2048,'maximumBasisBytes',768*1024^2,'saveBasisEvidence',true, ...
        'preconditionerEvidenceFolder',fullfile(folder,'preconditioner'),'failureEvidenceFile',fullfile(folder,'flexible_arnoldi_failure.mat'));
    [value,detail]=solveRocketFlexibleReductionKrylov(packed,reduction,coarseFactor,Fgraph,columnScale,rhs,options);
    image=applyRocketPackedGraph(packed,value);remaining=rhs-image;
    assert(abs(norm(remaining,'fro')/norm(rhs,'fro')-detail.relativeCResidual)<1e-8,'FlexibleReduction:FinalImage','独立最后完整C像与求解摘要不一致。');
    split=reduction.split;diagnostic=struct('FRowsRelativeToOriginalRHS',norm(remaining(split.fineOriginalCVs,:),'fro')/norm(rhs,'fro'), ...
        'CRowsRelativeToOriginalRHS',norm(remaining(split.coarseOriginalCVs,:),'fro')/norm(rhs,'fro'));
    save(fullfile(folder,'flexible_reduction_real_evidence.mat'),'rhs','value','image','remaining','detail','diagnostic','sourceBinding','sourceIdentity', ...
        'baseFile','assemblyFile','coarseFactorFile','packedFile','graphFile','scaleFile','measured','cacheBytes','previousFixedSpatial32StepCResidual','-v7');
    innerComplete=true;stage='bounded-complete-parent-C-flexible-study-complete';
catch ex
    if ~(startsWith(ex.identifier,'FlexibleReduction:')||startsWith(ex.identifier,'RecycledGMRES:')||startsWith(ex.identifier,'FullFLSQR:')|| ...
            startsWith(ex.identifier,'SpatialF:')||startsWith(ex.identifier,'CoupledPatches:')||startsWith(ex.identifier,'BlockGraph:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-bounded-complete-parent-C-flexible-right-preconditioning','stage',stage,'candidateRejected',rejected, ...
    'errorIdentifier',identifier,'message',message,'innerRunComplete',innerComplete,'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity, ...
    'actualInputCacheArrayBytes',cacheBytes,'processPeakMemoryClaimed',false,'previousFixedSpatial32StepCResidual',previousFixedSpatial32StepCResidual, ...
    'newRootMatrixAssemblies',0,'newFGraphExtractions',0,'newLocalInverseConstructions',0,'newApproximateSchurAssemblies',0, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if innerComplete
    report.innerReport=compactRocketFlexibleReductionReport(detail);report.internalCTargetMet=detail.internalCTargetMet;report.finalCDiagnosticApplications=1;
    report.realFullOuterArnoldiBasisSaved=true;report.realInnerLSQRBasisSaved=false;report.finalRowDiagnostic=diagnostic;
end
raw=jsonencode(report);assert(numel(raw)<128*1024,'FlexibleReduction:Summary','真实完整父C可变研究摘要超128KiB。');
fid=fopen(fullfile(folder,'flexible_reduction_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
