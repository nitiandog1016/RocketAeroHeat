function report=studyRocketSpatialSchurKrylov(folder)
%STUDYROCKETSPATIALSCHURKRYLOV 只读已验空间消元，区分粗S剩余与完整F行误差。
% 模块一：不重装F逆/近似S，粗S仅作首次32方向试算并保存全部小基。
% 模块二：代回完整原C后分别报告F/C行，S目标通过不能替代原C目标。
% 模块三：这是诊断，不更新原NS状态或把新向量并入接受历史。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SpatialSchurStudy:Overwrite','不得覆盖空间粗S诊断。');mkdir(folder);
assemblyFile=fullfile(root,'outputs','v83_validation','spatialreductionreal01','spatial_reduction_assembly_evidence.mat');
coarseFactorFile=fullfile(root,'outputs','v83_validation','spatialreductionreal01','spatial_coarse_factor_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(assemblyFile),'b04c5c8d604c52c9289ce5d08127270a41c31322eef5ac44cd0fda4ebc0c76ca'));
assert(strcmp(rocketStreamingFileSHA256(coarseFactorFile),'fc7fcdf76174f97829eab02cc7b0405f42cdef327abeac2cb9ea691be076efaf'));
a=load(assemblyFile,'reduction','packedFile','rhsFile','sourceBinding');reduction=a.reduction;packedFile=a.packedFile;rhsFile=a.rhsFile;sourceBinding=a.sourceBinding;clear a;
s=load(coarseFactorFile,'coarseFactor','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));coarseFactor=s.coarseFactor;clear s;
s=load(packedFile,'packed','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));packed=s.packed;clear s;
s=load(rhsFile,'remaining','sourceFile');assert(strcmp(s.sourceFile,sourceBinding.sourceFile));rhs=s.remaining;clear s;
sourceIdentity=struct;for key={'solveRocketCoupledPatchKrylov','solveRocketCoupledPatches','applyRocketPackedGraph','solveRocketRecycledGMRES', ...
        'compactRocketCoupledKrylovReport','solveRocketSpatialFBlocks','applyRocketRectangularCoupling'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
split=reduction.split;bF=rhs(split.fineOriginalCVs,:);bC=rhs(split.coarseOriginalCVs,:);first=solveRocketSpatialFBlocks(reduction.Ffactor,bF);
coarseRhs=bC-applyRocketRectangularCoupling(reduction.E,first);singleCoarse=solveRocketCoupledPatches(coarseFactor,coarseRhs);singleSImage=applyRocketPackedGraph(reduction.schur,singleCoarse);
singleSResidual=norm(coarseRhs-singleSImage,'fro')/norm(coarseRhs,'fro');
save(fullfile(folder,'spatial_schur_condensed_rhs.mat'),'rhs','bF','bC','coarseRhs','first','singleCoarse','singleSImage','singleSResidual','sourceBinding','-v7');
rejected=false;identifier='';message='';complete=false;stage='first-bounded-approximate-S-GMRES';
try
    options=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',256*1024^2,'saveBasisEvidence',true, ...
        'failureEvidenceFile',fullfile(folder,'spatial_schur_arnoldi_failure.mat'));
    [coarseValue,detail]=solveRocketCoupledPatchKrylov(reduction.schur,coarseFactor,coarseRhs,options);Simage=applyRocketPackedGraph(reduction.schur,coarseValue);
    value=zeros(size(rhs));value(split.coarseOriginalCVs,:)=coarseValue;
    value(split.fineOriginalCVs,:)=solveRocketSpatialFBlocks(reduction.Ffactor,bF-applyRocketRectangularCoupling(reduction.B,coarseValue));
    image=applyRocketPackedGraph(packed,value);remaining=rhs-image;denominator=norm(rhs,'fro');
    diagnostic=struct('relativeApproximateSResidual',norm(coarseRhs-Simage,'fro')/norm(coarseRhs,'fro'), ...
        'relativeFullCResidual',norm(remaining,'fro')/denominator,'fineRowsRelativeToOriginalRHS',norm(remaining(split.fineOriginalCVs,:),'fro')/denominator, ...
        'coarseRowsRelativeToOriginalRHS',norm(remaining(split.coarseOriginalCVs,:),'fro')/denominator, ...
        'coarseRowsCondensationConsistency',norm(remaining(split.coarseOriginalCVs,:)-(coarseRhs-Simage),'fro')/max(1,denominator), ...
        'fineRowsRemainingSquaredNormFraction',norm(remaining(split.fineOriginalCVs,:),'fro')^2/norm(remaining,'fro')^2, ...
        'originalCInternalTargetMet',norm(remaining,'fro')/denominator<=1e-3,'approximateSTargetMet',detail.krylov.flag==0);
    save(fullfile(folder,'spatial_schur_krylov_evidence.mat'),'coarseRhs','coarseValue','Simage','detail','rhs','value','image','remaining','diagnostic', ...
        'assemblyFile','coarseFactorFile','packedFile','rhsFile','sourceBinding','sourceIdentity','-v7');
    assert(diagnostic.coarseRowsCondensationConsistency<1e-8,'SpatialSchurStudy:Condensation','代回完整原C的粗行与凝聚剩余不符。');
    complete=true;stage='bounded-approximate-S-and-parent-C-diagnostic-complete';
catch ex
    if ~(startsWith(ex.identifier,'SpatialSchurStudy:')||startsWith(ex.identifier,'SpatialF:')||startsWith(ex.identifier,'SpatialCoupling:')|| ...
            startsWith(ex.identifier,'CoupledPatches:')||startsWith(ex.identifier,'RecycledGMRES:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-first-bounded-spatial-approximate-Schur-solve-and-parent-C-diagnostic','stage',stage,'candidateRejected',rejected,'errorIdentifier',identifier,'message',message, ...
    'numericalRunComplete',complete,'systemScope','approximate-spatial-Schur-S-with-all-five-fields-not-original-full-C', ...
    'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'singleApproximateSResidual',singleSResidual, ...
    'newFInverseConstructions',0,'newSchurAssemblies',0,'newSelectedJApplications',0,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if complete,report.diagnostic=diagnostic;report.innerApproximateSReport=compactRocketCoupledKrylovReport(detail);report.fullSArnoldiBasisSaved=true;report.newOriginalCDiagnosticApplications=1;report.extraSProbeApplications=2;end
raw=jsonencode(report);assert(numel(raw)<128*1024,'SpatialSchurStudy:Summary','粗S诊断摘要过大。');fid=fopen(fullfile(folder,'spatial_schur_krylov_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
