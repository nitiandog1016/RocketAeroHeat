function report=studyRocketSpatialReduction(folder)
%STUDYROCKETSPATIALREDUCTION 固定真实C/F代表，完整装配有预算的五量局部消元并对照原C。
% 模块一：全部原C、预检代表和结构按SHA绑定，不重选点或重做旧粗空间。
% 模块二：先保存F逆/近似S，再检查三条完整S作用；粗局部逆另存，失败不重装S。
% 模块三：完整C仅作一次32方向试算，原NS/接受/区域资格保持独立。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SpatialReductionStudy:Overwrite','不得覆盖真实空间消元。');mkdir(folder);
inputFile=fullfile(root,'outputs','v83_validation','spatialreductionpreflight01','spatial_reduction_input.mat');
assert(strcmp(rocketStreamingFileSHA256(inputFile),'f7a581f6467423799b0d771923718d56a678040b01b65e4afb9ce4aa83e5df20'));t=load(inputFile);
factorFile=fullfile(root,'outputs','v83_validation','coupledpatchreal01','coupled_patch_factor_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(factorFile),'09b04be0c77f76f8c04bd717745a0dd2fbcef7cdb6809d53bc403a07ed6349e0'));
f=load(factorFile,'factor','packedFile','rhsFile','sourceBinding');partition=f.factor.partition;packedFile=f.packedFile;rhsFile=f.rhsFile;sourceBinding=f.sourceBinding;clear f;
assert(strcmp(rocketStreamingFileSHA256(packedFile),'beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'));
s=load(packedFile,'packed','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));packed=s.packed;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256)&&strcmp(t.sourceEvidenceSHA256,sourceBinding.sourceFileSHA256));
assert(isequal(double(partition.groupID),double(t.originalGroupIDOneBased(:))));
s=load(rhsFile,'remaining','sourceFile');assert(strcmp(s.sourceFile,sourceBinding.sourceFile));rhs=s.remaining;clear s;
sourceIdentity=struct;for key={'buildRocketSpatialSplit','factorRocketSpatialFBlocks','solveRocketSpatialFBlocks','applyRocketRectangularCoupling', ...
        'buildRocketSpatialReduction','applyRocketSpatialSchur','applyRocketSpatialReduction','solveRocketSpatialReductionKrylov','compactRocketSpatialReductionReport', ...
        'certifyRocketLocalInverse','rocketCompensatedInverseResidual','buildRocketStrongPatches','factorRocketCoupledPatches','solveRocketCoupledPatches', ...
        'applyRocketPackedGraph','solveRocketRecycledGMRES'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
stage='all-F-inverses-and-complete-approximate-Schur';rejected=false;identifier='';message='';assemblyComplete=false;schurProbesComplete=false;coarseAvailable=false;singleComplete=false;innerComplete=false;
try
    options=struct('maximumReductionBytes',768*1024^2,'maximumFFFactorBytes',384*1024^2,'maximumConstructionBytes',3072*1024^2, ...
        'maximumSchurBlocks',1000000,'failureEvidenceFile',fullfile(folder,'spatial_ff_failure.mat'));
    reduction=buildRocketSpatialReduction(packed,partition,t.coarseOriginalCVsOneBased,options);assemblyComplete=true;
    assemblyFile=fullfile(folder,'spatial_reduction_assembly_evidence.mat');save(assemblyFile,'reduction','inputFile','factorFile','packedFile','rhsFile','sourceBinding','sourceIdentity','-v7');
    assert(isequal(double(reduction.schur.rowPointer),double(t.schurStructuralRowPointerZeroBased(:))+1)&& ...
        isequal(double(reduction.schur.column),double(t.schurStructuralColumnZeroBased(:))+1),'SpatialReductionStudy:PreflightMismatch','已冻结的完整Schur结构与实际装配不符。');
    stage='three-full-approximate-Schur-probes';schurProbes=cell(3,1);nc=reduction.split.coarseCells;
    for j=1:3
        if j==1,b=rhs(reduction.split.coarseOriginalCVs,:);else,b=reshape(sin((1:5*nc)'*(.217+.13*j)),nc,5);end
        implicit=applyRocketSpatialSchur(reduction,b);assembled=applyRocketPackedGraph(reduction.schur,b);difference=norm(implicit-assembled,'fro')/max(1,norm(implicit,'fro'));
        schurProbes{j}=struct('rhs',b,'implicit',implicit,'assembled',assembled,'relativeDifference',difference);
        assert(difference<1e-8,'SpatialReductionStudy:SchurAction','近似Schur装配与实际固定F映射不符。');
    end
    schurProbesComplete=true;save(fullfile(folder,'spatial_schur_probe_evidence.mat'),'schurProbes','assemblyFile','sourceBinding','-v7');
    stage='fixed-coarse-local-inverses';coarsePartition=buildRocketStrongPatches(reduction.schur,8);
    coarseFactor=factorRocketCoupledPatches(reduction.schur,coarsePartition,struct('maximumFactorBytes',128*1024^2, ...
        'maximumConstructionBytes',3072*1024^2,'failureEvidenceFile',fullfile(folder,'spatial_coarse_failure.mat')));coarseAvailable=true;
    combined=whos('reduction','coarseFactor');assert(sum([combined.bytes])<=896*1024^2,'SpatialReduction:Budget','完整消元及粗局部因子共同超预算。');
    coarseFactorFile=fullfile(folder,'spatial_coarse_factor_evidence.mat');save(coarseFactorFile,'coarseFactor','assemblyFile','sourceBinding','sourceIdentity','-v7');
    stage='single-complete-spatial-back-substitution';[single,singleWork]=applyRocketSpatialReduction(reduction,coarseFactor,rhs);singleImage=applyRocketPackedGraph(packed,single);singleComplete=true;
    singleSummary=struct('relativeCResidual',norm(rhs-singleImage,'fro')/norm(rhs,'fro'),'inverseGainSeconds',norm(single,'fro')/norm(rhs,'fro'));
    stage='bounded-original-C-GMRES';kopt=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2, ...
        'saveBasisEvidence',false,'failureEvidenceFile',fullfile(folder,'spatial_arnoldi_failure.mat'));
    [value,detail]=solveRocketSpatialReductionKrylov(packed,reduction,coarseFactor,rhs,kopt);innerComplete=true;image=applyRocketPackedGraph(packed,value);
    save(fullfile(folder,'spatial_reduction_real_evidence.mat'),'rhs','single','singleImage','singleWork','singleSummary','value','image','detail', ...
        'assemblyFile','coarseFactorFile','inputFile','packedFile','rhsFile','sourceBinding','sourceIdentity','-v7');stage='bounded-original-C-study-complete';
catch ex
    if ~(startsWith(ex.identifier,'Spatial')||startsWith(ex.identifier,'StrongPatches:')||startsWith(ex.identifier,'CoupledPatches:')|| ...
            startsWith(ex.identifier,'LocalInverse:')||startsWith(ex.identifier,'InverseResidual:')||startsWith(ex.identifier,'RecycledGMRES:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-bounded-real-five-field-spatial-local-F-Schur-reduction','stage',stage,'candidateRejected',rejected,'errorIdentifier',identifier,'message',message, ...
    'assemblyComplete',assemblyComplete,'schurProbesComplete',schurProbesComplete,'coarseFactorAvailable',coarseAvailable,'singleComplete',singleComplete, ...
    'innerRunComplete',innerComplete,'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'newRootMatrixAssemblies',0,'newSelectedJApplications',0, ...
    'newActualNSEvaluations',0,'newCoarseImageConstructions',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if assemblyComplete,report.assemblyReport=reduction.constructionReport;report.FFactorReport=reduction.Ffactor.factorizationReport;end
if schurProbesComplete,report.maximumSchurProbeRelativeDifference=max(cellfun(@(v)v.relativeDifference,schurProbes));end
if coarseAvailable,report.coarseFactorReport=coarseFactor.factorizationReport;end
if singleComplete,report.singleSummary=singleSummary;report.singleWork=singleWork;end
if innerComplete,report.innerReport=compactRocketSpatialReductionReport(detail);report.internalCTargetMet=detail.krylov.flag==0;report.extraRootCDiagnosticApplications=2;report.realFullArnoldiBasisSaved=false;end
raw=jsonencode(report);assert(numel(raw)<128*1024,'SpatialReductionStudy:Summary','空间研究摘要超过128KiB。');fid=fopen(fullfile(folder,'spatial_reduction_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
