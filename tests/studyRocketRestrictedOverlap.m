function report=studyRocketRestrictedOverlap(folder)
%STUDYROCKETRESTRICTEDOVERLAP 复用冻结核心/两邻点拓扑，对完整真实F作一次32方向对照。
% 模块一：完整F、右端和固定C值与此前两条真实F试算保持原字节身份。
% 模块二：先验完整数组预算，再认证所有扩展逆，失败记录保留，不重选邻点。
% 模块三：只作F/固定C父矩阵诊断，不改变旧S或原NS状态。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'OverlapFStudy:Overwrite','不得覆盖受限重叠真实研究。');mkdir(folder);
graphFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_graph_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(graphFile),'b8ac74781499c0054ea49af77f3bad41f6ae6600dea9045264bba0230c5e4a04'));
s=load(graphFile,'graph','packedFile','sourceBinding');graph=s.graph;packedFile=s.packedFile;sourceBinding=s.sourceBinding;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
topologyFile=fullfile(root,'outputs','v83_validation','fullFoverlappreflight01','full_F_overlap_input.mat');
assert(strcmp(rocketStreamingFileSHA256(topologyFile),'af0061652a19dddb4f80980ebd7e3e33f07e033d5901239ddb0efdcd93403910'));topology=load(topologyFile);
assert(strcmp(topology.graphEvidenceSHA256,rocketStreamingFileSHA256(graphFile)));
oldFactorFile=fullfile(root,'outputs','v83_validation','fullFstrengthreal02','full_F_strength_factor_evidence.mat');
assert(strcmp(topology.sourceStrengthFactorSHA256,rocketStreamingFileSHA256(oldFactorFile)));
inputFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_inputs.mat');
assert(strcmp(rocketStreamingFileSHA256(inputFile),'8d90e8f6a9416b189ecae56e029a500c889f47d413534d7bb31a73f53d06bbf6'));
s=load(inputFile,'rhsF','originalRHS','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));rhsF=s.rhsF;originalRHS=s.originalRHS;clear s;
baseFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_real_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(baseFile),'1bf2bfef8fd8f5641f0e2441ab358991b553cd3a5ab77c3a96d15a3bf4b576d5'));
s=load(baseFile,'value','diagnostic');baseValue=s.value;oldDiagnostic=s.diagnostic;clear s;
sourceIdentity=struct;for key={'factorRocketRestrictedOverlap','solveRocketRestrictedOverlap','solveRocketOverlapFullF','compactRocketFullFReport', ...
        'certifyRocketLocalInverse','rocketCompensatedInverseResidual','applyRocketPackedGraph','solveRocketRecycledGMRES'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
stage='all-frozen-extended-local-inverses';rejected=false;identifier='';message='';factorAvailable=false;singleComplete=false;innerComplete=false;parentComplete=false;
try
    options=struct('maximumFactorBytes',512*1024^2,'maximumConstructionBytes',1536*1024^2,'failureEvidenceFile',fullfile(folder,'restricted_overlap_inverse_failure.mat'));
    factor=factorRocketRestrictedOverlap(graph,topology,options);
    assert(factor.factorizationReport.originalFFBlockOccurrencesInsideOverlapSets==863948,'OverlapFStudy:Topology','实际扩展原块数与冻结预检不符。');
    measuredFactor=whos('factor');measuredGraph=whos('graph');cacheBytes=measuredFactor.bytes+measuredGraph.bytes;
    assert(cacheBytes<=768*1024^2,'OverlapF:Budget','完整F图与受限重叠缓存共同超768MiB预算。');
    factorFile=fullfile(folder,'restricted_overlap_factor_evidence.mat');save(factorFile,'factor','graphFile','topologyFile','inputFile','baseFile','sourceBinding','sourceIdentity','measuredFactor','measuredGraph','cacheBytes','-v7');factorAvailable=true;
    stage='single-fixed-core-restricted-map';single=solveRocketRestrictedOverlap(factor,rhsF);singleImage=applyRocketPackedGraph(graph,single);singleResidual=norm(rhsF-singleImage,'fro')/norm(rhsF,'fro');singleComplete=true;
    stage='bounded-complete-FF-GMRES';kopt=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2, ...
        'saveBasisEvidence',false,'failureEvidenceFile',fullfile(folder,'restricted_overlap_arnoldi_failure.mat'));
    [valueF,detail]=solveRocketOverlapFullF(graph,factor,rhsF,kopt);Fimage=applyRocketPackedGraph(graph,valueF);innerComplete=true;
    save(fullfile(folder,'restricted_overlap_inner_evidence.mat'),'rhsF','single','singleImage','singleResidual','valueF','Fimage','detail','factorFile','graphFile','sourceBinding','sourceIdentity','-v7');
    stage='fixed-C-parent-original-C-diagnostic';s=load(packedFile,'packed');packed=s.packed;clear s;value=baseValue;value(graph.originalPhysicalCVs,:)=valueF;
    image=applyRocketPackedGraph(packed,value);remaining=originalRHS-image;isF=false(size(value,1),1);isF(graph.originalPhysicalCVs)=true;
    diagnostic=struct('relativeFullFResidual',norm(rhsF-Fimage,'fro')/norm(rhsF,'fro'),'relativeParentOriginalCResidual',norm(remaining,'fro')/norm(originalRHS,'fro'), ...
        'FRowsRelativeToOriginalRHS',norm(remaining(isF,:),'fro')/norm(originalRHS,'fro'),'CRowsRelativeToOriginalRHS',norm(remaining(~isF,:),'fro')/norm(originalRHS,'fro'), ...
        'FRowsSubmatrixConsistency',norm(remaining(graph.originalPhysicalCVs,:)-(rhsF-Fimage),'fro')/max(1,norm(originalRHS,'fro')), ...
        'parentOriginalCTargetMet',norm(remaining,'fro')/norm(originalRHS,'fro')<=1e-3,'oldApproximateSchurIdentityClaimedForNewMap',false);
    save(fullfile(folder,'restricted_overlap_real_evidence.mat'),'rhsF','single','singleImage','singleResidual','valueF','Fimage','detail','value','image','remaining','originalRHS', ...
        'diagnostic','factorFile','graphFile','topologyFile','inputFile','baseFile','packedFile','sourceBinding','sourceIdentity','-v7');
    assert(diagnostic.FRowsSubmatrixConsistency<1e-8,'OverlapFStudy:ParentRows','完整F与原C的F行不符。');parentComplete=true;stage='bounded-restricted-overlap-study-complete';
catch ex
    if ~(startsWith(ex.identifier,'OverlapF')||startsWith(ex.identifier,'LocalInverse:')||startsWith(ex.identifier,'InverseResidual:')|| ...
            startsWith(ex.identifier,'RecycledGMRES:')||startsWith(ex.identifier,'BlockGraph:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-real-frozen-two-neighbor-restricted-core-overlap-F-preconditioner','stage',stage,'candidateRejected',rejected,'errorIdentifier',identifier,'message',message, ...
    'factorAvailable',factorAvailable,'singleComplete',singleComplete,'innerRunComplete',innerComplete,'parentDiagnosticComplete',parentComplete, ...
    'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'systemScope','complete-A_FF-only-with-separate-fixed-C-parent-diagnostic', ...
    'oldThirtyTwoStepFullFResidual',oldDiagnostic.relativeFullFResidual,'newFullFGraphExtractions',0,'newApproximateSchurAssemblies',0, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if factorAvailable,report.factorizationReport=factor.factorizationReport;report.actualFactorArrayBytes=measuredFactor.bytes;report.actualGraphArrayBytes=measuredGraph.bytes;report.actualCombinedCacheBytes=cacheBytes;end
if singleComplete,report.singleFullFResidual=singleResidual;end
if innerComplete,report.innerFullFReport=compactRocketFullFReport(detail);report.fullFInternalTargetMet=detail.krylov.flag==0;report.fullFArnoldiBasisSaved=false;report.extraFullFDiagnosticApplications=2;end
if parentComplete,report.parentDiagnostic=diagnostic;report.newOriginalCDiagnosticApplications=1;end
raw=jsonencode(report);assert(numel(raw)<128*1024,'OverlapFStudy:Summary','受限重叠摘要超过128KiB。');fid=fopen(fullfile(folder,'restricted_overlap_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
