function report=studyRocketCoupledPatches(folder)
%STUDYROCKETCOUPLEDPATCHES 在同Q完整五量C上作一次有界强连接局部块对照。
% 模块一：只读原C和同Q剩余，不构造原NS/J像或重复旧粗方向。
% 模块二：8CV局部组源于512MiB因子预算，保留全部五量和原C连接。
% 模块三：真实默认不复制完整V/Z到报告，原生Arnoldi护栏保持，最终C作用单独保存。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CoupledPatchStudy:Overwrite','不得覆盖完整五量试算。');mkdir(folder);
packedFile=fullfile(root,'outputs','v83_validation','blockilureal01','steady_graph_packed_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(packedFile),'beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'));
s=load(packedFile,'packed','sourceBinding');packed=s.packed;sourceBinding=s.sourceBinding;clear s;assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
rhsFile=fullfile(root,'outputs','v83_validation','innergraphreal01','inner_graph_real_evidence.mat');s=load(rhsFile,'remaining','sourceFile');assert(strcmp(s.sourceFile,sourceBinding.sourceFile));rhs=s.remaining;clear s;
sourceIdentity=struct;for key={'buildRocketStrongPatches','certifyRocketLocalInverse','factorRocketCoupledPatches','solveRocketCoupledPatches','solveRocketCoupledPatchKrylov','applyRocketPackedGraph','solveRocketRecycledGMRES','rocketCompensatedInverseResidual','compactRocketCoupledKrylovReport'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
stage='bounded-strength-partition';rejected=false;identifier='';message='';factorAvailable=false;innerComplete=false;
try
    partition=buildRocketStrongPatches(packed,8);stage='local-five-field-inverse-factorization';
    options=struct('maximumFactorBytes',512*1024^2,'maximumConstructionBytes',2048*1024^2,'failureEvidenceFile',fullfile(folder,'unsafe_local_inverse.mat'));
    factor=factorRocketCoupledPatches(packed,partition,options);factorAvailable=true;
    factorFile=fullfile(folder,'coupled_patch_factor_evidence.mat');save(factorFile,'factor','packedFile','rhsFile','sourceBinding','sourceIdentity','-v7');
    stage='one-fixed-full-C-action';single=solveRocketCoupledPatches(factor,rhs);singleImage=applyRocketPackedGraph(packed,single);
    singleSummary=struct('relativeCResidual',norm(rhs-singleImage,'fro')/norm(rhs,'fro'),'inverseGainSeconds',norm(single,'fro')/norm(rhs,'fro'));
    stage='bounded-full-five-field-GMRES';krylovOptions=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false,'failureEvidenceFile',fullfile(folder,'coupled_patch_arnoldi_failure.mat'));
    [value,detail]=solveRocketCoupledPatchKrylov(packed,factor,rhs,krylovOptions);innerComplete=true;image=applyRocketPackedGraph(packed,value);
    save(fullfile(folder,'coupled_patch_real_evidence.mat'),'rhs','single','singleImage','singleSummary','value','image','detail','factorFile','packedFile','rhsFile','sourceBinding','sourceIdentity','-v7');stage='bounded-full-five-field-study-complete';
catch ex
    if ~(startsWith(ex.identifier,'StrongPatches:')||startsWith(ex.identifier,'LocalInverse:')||startsWith(ex.identifier,'CoupledPatches:')||startsWith(ex.identifier,'InverseResidual:')||startsWith(ex.identifier,'RecycledGMRES:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-bounded-strong-connection-all-five-field-preconditioner-study','stage',stage,'candidateRejected',rejected,'errorIdentifier',identifier,'message',message, ...
    'factorAvailable',factorAvailable,'innerRunComplete',innerComplete,'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'maximumPatchCells',8, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'newCoarseImageConstructions',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
if factorAvailable,report.partitionReport=partition.constructionReport;report.factorizationReport=factor.factorizationReport;end
if innerComplete,report.singleSummary=singleSummary;report.innerReport=compactRocketCoupledKrylovReport(detail);report.internalCTargetMet=detail.krylov.flag==0;report.extraFullCDiagnosticApplications=2;report.realFullArnoldiBasisSaved=false;end
raw=jsonencode(report);assert(numel(raw)<128*1024,'CoupledPatchStudy:Summary','摘要过大。');fid=fopen(fullfile(folder,'coupled_patch_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
