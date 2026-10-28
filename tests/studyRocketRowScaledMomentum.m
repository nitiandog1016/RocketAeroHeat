function report=studyRocketRowScaledMomentum(folder)
%STUDYROCKETROWSCALEDMOMENTUM 复用冻结真实因子，作一次同右端等价行尺度32方向对照。
% 模块一：不重组法向列或原S，明确绑定已独立核对的A、因子和数值右端。
% 模块二：只改变内部代数坐标，原单位A剩余单列，与旧同Q未尺度结果比较。
% 模块三：完整V/Z仅保存在MAT；JSON白名单摘要，不把内部结果作为物理场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RowMomentumStudy:Overwrite','不得覆盖真实行尺度对照。');mkdir(folder);
factorFile=fullfile(root,'outputs','v83_validation','momentumlinereal03','momentum_line_factor_evidence.mat');
factorSHA256=rocketStreamingFileSHA256(factorFile);assert(strcmp(factorSHA256,'a5e0dcfaaa87125e13fd9d5d1a93d14c265c1b5be8f7a06ae320539c30c17cb0'),'RowMomentumStudy:FactorIdentity','因子不是已独立核对的原字节。');
f=load(factorFile,'factor','assemblyFile','sourceBinding');factor=f.factor;assemblyFile=f.assemblyFile;sourceBinding=f.sourceBinding;clear f;
assert(strcmp(rocketStreamingFileSHA256(assemblyFile),'efffcbd85ba71246dac8aacb735b8efd647f23bc2933332ab143c426d82db1e7'));
a=load(assemblyFile,'momentumMatrix','sourceBinding');assert(isequaln(a.sourceBinding,sourceBinding));matrix=a.momentumMatrix;clear a;
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
rhsFile=fullfile(root,'outputs','v83_validation','momentumlinereal03','momentum_line_krylov_evidence.mat');rhsSHA256=rocketStreamingFileSHA256(rhsFile);
assert(strcmp(rhsSHA256,'2538b3c65eed77a21b7483c5f9bf8f07c98e014f026f2f37acee4800ccf0abd3'));s=load(rhsFile,'rhs','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));rhs=s.rhs;clear s;
referenceFile=fullfile(root,'outputs','v83_validation','momentumlinereport02','momentum_lines_real_summary.json');reference=jsondecode(fileread(referenceFile));assert(isequaln(reference.sourceBinding,sourceBinding));
sourceIdentity=struct;for key={'rocketMomentumRowScales','solveRocketRowScaledMomentumKrylov','solveRocketMomentumLines','applyRocketMomentumGraph','solveRocketRecycledGMRES','compactRocketMomentumKrylovReport'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
options=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',true,'failureEvidenceFile',fullfile(folder,'row_momentum_arnoldi_failure.mat'));
rejected=false;identifier='';message='';complete=false;
try
    [value,detail,rowScale,image]=solveRocketRowScaledMomentumKrylov(matrix,factor,rhs,options);complete=true;
    save(fullfile(folder,'row_momentum_real_evidence.mat'),'rhs','value','detail','rowScale','image','factorFile','factorSHA256','assemblyFile','rhsFile','rhsSHA256','sourceBinding','sourceIdentity','referenceFile','-v7');
catch ex
    if ~(startsWith(ex.identifier,'RowMomentum:')||startsWith(ex.identifier,'MomentumScale:')||startsWith(ex.identifier,'MomentumGraph:')||startsWith(ex.identifier,'MomentumLines:')||startsWith(ex.identifier,'RecycledGMRES:')),rethrow(ex);end
    rejected=true;identifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-one-equivalent-row-scale-comparison-on-frozen-original-momentum-A','candidateRejected',rejected,'errorIdentifier',identifier,'message',message, ...
    'numericalRunComplete',complete,'factorFile',factorFile,'factorSHA256',factorSHA256,'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity, ...
    'unscaledReferenceOriginalAResidual',reference.innerReport.relativeAResidual,'unscaledReferenceArnoldiSteps',reference.innerReport.krylov.arnoldiSteps, ...
    'newFactorConstructions',0,'newSchurAssemblies',0,'newSelectedJApplications',0,'newActualNSEvaluations',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0, ...
    'wholeFiveFieldPreconditionerQualified',false,'physicalSolution',false,'engineProcessID',feature('getpid'));
if complete
    report.innerReport=compactRocketMomentumKrylovReport(detail);report.originalAResidualImproved=detail.relativeAResidual<reference.innerReport.relativeAResidual;
    report.originalUnitTargetMet=detail.originalUnitTargetMet;report.rowScaledTargetMet=detail.rowScaledTargetMet;
end
raw=jsonencode(report);assert(numel(raw)<128*1024,'RowMomentumStudy:SummarySize','研究摘要超限。');fid=fopen(fullfile(folder,'row_momentum_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
