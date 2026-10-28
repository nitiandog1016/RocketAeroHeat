function report=studyRocketParentLeastSquares(folder)
%STUDYROCKETPARENTLEASTSQUARES 直接完整父C的像列归一/伴随/有界LSQR真实同端点试算。
% 模块一：复用全部原CV五量、原矩阵与原认证局部块，不使用旧Schur或固定C代表。
% 模块二：尺度和原生坐标单独保存；完整C总残差、原生停止与原NS/物理资格分开。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
factorFile=fullfile(root,'outputs','v83_validation','coupledpatchreal01','coupled_patch_factor_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(factorFile),'09b04be0c77f76f8c04bd717745a0dd2fbcef7cdb6809d53bc403a07ed6349e0'));
s=load(factorFile,'factor','packedFile','sourceBinding');factor=s.factor;packedFile=s.packedFile;sourceBinding=s.sourceBinding;clear s;
assert(strcmp(rocketStreamingFileSHA256(packedFile),'beaa643c209d36ead3d7b859f5968f34c06435023cab2b1e676ac080c9bd0dcf'));
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
s=load(packedFile,'packed');packed=s.packed;clear s;assert(packed.cellCount==217272&&numel(packed.column)==1619664);
baseFile=fullfile(root,'outputs','v83_validation','spatialreductionreal01','spatial_reduction_real_evidence.mat');
s=load(baseFile,'rhs','value','detail','sourceBinding');assert(isequaln(s.sourceBinding,sourceBinding));rhs=s.rhs;initial=s.value;baselineRelativeCResidual=s.detail.relativeCResidual;clear s;
sourceIdentity=struct;for key={'buildRocketParentImageScale','solveRocketParentLeastSquares','solveRocketCoupledPatchesTranspose','solveRocketCoupledPatches','applyRocketPackedGraph','applyRocketPackedGraphTranspose'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
[columnScale,columnNorm,construction]=buildRocketParentImageScale(packed,factor);
scaleFile=fullfile(folder,'parent_image_scale_evidence.mat');save(scaleFile,'columnScale','columnNorm','construction','packedFile','factorFile','sourceBinding','sourceIdentity','-v7');
measured=whos('packed','factor','columnScale','columnNorm');cacheBytes=sum([measured.bytes]);assert(cacheBytes<=1024*1024^2,'ParentLSQR:Budget','完整父矩阵/局部因子/尺度共同缓存超1024MiB。');
before=applyRocketPackedGraph(packed,initial);assert(abs(norm(rhs-before,'fro')/norm(rhs,'fro')-baselineRelativeCResidual)<1e-12);clear before;
[value,detail,coordinates]=solveRocketParentLeastSquares(packed,factor,rhs,initial,struct('maximumIterations',128,'relativeTolerance',1e-3),columnScale);
image=applyRocketPackedGraph(packed,value);remaining=rhs-image;assert(abs(norm(remaining,'fro')/norm(rhs,'fro')-detail.relativeCResidual)<1e-12);
save(fullfile(folder,'parent_least_squares_real_evidence.mat'),'rhs','initial','value','image','remaining','coordinates','detail','construction','scaleFile', ...
    'baseFile','factorFile','packedFile','sourceBinding','sourceIdentity','measured','cacheBytes','baselineRelativeCResidual','-v7');
report=struct('status','tested-direct-complete-parent-C-image-scaled-LSQR-with-actual-original-norm','sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity, ...
    'construction',construction,'innerReport',detail,'actualInputCacheArrayBytes',cacheBytes,'scaleFileSHA256',rocketStreamingFileSHA256(scaleFile), ...
    'separateBaselineCReplayApplications',1,'separateFinalCDiagnosticApplications',1,'oldSchurUsed',false,'fixedCRepresentativesUsed',false, ...
    'newRootMatrixAssemblies',0,'newLocalInverseConstructions',0,'newSelectedJApplications',0,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'newProductionSolvedAreaM2',0,'physicalSolution',false,'processPeakMemoryClaimed',false,'engineProcessID',feature('getpid'));
raw=jsonencode(report);assert(numel(raw)<128*1024,'ParentLSQR:Summary','完整父矩阵LSQR摘要超128KiB。');
fid=fopen(fullfile(folder,'parent_least_squares_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',raw);fprintf('%s\n',raw);
end
