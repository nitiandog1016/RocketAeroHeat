function report=studyRocketFullFImageScale(folder)
%STUDYROCKETFULLFIMAGESCALE 构造真实完整像列归一，对同一F端点作显式有界LSQR并核对父C。
% 本研究使用新固定右侧归一，最大128步；不与原32步比较等工作量加速，不改旧S/NS。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
baseFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_real_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(baseFile),'1bf2bfef8fd8f5641f0e2441ab358991b553cd3a5ab77c3a96d15a3bf4b576d5'));
s=load(baseFile,'assemblyFile','graphFile','packedFile','sourceBinding','rhsF','valueF','Fimage','value','originalRHS');
assemblyFile=s.assemblyFile;graphFile=s.graphFile;packedFile=s.packedFile;sourceBinding=s.sourceBinding;rhsF=s.rhsF;initial=s.valueF;baseImage=s.Fimage;baseValue=s.value;originalRHS=s.originalRHS;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
assert(strcmp(rocketStreamingFileSHA256(graphFile),'b8ac74781499c0054ea49af77f3bad41f6ae6600dea9045264bba0230c5e4a04'));
s=load(graphFile,'graph');graph=s.graph;clear s;s=load(assemblyFile,'reduction');factor=s.reduction.Ffactor;clear s;
sourceIdentity=struct;for key={'buildRocketFullFImageScale','solveRocketFullFLeastSquares','applyRocketPackedGraph','applyRocketPackedGraphTranspose','solveRocketSpatialFBlocks','solveRocketSpatialFBlocksTranspose'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
[columnScale,columnNorm,construction]=buildRocketFullFImageScale(graph,factor);
scaleFile=fullfile(folder,'full_F_image_scale_evidence.mat');save(scaleFile,'columnScale','columnNorm','construction','graphFile','assemblyFile','sourceBinding','sourceIdentity','-v7');
check=applyRocketPackedGraph(graph,initial);assert(norm(check-baseImage,'fro')/max(1,norm(baseImage,'fro'))<1e-12);clear check;
[valueF,detail]=solveRocketFullFLeastSquares(graph,factor,rhsF,initial,struct('maximumIterations',128,'relativeTolerance',1e-3),columnScale);
Fimage=applyRocketPackedGraph(graph,valueF);assert(abs(norm(rhsF-Fimage,'fro')/norm(rhsF,'fro')-detail.relativeFullFResidual)<1e-12);
save(fullfile(folder,'full_F_image_scaled_inner_evidence.mat'),'valueF','Fimage','rhsF','initial','detail','scaleFile','baseFile','graphFile','assemblyFile','sourceBinding','sourceIdentity','-v7');
clear factor;s=load(packedFile,'packed');packed=s.packed;clear s;value=baseValue;isF=false(size(value,1),1);isF(graph.originalPhysicalCVs)=true;value(graph.originalPhysicalCVs,:)=valueF;
image=applyRocketPackedGraph(packed,value);remaining=originalRHS-image;
diagnostic=struct('relativeParentOriginalCResidual',norm(remaining,'fro')/norm(originalRHS,'fro'), ...
    'FRowsRelativeToOriginalRHS',norm(remaining(isF,:),'fro')/norm(originalRHS,'fro'),'CRowsRelativeToOriginalRHS',norm(remaining(~isF,:),'fro')/norm(originalRHS,'fro'), ...
    'FRowsSubmatrixConsistency',norm(remaining(graph.originalPhysicalCVs,:)-(rhsF-Fimage),'fro')/max(1,norm(originalRHS,'fro')));
assert(diagnostic.FRowsSubmatrixConsistency<1e-8);
save(fullfile(folder,'full_F_image_scaled_real_evidence.mat'),'valueF','Fimage','rhsF','initial','detail','construction','value','image','remaining','originalRHS','diagnostic', ...
    'scaleFile','baseFile','graphFile','assemblyFile','packedFile','sourceBinding','sourceIdentity','-v7');
report=struct('status','tested-real-complete-image-column-scaled-F-LSQR-with-separate-fixed-C-parent-diagnostic','sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity, ...
    'construction',construction,'fullFReport',detail,'parentDiagnostic',diagnostic,'scaleFileSHA256',rocketStreamingFileSHA256(scaleFile), ...
    'separateBaselineReplayFApplications',1,'separateFinalDiagnosticFApplications',1,'originalCDiagnosticApplications',1, ...
    'newLocalInverseConstructions',0,'newApproximateSchurAssemblies',0,'newSelectedJApplications',0,'newActualNSEvaluations',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
fid=fopen(fullfile(folder,'full_F_image_scaled_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
