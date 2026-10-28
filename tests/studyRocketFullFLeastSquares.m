function report=studyRocketFullFLeastSquares(folder)
%STUDYROCKETFULLFLEASTSQUARES 从同Q原最佳F端点做一次32步LSQR，并单列固定C父诊断。
% 模块一：复用原完整F/认证局部逆/右端，不重新选分组或更改旧近似S。
% 模块二：保存实际原F/父C总方向，不将子系统改善或程序退出0当成物理合格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);
baseFile=fullfile(root,'outputs','v83_validation','fullFreal01','full_F_real_evidence.mat');
assert(strcmp(rocketStreamingFileSHA256(baseFile),'1bf2bfef8fd8f5641f0e2441ab358991b553cd3a5ab77c3a96d15a3bf4b576d5'));
s=load(baseFile,'assemblyFile','graphFile','packedFile','sourceBinding','rhsF','valueF','Fimage','value','originalRHS');
assemblyFile=s.assemblyFile;graphFile=s.graphFile;packedFile=s.packedFile;sourceBinding=s.sourceBinding;
rhsF=s.rhsF;initial=s.valueF;baseImage=s.Fimage;baseValue=s.value;originalRHS=s.originalRHS;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceBinding.sourceFile),sourceBinding.sourceFileSHA256));
assert(strcmp(rocketStreamingFileSHA256(graphFile),'b8ac74781499c0054ea49af77f3bad41f6ae6600dea9045264bba0230c5e4a04'));
s=load(graphFile,'graph');graph=s.graph;clear s;s=load(assemblyFile,'reduction');factor=s.reduction.Ffactor;clear s;
sourceIdentity=struct;for key={'applyRocketPackedGraphTranspose','solveRocketSpatialFBlocksTranspose','solveRocketFullFLeastSquares','applyRocketPackedGraph','solveRocketSpatialFBlocks'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
check=applyRocketPackedGraph(graph,initial);assert(norm(check-baseImage,'fro')/max(1,norm(baseImage,'fro'))<1e-12);clear check;
[valueF,detail]=solveRocketFullFLeastSquares(graph,factor,rhsF,initial,struct('maximumIterations',32,'relativeTolerance',1e-3));
Fimage=applyRocketPackedGraph(graph,valueF);assert(abs(norm(rhsF-Fimage,'fro')/norm(rhsF,'fro')-detail.relativeFullFResidual)<1e-12);
save(fullfile(folder,'full_F_least_squares_inner_evidence.mat'),'valueF','Fimage','rhsF','initial','detail','baseFile','graphFile','assemblyFile','sourceBinding','sourceIdentity','-v7');
clear factor;s=load(packedFile,'packed');packed=s.packed;clear s;value=baseValue;isF=false(size(value,1),1);isF(graph.originalPhysicalCVs)=true;
value(graph.originalPhysicalCVs,:)=valueF;image=applyRocketPackedGraph(packed,value);remaining=originalRHS-image;
diagnostic=struct('relativeParentOriginalCResidual',norm(remaining,'fro')/norm(originalRHS,'fro'), ...
    'FRowsRelativeToOriginalRHS',norm(remaining(isF,:),'fro')/norm(originalRHS,'fro'),'CRowsRelativeToOriginalRHS',norm(remaining(~isF,:),'fro')/norm(originalRHS,'fro'), ...
    'FRowsSubmatrixConsistency',norm(remaining(graph.originalPhysicalCVs,:)-(rhsF-Fimage),'fro')/max(1,norm(originalRHS,'fro')));
assert(diagnostic.FRowsSubmatrixConsistency<1e-8);
save(fullfile(folder,'full_F_least_squares_real_evidence.mat'),'valueF','Fimage','rhsF','initial','detail','value','image','remaining','originalRHS','diagnostic', ...
    'baseFile','graphFile','assemblyFile','packedFile','sourceBinding','sourceIdentity','-v7');
report=struct('status','tested-real-bounded-complete-F-LSQR-with-fixed-C-parent-diagnostic','sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity, ...
    'fullFReport',detail,'parentDiagnostic',diagnostic,'separateBaselineReplayFApplications',1,'separateFinalDiagnosticFApplications',1, ...
    'originalCDiagnosticApplications',1,'newLocalInverseConstructions',0,'newApproximateSchurAssemblies',0,'newSelectedJApplications',0, ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'newProductionSolvedAreaM2',0,'physicalSolution',false,'engineProcessID',feature('getpid'));
fid=fopen(fullfile(folder,'full_F_least_squares_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
