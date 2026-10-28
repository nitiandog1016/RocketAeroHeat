function report=refineRocketSchurDiagonalDiagnostics(folder,sourceFile)
%REFINEROCKETSCHURDIAGONALDIAGNOSTICS 仅修复保存的局部奇异值诊断，不重演任何回代。
% 模块一：保留原失败MAT与全部Schur/动量中间量，只在新目录写修复诊断。
% 模块二：记录旧诊断和来源，原方程/原因子/拒绝状态及物理场不改变。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FieldRatioRefine:Overwrite','不得覆盖诊断修复证据。');mkdir(folder);s=load(sourceFile);sourceSHA256=rocketStreamingFileSHA256(sourceFile);
old=s.report.localDiagnostics;s.sigmaRatio=rocketTwoFieldSingularRatio(s.diagonal);[sorted,ids]=sort(s.sigmaRatio);worst=ids(1:min(12,numel(ids)));cv=s.algebraicOrder(worst);nf=s.sourceBinding.sourceFaces;
s.report.localDiagnostics=struct('scaledDiagonalSingularRatioRange',[min(s.sigmaRatio) max(s.sigmaRatio)],'ratiosBelow1e8',nnz(s.sigmaRatio<1e-8), ...
    'ratiosBelow1e12',nnz(s.sigmaRatio<1e-12),'worstAlgebraicRows',worst,'worstOriginalCV',cv,'worstOriginalFace',mod(cv-1,nf)+1, ...
    'worstLayers',floor((cv-1)/nf)+1,'worstRatios',sorted(1:numel(worst)),'globalNumericalRankProven',false);
s.refinementReport=struct('status','saved-diagonal-diagnostic-only-refinement-no-momentum-or-Schur-replay', ...
    'sourceFile',sourceFile,'sourceSHA256',sourceSHA256,'originalLocalDiagnostics',old, ...
    'helperSHA256',rocketStreamingFileSHA256(which('rocketTwoFieldSingularRatio')),'newMomentumSolves',0,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'consistent_schur_real_evidence.mat'),'-struct','s','-v7');report=s.report;
fid=fopen(fullfile(folder,'consistent_schur_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(s.refinementReport));
end
