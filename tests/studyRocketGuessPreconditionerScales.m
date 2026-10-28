function report=studyRocketGuessPreconditionerScales(folder,sourceFolder)
%STUDYROCKETGUESSPRECONDITIONERSCALES 用同一保存剩余比较预条件时间尺度，不更新Q或重新构造粗像。
% 模块一：实际方程始终是-J，tau只改变图扫描；原大tau放大证据不能直接外推新粗空间。
% 模块二：每候选一次预条件模型加一次当前完整作用，以最优标量投影测一方向收益。
% 模块三：图残差与完整空间投影分列，失败候选不补零；本研究没有新NS或接受步。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'GuessScale:Overwrite','不得覆盖预条件尺度探针。');mkdir(folder);
sourceFile=fullfile(sourceFolder,'wall_velocity_krylov_evidence.mat');s=load(sourceFile,'result','gas','boundary');r=s.result;gas=s.gas;boundary=s.boundary;
assert(r.acceptedNonlinearIterations==0&&isequaln(r.Q,r.initialQ));Q=r.Q;scale=r.residualScale;rhs=r.initialRate./scale;
remaining=rhs+r.linearLog{end}.steadyDirectionEvidence.value;coarseFile=r.options.coarseOperatorFile;
t=load(coarseFile,'coarse');space=t.coarse.space;primitiveSpace=t.coarse.mixedLift.primitiveSpace;primitiveSpace.geometry=space.geometry;cut=space.geometry;
info=struct;for key={'primitive','gradient','viscousFaceFlux'},info.(key{1})=r.initialInfo.(key{1});end;clear r s t;
cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);clear info;
columns=struct;options=load(sourceFile,'result');options=options.result.options;columns=options.coarseLocalColumns;
[coarse,reuseReport]=loadRocketLocalCoarseOperator(coarseFile,cache,space,primitiveSpace,options.coarseMixedColumns,columns,256*1024^2);clear options primitiveSpace space;
taus=[1e-8 3e-8 1e-7 2e-7 3e-7 1e-6];directions=cell(size(taus));products=directions;details=directions;errors=repmat({''},size(taus));history=nan(numel(taus),7);modelApplications=0;
for j=1:numel(taus)
    timer=tic;
    try
        graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,taus(j),1);
        [direction,detail]=solveRocketMinResidualCoarsePreconditioner(graph,coarse,cache,remaining);modelApplications=modelApplications+1;
        action=-applyCutCellSelectedDirection(cache,direction);modelApplications=modelApplications+1;
        coefficient=dot(remaining(:),action(:))/dot(action(:),action(:));prediction=remaining-coefficient*action;
        relative=norm(prediction,'fro')/norm(rhs,'fro');assert(isfinite(relative)&&relative<=norm(remaining,'fro')/norm(rhs,'fro')+1e-12);
        history(j,:)=[taus(j) detail.graphRelativeResidual coefficient relative 1-norm(prediction,'fro')^2/norm(remaining,'fro')^2 norm(direction,Inf) toc(timer)];
        directions{j}=direction;products{j}=action;details{j}=detail;
        fprintf('GUESS_SCALE tau=%.8g graph=%.8g oneDirectionDefect=%.9g gain=%.8g\n',taus(j),detail.graphRelativeResidual,relative,history(j,5));
    catch ex
        errors{j}=ex.identifier;fprintf('GUESS_SCALE_REJECT tau=%.8g error=%s\n',taus(j),ex.identifier);
    end
    clear graph;
end
[best,bestIndex]=min(history(:,4));assert(isfinite(best),'GuessScale:NoCandidate','没有有限候选。');
report=struct('status','tested-same-remaining-preconditioner-scale-one-direction-model-screen-only', ...
    'initialRemainingRelativeDefect',norm(remaining,'fro')/norm(rhs,'fro'),'bestPseudoDt',taus(bestIndex), ...
    'bestOneDirectionPredictedDefect',best,'modelApplications',modelApplications,'reuseVerificationModelApplications',1, ...
    'newCoarseConstructionApplications',0,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'noPhysicalEquationOrStateChange',true,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'guess_scale_evidence.mat'),'report','history','taus','directions','products','details','errors','remaining','rhs','sourceFile','coarseFile','reuseReport','-v7');
fid=fopen(fullfile(folder,'guess_scale_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
