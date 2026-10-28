function report=studyRocketCoarsePostsmoothing(folder)
%STUDYROCKETCOARSEPOSTSMOOTHING 复用已核对的90维粗矩阵，检查固定后平滑的细层模型影响。
% 模块一：不重建RAP、不调用原NS；保持真实原Q和同一完整空间模型。
% 模块二：逐次记录粗/细残差，选项只能为整次试验固定次数，不按每个右端自选。
% 模块三：正反馈或粗残差变化只作预条件研究，不能视为流场稳态。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'CoarsePost:Overwrite','不得覆盖后平滑研究。');
sourceFile=fullfile(root,'outputs','v83_validation','coarserocket01','volume_coarse_real_evidence.mat');
s=load(sourceFile,'coarse','rhs','update','afterCorrection');coarse=s.coarse;rhs=s.rhs;initial=s.update;priorRemaining=s.afterCorrection;clear s;
d=load(fullfile(root,'outputs','v83_validation','krylov08','wall_velocity_krylov_evidence.mat'),'result','gas','boundary');
Q=d.result.Q;scale=d.result.residualScale;gas=d.gas;boundary=d.boundary;info=struct;
for key={'primitive','gradient','viscousFaceFlux'},info.(key{1})=d.result.info.(key{1});end;clear d;
assert(isequaln(Q,coarse.baseQ));cut=coarse.space.geometry;t=tic;cache=buildCutCellSelectedDirectionCache(Q,info,cut,gas,boundary,scale);cacheSeconds=toc(t);
graph=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-7,1);modelCalls=0;
remaining=rhs+applyCutCellSelectedDirection(cache,initial);modelCalls=modelCalls+1;
assert(norm(remaining-priorRemaining,'fro')/max(1,norm(priorRemaining,'fro'))<1e-11,'CoarsePost:Identity','没有复现原粗校正后的完整模型。');
updates=cell(5,1);updates{1}=initial;history=zeros(5,3);times=zeros(4,1);history(1,:)=measure(remaining);update=initial;
for sweep=1:4
    t=tic;correction=solveRocketShellGraphPreconditioner(graph,remaining);update=update+correction;
    remaining=rhs+applyCutCellSelectedDirection(cache,update);modelCalls=modelCalls+1;times(sweep)=toc(t);
    updates{sweep+1}=update;history(sweep+1,:)=measure(remaining);
end
% 最后一次固定次数的高层入口须与逐步重放一致。
replayed=solveRocketTwoLevelPreconditioner(graph,coarse,cache,rhs,4);modelCalls=modelCalls+5;
replayError=norm(replayed-update,'fro')/max(1,norm(update,'fro'));assert(replayError<1e-11,'CoarsePost:Replay','固定次数入口不符合重放。');
report=struct('status','measured-fixed-postsmoothing-reusing-real-Galerkin-coarse-matrix','coarseUnknowns',coarse.space.coarseUnknowns, ...
    'historyColumns',{{'fineModelRelativeResidual','restrictedModelRelativeResidual','maximumAbsoluteModelResidual'}}, ...
    'fixedPostSweeps',(0:4)','history',history,'perPostSeconds',times,'directionCacheBuildSeconds',cacheSeconds, ...
    'fixedFourPostReplayRelativeDifference',replayError,'newCoarseConstructionApplications',0,'modelApplications',modelCalls, ...
    'actualOriginalNSEvaluations',0,'newAcceptedIterations',0,'physicalMeshChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
mkdir(folder);save(fullfile(folder,'coarse_post_evidence.mat'),'report','sourceFile','rhs','initial','updates','history','replayed','-v7');
fid=fopen(fullfile(folder,'coarse_post_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    function row=measure(value)
        restricted=coarse.space.restriction*value;
        row=[norm(value,'fro')/norm(rhs,'fro'),norm(restricted,'fro')/max(1,norm(coarse.space.restriction*rhs,'fro')),max(abs(value),[],'all')];
    end
end
