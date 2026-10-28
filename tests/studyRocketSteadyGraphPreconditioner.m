function report=studyRocketSteadyGraphPreconditioner(folder,sourceCase)
%STUDYROCKETSTEADYGRAPHPRECONDITIONER 在当前完整模型剩余上有界评估整体去正则C。
% 模块一：同一Q、原网格及粗像；只改变预条件矩阵，时间与物理方程不变。
% 模块二：分解失败保留明确拒绝记录；成功才比较一次扫描和内部C求解。
% 模块三：模型预测不替代原NS总方向和候选验收，不更新压力/热流。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SteadyGraphStudy:Overwrite','不得覆盖真实稳态图研究。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation',sourceCase,'wall_velocity_krylov_evidence.mat');s=load(sourceFile,'result','gas','boundary');r=s.result;
Q=r.Q;scale=r.residualScale;rhs=r.initialRate./scale;gas=s.gas;boundary=s.boundary;original=r.firstOperator;
assert(r.acceptedNonlinearIterations==0&&isequaln(Q,r.initialQ));direction=r.linearLog{end}.steadyDirectionEvidence.directionU;
coarseFile=r.options.coarseOperatorFile;optionsSource=r.options;baseInfo=struct;
for key={'primitive','gradient','viscousFaceFlux'},baseInfo.(key{1})=r.initialInfo.(key{1});end;clear r s;
try,operator=buildRocketSteadyGraphPreconditioner(original);
catch ex
    report=struct('status','rejected-real-unshifted-line-factor-before-model-or-flow-update','errorIdentifier',ex.identifier, ...
        'newSelectedJApplications',0,'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
    save(fullfile(folder,'steady_graph_factor_rejection.mat'),'report','sourceFile','-v7');writeReport;return;
end
clear original;
t=load(coarseFile,'coarse');space=t.coarse.space;primitiveSpace=t.coarse.mixedLift.primitiveSpace;primitiveSpace.geometry=space.geometry;cut=space.geometry;clear t;
cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);clear baseInfo;
initialModel=applyCutCellSelectedDirection(cache,direction);remaining=rhs+initialModel;
[coarse,reuseReport]=loadRocketLocalCoarseOperator(coarseFile,cache,space,primitiveSpace,optionsSource.coarseMixedColumns,optionsSource.coarseLocalColumns,256*1024^2);clear space primitiveSpace optionsSource;
options=struct('maximumIterations',32,'relativeTolerance',1e-3,'maximumBasisBytes',1024*1024^2,'saveBasisEvidence',false);
[single,singleReport]=solveRocketWlsNormalPreconditioner(operator,remaining);[inner,innerReport]=solveRocketInnerGraphKrylov(operator,remaining,options);
values={single,inner};names={'steadySingleScan','steadyInner32'};evidence=cell(2,1);summaries=cell(2,1);
for j=1:2
    value=values{j};Cvalue=applyRocketWlsNormalOperator(operator,value);Avalue=-applyCutCellSelectedDirection(cache,value);
    coarseRight=remaining-Avalue;weights=coarse.upper\(coarse.upper'\(coarse.image'*coarseRight(:)));
    [correction,coarseReport]=solveRocketLocalCoarseSystem(coarse,coarseRight);total=value+correction;totalAction=-applyCutCellSelectedDirection(cache,total);
    predicted=Avalue+reshape(coarse.image*weights,size(Q));assert(norm(predicted-totalAction,'fro')/max(1,norm(totalAction,'fro'))<1e-9);
    alpha=dot(remaining(:),totalAction(:))/dot(totalAction(:),totalAction(:));best=remaining-alpha*totalAction;
    summary=struct('name',names{j},'relativeCResidual',norm(remaining-Cvalue,'fro')/norm(remaining,'fro'), ...
        'fullADirectionRelativeResidual',norm(remaining-Avalue,'fro')/norm(remaining,'fro'), ...
        'afterCoarseRelativeResidual',norm(remaining-totalAction,'fro')/norm(remaining,'fro'), ...
        'bestScalarAfterCoarse',alpha,'predictedTotalDefectAfterBestScalar',norm(best,'fro')/norm(rhs,'fro'), ...
        'remainingSquaredNormRemovedByBestScalar',1-norm(best,'fro')^2/norm(remaining,'fro')^2, ...
        'steadyGraphToFullAActionRelativeDifference',norm(Cvalue-Avalue,'fro')/max(1,norm(Avalue,'fro')), ...
        'originalNSAcceptanceClaimed',false,'fixedOuterLinearMapClaimed',j==1);
    evidence{j}=struct('value',value,'Cvalue',Cvalue,'Avalue',Avalue,'correction',correction,'coarseWeights',weights, ...
        'total',total,'totalAction',totalAction,'coarseReport',coarseReport,'summary',summary);summaries{j}=summary;
    fprintf('STEADY_GRAPH %s C=%.8g Afull=%.8g afterCoarse=%.8g bestTotal=%.9g\n',names{j},summary.relativeCResidual,summary.fullADirectionRelativeResidual,summary.afterCoarseRelativeResidual,summary.predictedTotalDefectAfterBestScalar);
end
report=struct('status','tested-real-current-remaining-all-five-unshifted-steady-graph-model-candidate', ...
    'sourceCase',sourceCase,'cellCount',size(Q,1),'initialCurrentModelDefect',norm(remaining,'fro')/norm(rhs,'fro'), ...
    'factorReport',operator.steadyGraphPreconditionerReport,'summaries',{summaries},'newSelectedJApplications',5, ...
    'coarseReuseProbeModelApplications',1,'newCoarseConstructionApplications',0,'newActualNSEvaluations',0, ...
    'newFlowAcceptances',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_graph_real_evidence.mat'),'report','sourceFile','coarseFile','Q','scale','rhs','direction','initialModel','remaining', ...
    'options','singleReport','innerReport','reuseReport','evidence','-v7');writeReport;
    function writeReport
        fid=fopen(fullfile(folder,'steady_graph_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
    end
end
