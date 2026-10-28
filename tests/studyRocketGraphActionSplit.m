function report=studyRocketGraphActionSplit(folder,sourceCase,parentCase)
%STUDYROCKETGRAPHACTIONSPLIT 用真实保存方向量化当前完整J与Jgraph的耦合遗漏。
% 模块一：保留同一Q及全部原面/CV；读取已存图矩阵，不重算792列粗像。
% 模块二：分别检查总Newton方向、最新增强增量和归一化实际模型剩余。
% 模块三：记录Euler/特征边界、应力、导热和两项功的误差与相消，不宣称完整谱结论。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'),fullfile(root,'tests'));setupCutCellFlow;
assert(~exist(folder,'dir'),'GraphSplitStudy:Overwrite','不得覆盖分项研究。');mkdir(folder);
assert(~isempty(regexp(sourceCase,'^krylov[0-9]+$','once'))&&~isempty(regexp(parentCase,'^krylov[0-9]+$','once')));
sourceFile=fullfile(root,'outputs','v83_validation',sourceCase,'wall_velocity_krylov_evidence.mat');
parentFile=fullfile(root,'outputs','v83_validation',parentCase,'wall_velocity_krylov_evidence.mat');
s=load(sourceFile,'result','gas','boundary','report');r=s.result;gas=s.gas;boundary=s.boundary;
p=load(parentFile,'result');assert(isequaln(r.Q,r.initialQ)&&isequaln(p.result.Q,r.Q)&&r.acceptedNonlinearIterations==0&&p.result.acceptedNonlinearIterations==0, ...
    'GraphSplitStudy:Q','只允许同Q且没有接受新状态的两个端点。');
Q=r.Q;scale=r.residualScale;rhs=r.initialRate./scale;operator=r.firstOperator;sourceActualNS=s.report.totalNewEquationResidualEvaluations;
savedOriginalNSAction=r.linearLog{end}.steadyDirectionEvidence.value;
directions={r.linearLog{end}.steadyDirectionEvidence.directionU, ...
    r.linearLog{end}.steadyDirectionEvidence.directionU-p.result.linearLog{end}.steadyDirectionEvidence.directionU, ...
    rhs+savedOriginalNSAction};directionNames={'totalNewton','lastRecycledIncrement','normalizedRemaining'};
directions{3}=directions{3}/max(abs(directions{3}),[],'all');
assert(isequaln(p.result.residualScale,scale)&&isequaln(p.result.initialRate,r.initialRate));
baseInfo=struct;for key={'primitive','gradient','viscousFaceFlux'},baseInfo.(key{1})=r.initialInfo.(key{1});end;clear p r s;
d=load(fullfile(root,'outputs','v80_validation','domain02','shell_domain_evidence.mat'),'domain');cut=d.domain;clear d;
cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);nf=operator.sourceFaces;nl=operator.layers;
save(fullfile(folder,'graph_split_context.mat'),'sourceFile','parentFile','Q','cut','gas','boundary','scale','baseInfo','rhs','savedOriginalNSAction','sourceActualNS','-v7');
summaries=cell(numel(directions),1);
for j=1:numel(directions)
    timer=tic;directionU=directions{j};[full,graph,detail]=splitRocketGraphDirection(cache,operator,directionU);
    fullAction=detail.fullTotal;graphAction=detail.graphTotal;error=graphAction-fullAction;names=fieldnames(full);
    errorNorm2=sum(error.^2,'all');assert(errorNorm2>0);componentGram=zeros(numel(names));componentRows=zeros(numel(names),5);
    for a=1:numel(names)
        ea=graph.(names{a})-full.(names{a});
        componentRows(a,:)=[norm(full.(names{a}),'fro'),norm(graph.(names{a}),'fro'),norm(ea,'fro'),sum(ea.*error,'all')/errorNorm2, ...
            norm(error-ea,'fro')/norm(fullAction,'fro')];
        for b=1:numel(names),eb=graph.(names{b})-full.(names{b});componentGram(a,b)=sum(ea.*eb,'all')/errorNorm2;end
    end
    layerError=zeros(nl,1);for layer=1:nl,ids=(layer-1)*nf+(1:nf);layerError(layer)=sum(error(ids,:).^2,'all')/errorNorm2;end
    [~,order]=sort(sum(error.^2,2),'descend');top=order(1:20);
    summary=struct('name',directionNames{j},'graphActionRelativeDifference',norm(error,'fro')/norm(fullAction,'fro'), ...
        'componentNames',{names},'componentColumns',{{'fullJNorm','graphJNorm','differenceNorm','signedContributionToTotalDifferenceSquaredNorm','remainingDifferenceRelativeToFullAfterExactReplacement'}}, ...
        'componentRows',componentRows,'componentDifferenceGram',componentGram,'equationDifferenceRelativeToFull',vecnorm(error,2,1)./vecnorm(fullAction,2,1), ...
        'layerDifferenceSquaredNormFraction',layerError,'largestDifferenceCellIDs',top,'largestDifferenceOriginalFaceIDs',mod(top-1,nf)+1, ...
        'largestDifferenceLayerIDs',floor((top-1)/nf)+1,'graphAssemblyRelativeDifference',detail.graphAssemblyRelativeDifference, ...
        'fullModelRelativeDefect',norm(rhs+fullAction,'fro')/norm(rhs,'fro'),'graphModelRelativeDefect',norm(rhs+graphAction,'fro')/norm(rhs,'fro'), ...
        'elapsedSeconds',toc(timer),'singleDirectionOnly',true,'physicalSolution',false,'newProductionSolvedAreaM2',0);
    if j==1
        summary.savedOriginalNSDirectionRelativeDifference=norm(fullAction-savedOriginalNSAction,'fro')/norm(savedOriginalNSAction,'fro');
        assert(summary.savedOriginalNSDirectionRelativeDifference<2e-6,'GraphSplitStudy:NS','当前完整作用没有复现原NS方向。');
    end
    summaries{j}=summary;fprintf('GRAPH_SPLIT %s relative=%.9g assembly=%.4g\n',summary.name,summary.graphActionRelativeDifference,summary.graphAssemblyRelativeDifference);
    save(fullfile(folder,['graph_split_' directionNames{j} '.mat']),'directionU','full','graph','detail','summary','-v7');
    clear full graph detail fullAction graphAction error;
end
report=struct('status','measured-current-real-J-versus-graph-components-on-three-same-Q-directions', ...
    'sourceCase',sourceCase,'parentCase',parentCase,'sourceFaces',nf,'cellCount',size(Q,1),'directionCount',numel(directions), ...
    'completeSelectedActionApplications',numel(directions),'newGraphBuilds',0,'newCoarseColumnApplications',0, ...
    'newActualNSEvaluations',0,'newFlowAcceptances',0,'mainActualNSUnchanged',sourceActualNS,'physicalSolution',false, ...
    'newProductionSolvedAreaM2',0,'componentInterpretation','Signed contributions include cancellation; exact replacement is only an action diagnostic, not a solved correction or conditioning result.', ...
    'summaries',{summaries},'engineProcessID',feature('getpid'));
fid=fopen(fullfile(folder,'graph_split_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(rmfield(report,'summaries')));
end
