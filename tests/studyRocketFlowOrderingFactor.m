function report=studyRocketFlowOrderingFactor(folder)
%STUDYROCKETFLOWORDERINGFACTOR 冻结真实C按给定流向重排，先只检验分解和三条右端。
% 模块一：完整C值与物理CV保持，几何来源与主链端点均有明确字节绑定。
% 模块二：先量化固定三角解的放大/剩余，不盲目重复32方向内部迭代。
% 模块三：所有结果是代数诊断，不调用新的完整J/NS或推进物理状态。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FlowFactor:Overwrite','不得覆盖真实流向排序研究。');mkdir(folder);
packedFile=fullfile(root,'outputs','v83_validation','blockilureal01','steady_graph_packed_evidence.mat');s=load(packedFile,'packed','sourceFile','sourceBinding');
packed=s.packed;sourceFile=s.sourceFile;sourceBinding=s.sourceBinding;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceFile),sourceBinding.sourceFileSHA256));
geometryFile=fullfile(root,'outputs','v83_validation','flowordergeometry01','flow_ordering_geometry.mat');g=load(geometryFile,'report','centroid','unitFlowDirection');
assert(g.report.cellCount==packed.cellCount&&g.report.sourceFaces==packed.sourceFaces&&g.report.layers==packed.layers&& ...
    strcmp(g.report.sourceSHA256,rocketStreamingFileSHA256(g.report.sourceFile)));
anchorFile=fullfile(root,'outputs','v83_validation','innergraphreal01','inner_graph_real_evidence.mat');s=load(anchorFile,'sourceFile','remaining');
assert(strcmp(sourceFile,s.sourceFile));remaining=s.remaining;clear s;
stage='ordering';factorAvailable=false;rejected=false;errorIdentifier='';cases={};orderingReport=struct;
try
    ordered=reorderRocketBlockGraphByFlow(packed,g.centroid,g.unitFlowDirection);clear g;
    orderingReport=ordered.flowReorderingReport;pattern=struct('rowPointer',ordered.rowPointer,'column',ordered.column, ...
        'diagonalPosition',ordered.diagonalPosition,'order',ordered.order,'inverseOrder',ordered.inverseOrder,'flowReorderingReport',orderingReport);
    save(fullfile(folder,'flow_ordered_pattern_evidence.mat'),'pattern','packedFile','sourceFile','sourceBinding','geometryFile','-v7');
    probe=reshape(sin((1:numel(remaining))'*.37),size(remaining));reference=applyRocketPackedGraph(packed,probe);actual=applyRocketPackedGraph(ordered,probe);
    orderingActionError=norm(actual-reference,'fro')/max(1,norm(reference,'fro'));assert(orderingActionError<2e-10,'FlowFactor:C','流向重排改变了完整物理C作用。');clear packed probe actual reference;
    stage='factorization';ordered=factorRocketBlockILU0(ordered);factor=ordered;clear ordered;factorAvailable=true;
    factorIdentity=struct;for key={'reorderRocketBlockGraphByFlow','factorRocketBlockILU0','solveRocketBlockILU0'},factorIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
    factorFile=fullfile(folder,'flow_ordered_factor_evidence.mat');
    save(factorFile,'factor','factorIdentity','packedFile','anchorFile','geometryFile','sourceFile','sourceBinding','-v7');
    s=load(packedFile,'packed');packed=s.packed;clear s;stage='three-triangular-action-probes';
    summaries=cell(3,1);cases=cell(3,1);names={'sameQRemaining','sineProbe','cosineProbe'};
    for j=1:3
        if j==1,rhs=remaining;elseif j==2,rhs=reshape(sin((1:numel(remaining))'*.37),size(remaining));else,rhs=reshape(cos((1:numel(remaining))'*.23),size(remaining));end
        timer=tic;[value,detail]=solveRocketBlockILU0(factor,rhs);seconds=toc(timer);image=applyRocketPackedGraph(packed,value);
        summary=struct('name',names{j},'relativeCResidual',norm(rhs-image,'fro')/norm(rhs,'fro'), ...
            'directionNormToRhsSeconds',norm(value,'fro')/norm(rhs,'fro'),'triangularSeconds',seconds,'originalNSAcceptanceClaimed',false);
        summaries{j}=summary;cases{j}=struct('rhs',rhs,'value',value,'Cvalue',image,'summary',summary,'detail',detail);
    end
    stage='bounded-action-probes-complete';
catch ex
    if ~(startsWith(ex.identifier,'BlockFlowOrder:')||startsWith(ex.identifier,'BlockILU:')),rethrow(ex);end
    rejected=true;errorIdentifier=ex.identifier;
end
report=struct('status','tested-real-prescribed-flow-ordering-factor-and-fixed-triangular-actions', ...
    'stage',stage,'factorAvailable',factorAvailable,'factorRejected',rejected,'errorIdentifier',errorIdentifier, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'newCoarseConstructions',0,'newFlowAcceptances',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
if isfield(orderingReport,'method'),report.orderingReport=rmfield(orderingReport,'newToOldBlockPosition');end
if factorAvailable,report.factorizationReport=factor.factorizationReport;report.orderingActionRelativeDifference=orderingActionError;end
if ~rejected
    report.summaries=summaries;
    save(fullfile(folder,'flow_ordering_actions_evidence.mat'),'report','cases','packedFile','factorFile','geometryFile','sourceFile','sourceBinding','anchorFile','-v7');
end
fid=fopen(fullfile(folder,'flow_ordering_factor_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
