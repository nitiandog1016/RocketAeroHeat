function report=studyRocketFieldSchurActions(folder)
%STUDYROCKETFIELDSCHURACTIONS 在真实冻结C上先验证五量Schur构造和三条固定右端。
% 模块一：复用原symrcm完整C和同Q实际剩余，原NS/粗像不重新构造。
% 模块二：局部逆、因子及预算失败明确保存身份；不盲目开始内部迭代。
% 模块三：所有有限结果仍是模型诊断，不更新主链压力/热流或已解面积。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'FieldSchurStudy:Overwrite','不得覆盖真实Schur研究。');mkdir(folder);
packedFile=fullfile(root,'outputs','v83_validation','blockilureal01','steady_graph_packed_evidence.mat');s=load(packedFile,'packed','sourceFile','sourceBinding');packed=s.packed;sourceFile=s.sourceFile;sourceBinding=s.sourceBinding;clear s;
assert(strcmp(rocketStreamingFileSHA256(sourceFile),sourceBinding.sourceFileSHA256));
anchorFile=fullfile(root,'outputs','v83_validation','innergraphreal01','inner_graph_real_evidence.mat');s=load(anchorFile,'sourceFile','remaining');assert(strcmp(sourceFile,s.sourceFile));remaining=s.remaining;clear s;
stage='Schur-construction-and-factorization';rejected=false;errorIdentifier='';message='';preAvailable=false;summaries={};
try
    pre=buildRocketFieldSchurPreconditioner(packed);preAvailable=true;preFile=fullfile(folder,'field_schur_preconditioner_evidence.mat');sourceIdentity=struct;
    for key={'buildRocketFieldSchurPreconditioner','factorRocketFieldILU0','solveRocketFieldILU0','solveRocketFieldSchurPreconditioner'},sourceIdentity.(key{1})=rocketStreamingFileSHA256(which(key{1}));end
    save(preFile,'pre','sourceIdentity','packedFile','anchorFile','sourceFile','sourceBinding','-v7');
    stage='three-fixed-Schur-action-probes';cases=cell(3,1);summaries=cell(3,1);names={'sameQRemaining','sineProbe','cosineProbe'};
    for j=1:3
        if j==1,rhs=remaining;elseif j==2,rhs=reshape(sin((1:numel(remaining))'*.37),size(remaining));else,rhs=reshape(cos((1:numel(remaining))'*.23),size(remaining));end
        timer=tic;[value,detail]=solveRocketFieldSchurPreconditioner(pre,rhs);seconds=toc(timer);image=applyRocketPackedGraph(packed,value);
        summary=struct('name',names{j},'relativeCResidual',norm(rhs-image,'fro')/norm(rhs,'fro'), ...
            'directionNormToRhsSeconds',norm(value,'fro')/norm(rhs,'fro'),'fixedSchurSeconds',seconds,'originalNSAcceptanceClaimed',false);
        summaries{j}=summary;cases{j}=struct('rhs',rhs,'value',value,'Cvalue',image,'summary',summary,'detail',detail);
    end
    stage='bounded-Schur-action-probes-complete';
    save(fullfile(folder,'field_schur_actions_evidence.mat'),'cases','preFile','packedFile','sourceFile','sourceBinding','anchorFile','-v7');
catch ex
    if ~(startsWith(ex.identifier,'FieldSchur:')||startsWith(ex.identifier,'FieldILU:')),rethrow(ex);end
    rejected=true;errorIdentifier=ex.identifier;message=ex.message;
end
report=struct('status','tested-real-five-conservation-field-Schur-construction-and-fixed-action-study', ...
    'stage',stage,'preconditionerAvailable',preAvailable,'candidateRejected',rejected,'errorIdentifier',errorIdentifier,'message',message, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'newCoarseConstructions',0,'newFlowAcceptances',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
if preAvailable,report.constructionReport=pre.constructionReport;report.momentumFactorReport=pre.momentumFactor.factorizationReport;report.schurFactorReport=pre.schurFactor.factorizationReport;end
if ~rejected,report.summaries=summaries;end
fid=fopen(fullfile(folder,'field_schur_real_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
