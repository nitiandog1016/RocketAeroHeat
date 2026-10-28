function report=studyRocketRegularizedParentDirection(folder)
%STUDYROCKETREGULARIZEDPARENTDIRECTION 在同一Q20/C上限制局部更新，不重建旧筛查或改变接受历史。
% 模块一：复用已验当前Q20矩阵/因子，从同一新鲜方向作总方向正则化试算。
% 模块二：先核对完整解析模型；通过0.1后才另作原NS探针，不自动提交流场。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'));mkdir(folder);base=fullfile(root,'outputs','v83_validation','currentparent01');
p=load(fullfile(base,'current_packed_evidence.mat'),'packed','sourceBinding','sourceIdentity');f=load(fullfile(base,'current_factor_evidence.mat'),'factor');
d=load(fullfile(base,'current_parent_candidate_evidence.mat'),'rhs','value','cacheFile');c=load(fullfile(base,'current_image_scale_evidence.mat'),'columnScale');
assert(isequaln(p.sourceIdentity,rocketSelectedDirectionSourceIdentity)&& ...
    strcmp(rocketStreamingFileSHA256(p.sourceBinding.sourceFile),p.sourceBinding.sourceFileSHA256));
s=load(d.cacheFile,'Q','scale','gas');metric=buildRocketPrimitiveUpdateMetric(s.Q,s.scale,s.gas);clear s;
[value,detail,coordinates]=solveRocketRegularizedParentDirection(p.packed,f.factor,d.rhs,d.value,c.columnScale,metric, ...
    struct('maximumIterations',128,'regularizationStrength',.01));
sourceBinding=p.sourceBinding;sourceIdentity=p.sourceIdentity;rhs=d.rhs;initial=d.value;cacheFile=d.cacheFile;
save(fullfile(folder,'regularized_parent_candidate_evidence.mat'),'value','detail','coordinates','metric','rhs','initial', ...
    'sourceBinding','sourceIdentity','base','cacheFile','-v7');clear p f c d;
s=load(cacheFile,'cache');[fullModel,fullDetail]=applyCutCellSelectedDirection(s.cache,value);clear s;
fullResidual=norm(rhs+fullModel,'fro')/norm(rhs,'fro');
report=struct('status','bounded-current-Q20-total-primitive-update-regularized-parent-direction', ...
    'sourceBinding',sourceBinding,'sourceIdentity',sourceIdentity,'linearReport',detail, ...
    'completeSelectedModelRelativeDefect',fullResidual,'modelPointOneMet',fullResidual<=.1, ...
    'newSelectedModelApplications',1,'newActualNSEvaluations',0,'newFlowAcceptances',0, ...
    'physicalStateChanged',false,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'regularized_parent_model_evidence.mat'),'value','fullModel','fullDetail','report','-v7');
fid=fopen(fullfile(folder,'regularized_parent_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
