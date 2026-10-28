function report=summarizeRocketBlockStudies(path)
%SUMMARIZEROCKETBLOCKSTUDIES 汇集同Q实际拒绝证据供科研工作台显示，不修改场量。
% 模块一：逐个读取原拒绝端点，实际C剩余、正交/关系及调用数来自保存向量。
% 模块二：绑定当前原NS文件字节，所有研究均为0新NS/0物理解更新。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(path,'file'),'BlockStudySummary:Overwrite','不得覆盖研究摘要。');val=fullfile(root,'outputs','v83_validation');
sourceFile=fullfile(val,'krylov17','wall_velocity_krylov_evidence.mat');sourceSHA=rocketStreamingFileSHA256(sourceFile);
folders={'blockilufailurediag01','scaledilureal01','blockfillreal01'};
names={'原模式块ILU0 / 两遍正交','原模式块ILU0 / 缩放复正交','一级填充块ILU / 缩放复正交'};
files={'arnoldi_failure_evidence.mat','scaled_arnoldi_failure_evidence.mat','scaled_arnoldi_failure_evidence.mat'};studies=struct([]);
for j=1:3
    evidenceFile=fullfile(val,folders{j},files{j});s=load(evidenceFile,'failure');f=s.failure;clear s;
    if j==1,context=fullfile(val,folders{j},'single_ilu_failure_context.mat');else,context=fullfile(val,folders{j},'input_binding.mat');end
    b=load(context,'sourceBinding');assert(strcmp(b.sourceBinding.sourceFileSHA256,sourceSHA));
    assert(~f.usableDirectionClaimed&&f.orthogonalityThreshold==1e-8&&f.modelRelationThreshold==1e-8);
    item=struct('label',names{j},'caseID',folders{j},'linearSystem','C=-Jgraph', ...
        'arnoldiSteps',f.arnoldiSteps,'CApplications',f.modelApplications,'ILUApplications',f.preconditionerApplications, ...
        'innerRelativeCResidual',norm(f.actualRemaining)/norm(f.b),'basisOrthogonalityError',f.orthogonalityError, ...
        'modelRelationRelativeDifference',f.modelRelationRelativeDifference,'orthogonalityThreshold',1e-8,'modelRelationThreshold',1e-8, ...
        'candidateUsable',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0, ...
        'evidenceFile',evidenceFile,'evidenceSHA256',rocketStreamingFileSHA256(evidenceFile));
    if isempty(studies),studies=item;else,studies(end+1)=item;end %#ok<AGROW>
    clear f;
end
report=struct('status','source-bound-rejected-internal-linear-studies-for-scientific-display', ...
    'sourceCase','krylov17','sourceEvidenceSHA256',sourceSHA,'studies',studies,'physicalStateChanged',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0,'newActualNSEvaluations',0,'browserVisualOrWebGLRuntimeVerified',false);
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
end
