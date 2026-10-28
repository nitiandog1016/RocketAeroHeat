function report=testRocketGraphActionSplit(folder)
%TESTROCKETGRAPHACTIONSPLIT 在已验证小域检查分项恒等式、原完整作用及原NS中心差分。
% 模块一：读取现有32CV同Q；不重新构造160列Jacobian或任何粗空间像。
% 模块二：确定性组合方向和剩余方向均对照保存的完整Jacobian及装配图矩阵。
% 模块三：另作两组原NS中心差分，并保存通量/梯度原始证据供独立语言复算。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'),fullfile(root,'tests'));setupCutCellFlow;
assert(~exist(folder,'dir'),'GraphSplitTest:Overwrite','不得覆盖分项单位校核。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','recycledkernel01','recycled_kernel_evidence.mat');
s=load(sourceFile,'Q','cut','gas','boundary','scale','unitA');Q=s.Q;cut=s.cut;gas=s.gas;boundary=s.boundary;scale=s.scale;unitA=s.unitA;clear s;
[rate,baseInfo]=cutCellNavierStokesResidual(Q,cut,gas,boundary);cache=buildCutCellSelectedDirectionCache(Q,baseInfo,cut,gas,boundary,scale);
operator=buildRocketShellGraphPreconditioner(Q,cut,gas,boundary,scale,1e-4,1);h=1e-6;errors=zeros(2,3);evidence=cell(2,1);
for j=1:2
    if j==1,directionU=reshape(sin((1:numel(Q))*.31),size(Q));else,directionU=rate./scale;end
    directionU=directionU/max(abs(directionU),[],'all');
    [full,graph,detail]=splitRocketGraphDirection(cache,operator,directionU);model=applyCutCellSelectedDirection(cache,directionU);
    dense=reshape(-unitA*directionU(:),size(Q));plus=cutCellNavierStokesResidual(Q+h*directionU.*scale,cut,gas,boundary);
    minus=cutCellNavierStokesResidual(Q-h*directionU.*scale,cut,gas,boundary);measured=(plus-minus)./scale/(2*h);
    errors(j,:)=[norm(detail.fullTotal-model,'fro')/max(1,norm(model,'fro')), ...
        norm(detail.fullTotal-dense,'fro')/max(1,norm(dense,'fro')),norm(detail.fullTotal-measured,'fro')/max(1,norm(measured,'fro'))];
    assert(max(errors(j,1:2))<2e-13&&errors(j,3)<2e-6,'GraphSplitTest:NS','分项或原NS方向核对失败。');
    evidence{j}=struct('directionU',directionU,'full',full,'graph',graph,'detail',detail,'plus',plus,'minus',minus);
end
report=struct('status','passed-small-current-J-and-graph-component-sums-dense-reference-and-original-NS-directions', ...
    'cellCount',size(Q,1),'directionCount',2,'directionErrors',errors,'completeSelectedActionApplications',4, ...
    'newCurrentJacobianColumnApplications',0,'actualOriginalNSEvaluations',5,'newFlowAcceptances',0,'physicalSolution',false, ...
    'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'graph_split_unit_evidence.mat'),'sourceFile','Q','cut','gas','boundary','scale','unitA','baseInfo','rate','operator','evidence','h','report','-v7');
fid=fopen(fullfile(folder,'graph_split_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
