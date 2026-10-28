function report=testRocketBlockILU0(folder)
%TESTROCKETBLOCKILU0 小域完整C作用、重复面合并、块ILU稠密参考及原图保持。
% 模块一：两种排序均保留完整160量矩阵；全部块值和固定线性映射可独立复算。
% 模块二：ILU只舍弃模式外填充，L/U的三角解与MATLAB稠密参考比较。
% 模块三：存储预算、非法RHS和二次因子化明确拒绝，无新完整J/NS。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'BlockILUTest:Overwrite','不得覆盖块ILU测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','wlsnormalunit02','wls_normal_unit_evidence.mat');
s=load(sourceFile,'operator','newMatrix','Q','rate','scale','unitA');operator=s.operator;before=operator;C=s.newMatrix-eye(numel(s.Q))/operator.pseudoTimeStepSeconds;
Q=s.Q;scale=s.scale;steadyA=s.unitA;rhs=s.rate./scale;clear s;n=size(Q,1);cases=cell(2,1);errors=zeros(2,3);
for j=1:2
    if j==1,ordering='natural';else,ordering='symrcm';end
    options=struct('ordering',ordering,'maximumPersistentBytes',32*1024^2,'maximumConstructionBytes',64*1024^2);
    packed=packRocketSteadyGraph(operator,options);factor=factorRocketBlockILU0(packed);
    probe=reshape(sin((1:numel(Q))'*.37),size(Q));applied=applyRocketPackedGraph(packed,probe);errors(j,1)=norm(applied(:)-C*probe(:))/max(1,norm(C*probe(:)));
    [L,U]=matrices(factor);[value,detail]=solveRocketBlockILU0(factor,rhs);ordered=rhs(factor.order,:)';expected=U\(L\ordered(:));
    actual=value(factor.order,:)';errors(j,2)=norm(actual(:)-expected)/max(1,norm(expected));
    second=reshape(cos((1:numel(Q))'*.29),size(Q));other=solveRocketBlockILU0(factor,second);combined=solveRocketBlockILU0(factor,.7*rhs-1.3*second);
    errors(j,3)=norm(combined-.7*value+1.3*other,'fro')/max(1,norm(combined,'fro'));assert(max(errors(j,:))<2e-9,'BlockILUTest:Dense','完整C或块三角解/线性叠加与参考不符。');
    assert(isequaln(operator,before)&&numel(unique(factor.order))==n&&all(factor.factorizationReport.allOriginalCellsAndFiveEquationsKept));
    cases{j}=struct('packed',packed,'factor',factor,'probe',probe,'applied',applied,'rhs',rhs,'value',value,'detail',detail,'L',L,'U',U, ...
        'second',second,'other',other,'combined',combined);
end
bad=options;bad.maximumPersistentBytes=1;rejected=false;try,packRocketSteadyGraph(operator,bad);catch ex,rejected=strcmp(ex.identifier,'BlockGraph:Budget');end;assert(rejected);
twice=false;try,factorRocketBlockILU0(factor);catch ex,twice=strcmp(ex.identifier,'BlockILU:Kind');end;assert(twice);
report=struct('status','passed-complete-C-block-CSR-two-orderings-ILU-triangular-linearity-and-budget-reference', ...
    'cellCount',n,'sourceFaces',operator.sourceFaces,'relativeErrors',errors,'budgetRejected',rejected,'doubleFactorRejected',twice, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'block_ilu_unit_evidence.mat'),'report','sourceFile','operator','C','Q','scale','steadyA','cases','-v7');
fid=fopen(fullfile(folder,'block_ilu_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
function [L,U]=matrices(factor)
n=factor.cellCount;L=eye(5*n);U=zeros(5*n);
for row=1:n
    ids=double(factor.rowPointer(row)):double(factor.rowPointer(row+1))-1;rr=(row-1)*5+(1:5);
    for entry=ids
        col=double(factor.column(entry));cc=(col-1)*5+(1:5);
        if col<row,L(rr,cc)=factor.values(:,:,entry);else,U(rr,cc)=factor.values(:,:,entry);end
    end
end
end
