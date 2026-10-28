function report=testRocketSteadyGraphPreconditioner(folder)
%TESTROCKETSTEADYGRAPHPRECONDITIONER 核对整体五量去正则、同线新因子及稠密内求解。
% 模块一：C与旧B减I/tau逐值一致，原B及旧因子保持原字节/数值。
% 模块二：新固定扫描逐色对照稠密解，不能沿用旧因子或忽略第二邻居。
% 模块三：内C求解仅为模型研究，不提交非线性状态或物理资格。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'SteadyGraphTest:Overwrite','不得覆盖稳态图测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','wlsnormalunit02','wls_normal_unit_evidence.mat');
s=load(sourceFile,'operator','newMatrix','Q','rate','scale','unitA');original=s.operator;before=original;B=s.newMatrix;Q=s.Q;scale=s.scale;steadyA=s.unitA;rhs=s.rate./scale;clear s;
operator=buildRocketSteadyGraphPreconditioner(original);assert(isequaln(original,before),'SteadyGraphTest:Mutation','新分解修改了旧算子。');
n=size(Q,1);nf=operator.sourceFaces;C=B-eye(numel(Q))/operator.pseudoTimeStepSeconds;
[value,detail]=solveRocketWlsNormalPreconditioner(operator,rhs);expected=zeros(numel(Q),1);groups=operator.colorGroups;nc=numel(groups);order=[1:nc nc-1:-1:1];
for color=order
    cells=groups{color}.cellIDs;ids=reshape(cells+(0:4)*n,[],1);other=setdiff((1:numel(Q))',ids);
    expected(ids)=C(ids,ids)\(rhs(ids)-C(ids,other)*expected(other));
end
scanError=norm(value(:)-expected)/max(1,norm(expected));applied=applyRocketWlsNormalOperator(operator,value);
actionError=norm(applied(:)-C*value(:))/max(1,norm(applied(:)));assert(max(scanError,actionError)<2e-10,'SteadyGraphTest:Dense','C五带作用/扫描与稠密参考不符。');
options=struct('maximumIterations',64,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
[inner,innerReport]=solveRocketInnerGraphKrylov(operator,rhs,options);direct=C\rhs(:);innerError=norm(C*inner(:)-rhs(:))/norm(rhs(:));
assert(innerReport.krylov.flag==0&&innerError<1.1e-10&&norm(inner(:)-direct)/max(1,norm(direct))<2e-8,'SteadyGraphTest:Inner','内部C求解与稠密参考不符。');
twice=false;try,buildRocketSteadyGraphPreconditioner(operator);catch ex,twice=strcmp(ex.identifier,'SteadyGraph:AlreadyUnshifted');end;assert(twice);
report=struct('status','passed-all-five-unshifted-C-new-factor-dense-scan-inner-reference-and-old-B-preservation', ...
    'cellCount',n,'sourceFaces',nf,'denseScanRelativeDifference',scanError,'denseActionRelativeDifference',actionError, ...
    'innerCRelativeResidual',innerError,'doubleUnshiftRejected',twice,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'steady_graph_unit_evidence.mat'),'report','sourceFile','original','operator','B','C','steadyA','Q','scale','rhs', ...
    'value','detail','expected','options','inner','innerReport','direct','-v7');
fid=fopen(fullfile(folder,'steady_graph_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
