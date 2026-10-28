function report=testRocketInnerGraphKrylov(folder)
%TESTROCKETINNERGRAPHKRYLOV 用已有小域完整B核对内层求逆、预算和当前稳态作用。
% 模块一：稠密B直接解与物理Arnoldi方向/像共同验证，不调用新的原NS。
% 模块二：构造一次扫描的完整线性映射，分开记录B与稳态A=-J的性质。
% 模块三：零右端和存储拒绝保持明确；本模块不进行非线性接受。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'InnerGraphTest:Overwrite','不得覆盖内部图求逆测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','wlsnormalunit02','wls_normal_unit_evidence.mat');
s=load(sourceFile,'operator','newMatrix','Q','rate','scale','unitA');operator=s.operator;B=s.newMatrix;Q=s.Q;scale=s.scale;steadyA=s.unitA;rate=s.rate;clear s;
n=numel(Q);P=zeros(n);for column=1:n,basis=zeros(size(Q));basis(column)=1;value=solveRocketWlsNormalPreconditioner(operator,basis);P(:,column)=value(:);end
options=struct('maximumIterations',32,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
inputs={reshape(sin((1:n)'*.29),size(Q)),rate./scale};evidence=cell(2,1);errors=zeros(2,2);
for j=1:2
    rhs=inputs{j};[value,detail]=solveRocketInnerGraphKrylov(operator,rhs,options);expected=B\rhs(:);
    errors(j,:)=[norm(B*value(:)-rhs(:))/norm(rhs(:)),norm(value(:)-expected)/max(1,norm(expected))];
    assert(detail.krylov.flag==0&&errors(j,1)<=1.1e-10&&errors(j,2)<2e-9,'InnerGraphTest:Dense','内部B求解不等于独立稠密参考。');
    assert(~detail.outerFixedLinearMapClaimed&&detail.requiresFlexibleOuterMethod&&detail.newActualNSEvaluations==0);
    evidence{j}=struct('rhs',rhs,'value',value,'directBReference',expected,'detail',detail);
end
[zero,zeroReport]=solveRocketInnerGraphKrylov(operator,zeros(size(Q)),options);assert(all(zero==0,'all')&&zeroReport.totalBApplicationsIncludingScanDiagnostics==0);
bad=options;bad.maximumBasisBytes=1;rejected=false;try,solveRocketInnerGraphKrylov(operator,inputs{1},bad);catch ex,rejected=strcmp(ex.identifier,'RecycledGMRES:Budget');end;assert(rejected);
report=struct('status','passed-small-dense-B-inner-GMRES-Arnoldi-budget-and-zero-RHS-reference', ...
    'cellCount',size(Q,1),'rightHandSides',2,'denseReferenceRelativeErrors',errors,'budgetRejected',rejected, ...
    'oneScanMatrixConstructionApplications',n,'oneScanConstructionDiagnosticBApplications',n, ...
    'newSelectedJApplications',0,'newActualNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'inner_graph_unit_evidence.mat'),'report','sourceFile','B','P','steadyA','Q','scale','rate','options','evidence','zero','zeroReport','-v7');
fid=fopen(fullfile(folder,'inner_graph_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
