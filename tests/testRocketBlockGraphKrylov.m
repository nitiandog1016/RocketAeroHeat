function report=testRocketBlockGraphKrylov(folder)
%TESTROCKETBLOCKGRAPHKRYLOV 在完整小域C上核对块ILU内部Krylov和独立稠密解。
% 模块一：使用上一单位测试的已验证因子，不重新构造当前Jacobian。
% 模块二：原B及五带时间对角保持；内C物理方向和剩余可供独立复算。
% 模块三：零右端/内存拒绝不生成流场解，本测试新增J/NS均为零。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'BlockKrylovTest:Overwrite','不得覆盖块图内部求解测试。');mkdir(folder);
sourceFile=fullfile(root,'outputs','v83_validation','blockiluunit01','block_ilu_unit_evidence.mat');
s=load(sourceFile);operator=s.operator;factor=s.cases{2}.factor;C=s.C;Q=s.Q;scale=s.scale;steadyA=s.steadyA;rhs=s.cases{2}.rhs;clear s;before=operator;
options=struct('maximumIterations',64,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
[value,detail]=solveRocketBlockGraphKrylov(operator,factor,rhs,options);direct=C\rhs(:);relative=norm(C*value(:)-rhs(:))/norm(rhs(:));
assert(detail.krylov.flag==0&&relative<1.1e-10&&norm(value(:)-direct)/max(1,norm(direct))<2e-8&&isequaln(operator,before));
[zero,zeroReport]=solveRocketBlockGraphKrylov(operator,factor,zeros(size(rhs)),options);assert(all(zero==0,'all')&&zeroReport.CApplications==0);
bad=options;bad.maximumBasisBytes=1;rejected=false;try,solveRocketBlockGraphKrylov(operator,factor,rhs,bad);catch ex,rejected=strcmp(ex.identifier,'RecycledGMRES:Budget');end;assert(rejected);
report=struct('status','passed-full-small-C-block-ILU-inner-Krylov-dense-reference-and-budget', ...
    'cellCount',size(Q,1),'relativeCResidual',relative,'budgetRejected',rejected,'newSelectedJApplications',0,'newActualNSEvaluations',0, ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'block_krylov_unit_evidence.mat'),'report','sourceFile','operator','factor','C','Q','scale','steadyA','rhs','options','value','detail','direct','zero','zeroReport','-v7');
fid=fopen(fullfile(folder,'block_krylov_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
