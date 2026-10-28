function report=testRocketRowScaledMomentum(folder)
%TESTROCKETROWSCALEDMOMENTUM 用独立稠密R*A与P*R^{-1}核对两种排序和残差身份。
% 模块一：原单位解必须与原A直接解一致，保存全部物理Z检查尺度/预条件回映。
% 模块二：原A目标与尺度目标分开，零右端/预算/零行不能制造伪合格状态。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'internalFlow'));setupCutCellFlow;
assert(~exist(folder,'dir'),'RowMomentumTest:Overwrite','不得覆盖行尺度测试。');mkdir(folder);
parentFile=fullfile(root,'outputs','v83_validation','momentumlineunit05','momentum_lines_unit_evidence.mat');s=load(parentFile,'cases');cases=cell(2,1);
options=struct('maximumIterations',96,'relativeTolerance',1e-10,'maximumBasisBytes',64*1024^2,'saveBasisEvidence',true);
for j=1:2
    c=s.cases{j};matrix=c.matrix;factor=c.factor;rhs=c.rhs;n=matrix.cellCount;[value,detail,rowScale,image]=solveRocketRowScaledMomentumKrylov(matrix,factor,rhs,options);
    expectedScale=1./max(abs(c.A),[],2);assert(norm(rowScale(:)-expectedScale,inf)<1e-13*norm(expectedScale,inf));
    B=rowScale(:).*c.A;P=c.P./rowScale(:)';direct=c.A\rhs(:);
    assert(detail.rowScaledTargetMet&&norm(value(:)-direct)/max(1,norm(direct))<1e-8&&detail.relativeAResidual<1e-8);
    assert(norm(image(:)-c.A*value(:))/max(1,norm(image(:)))<1e-11);
    assert(detail.originalUnitTargetMet==(detail.relativeAResidual<=options.relativeTolerance)&&detail.rowScaledTargetMet==(detail.krylov.flag==0));
    cycle=detail.krylov.cycles{1};assert(norm(B*cycle.Z-cycle.V*cycle.smallH,'fro')/max(1,norm(B*cycle.Z,'fro'))<1e-8);
    for k=1:size(cycle.Z,2),expected=P*cycle.V(:,k)/cycle.imageNormalizationScales(k);assert(norm(expected-cycle.Z(:,k))/max(1,norm(expected))<1e-10);end
    compact=compactRocketMomentumKrylovReport(detail);assert(compact.relativeAResidual==detail.relativeAResidual&&compact.relativeRowScaledResidual==detail.relativeRowScaledResidual&&numel(jsonencode(compact))<128*1024);
    [zero,z]=solveRocketRowScaledMomentumKrylov(matrix,factor,zeros(n,3),options);assert(all(zero==0,'all')&&z.solverAApplications==0&&z.originalUnitDiagnosticAApplications==1&&z.lineInverseApplications==0&&z.originalUnitTargetMet);
    bad=options;bad.maximumBasisBytes=1;rejected=false;try,solveRocketRowScaledMomentumKrylov(matrix,factor,rhs,bad);catch ex,rejected=strcmp(ex.identifier,'RecycledGMRES:Budget');end;assert(rejected);
    badMatrix=matrix;ids=double(matrix.rowPointer(1)):double(matrix.rowPointer(2))-1;badMatrix.values(1,:,ids)=0;zeroRowRejected=false;try,rocketMomentumRowScales(badMatrix);catch ex,zeroRowRejected=strcmp(ex.identifier,'MomentumScale:ZeroRow');end;assert(zeroRowRejected);
    cases{j}=struct('A',c.A,'P',c.P,'B',B,'PRinverse',P,'matrix',matrix,'factor',factor,'rhs',rhs,'value',value,'image',image,'rowScale',rowScale,'detail',detail,'direct',direct);
end
report=struct('status','passed-two-orderings-equivalent-row-system-right-preconditioner-original-units-and-dual-residuals','cases',2,'budgetRejected',true,'zeroRowRejected',true,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0,'engineProcessID',feature('getpid'));
save(fullfile(folder,'row_momentum_unit_evidence.mat'),'report','cases','parentFile','-v7');fid=fopen(fullfile(folder,'row_momentum_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));fprintf('%s\n',jsonencode(report));
end
