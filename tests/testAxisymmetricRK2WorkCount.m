function report=testAxisymmetricRK2WorkCount(folder)
%TESTAXISYMMETRICRK2WORKCOUNT 用实际分析器核对RK2残差计数，EOS检查不混入。
% 模块一：小型非均匀周期环隙短推进，不重跑历史长稳态。
% 模块二：分析器实际入口次数与推进器自记次数及无拒绝1+2N独立对照。
root=fileparts(fileparts(mfilename('fullpath')));addpath(fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','rk2_cost');end
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
[X,R]=ndgrid([0 .02 .05 .08],[.03 .035 .045 .06]);
mesh=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));n=mesh.cellCount;
t=(1:n)';Q=axisymmetricFVConservative([1+.01*sin(t),-.4+.02*cos(t),.01*sin(.7*t),300+.1*cos(t)],gas);
per=struct('type','periodic');wall=struct('type','wall','temperature',300,'velocity',[0 0]);
boundary=struct('iMin',per,'iMax',per,'jMin',wall,'jMax',wall);
options=struct('maxSteps',7,'convectiveCFL',.1,'diffusiveCFL',.1);
profile clear;profile on;cleanup=onCleanup(@()profile('off'));
result=advanceAxisymmetricNavierStokes(Q,mesh,gas,boundary,options);
profile off;observed=profile('info');clear cleanup;
table=observed.FunctionTable;entry=table(strcmp({table.FunctionName},'axisymmetricNavierStokesResidual'));
assert(numel(entry)==1&&result.acceptedSteps==7&&result.rejectedAttempts==0&& ...
    entry.NumCalls==result.residualEvaluations&&result.residualEvaluations==15&& ...
    isequal(result.history.residualEvaluations,1:2:15),'AxisFVTest:WorkCount','实际空间算子次数与自记/1+2N不一致。');
assert(result.globalBalanceRelative<1e-10&&isfield(result,'time')&&~isfield(result,'pseudoTime'));
report=struct('status','passed-profiler-observed-RK2-residual-work-count', ...
    'acceptedSteps',7,'profileObservedResidualCalls',entry.NumCalls,'solverRecordedResidualCalls',result.residualEvaluations, ...
    'noRejectFormula','1+2*N','physicalSolution',false,'newProductionSolvedAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'rk2_cost_evidence.mat'),'result','mesh','gas','boundary','report','observed');
fid=fopen(fullfile(folder,'rk2_cost_summary.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report));fclose(fid);
fprintf('Actual profiler RK2 cost passed: %d calls for 7 accepted steps; EOS excluded.\n',entry.NumCalls);
end
