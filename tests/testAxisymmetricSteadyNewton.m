function report=testAxisymmetricSteadyNewton(folder)
%TESTAXISYMMETRICSTEADYNEWTON 用解析均匀开放层流解检查直接定常Newton，不作实验认证。
% 模块一：同速等温墙及完整Riemann储库，常压几何源平衡已知；从扰动场实际求解。
% 模块二：检查完整四量残差/边界账、解析解误差和闭域/非法输入拒绝。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v75_validation','steady_root_unit');end
assert(~exist(folder,'dir'),'RootUnit:Overwrite','不能覆盖测试。');
gas=axisymmetricFVGas(struct('gamma',1.4,'R',1,'mu',.05,'Pr',.72,'temperatureRange',[200 800]));
[x,r]=ndgrid([0 .02 .05],[.03 .05 .07]);mesh=buildAxisymmetricFVGeometry(x+.07*(r-.03),r);
n=mesh.cellCount;t=(1:n)';exactP=repmat([1 -2 0 500],n,1);P=exactP;
P(:,1)=P(:,1).*(1+.001*sin(t));P(:,2)=P(:,2)+.001*cos(t);P(:,4)=P(:,4)+.01*sin(t);
Q=axisymmetricFVConservative(P,gas);trace=struct('type','riemannReservoir','primitive',[1 -2 0 500]);
wall=struct('type','wall','temperature',500,'velocity',[-2 0]);
bc=struct('iMin',trace,'iMax',trace,'jMin',wall,'jMax',wall,'pressureSourceQuadrature','cell-constant','interiorConvectiveFlux','hlle-einfeldt');
options=struct('residualScale',[1 20 20 2000]);solution=solveAxisymmetricSteadyNewton(Q,mesh,gas,bc,options);
assert(solution.spatialRootConverged&&solution.acceptedNewtonIterations>=1&&solution.residualNorm<=2e-6&&solution.boundaryFluxBalanceRelative<=5e-10);
actual=axisymmetricFVState(solution.Q,gas);error=max(abs(actual(:,1:4)-exactP)./max(abs(exactP),1),[],'all');assert(error<2e-7);
refusals=0;closed=bc;closed.iMin=wall;closed.iMax=wall;reject(@()solveAxisymmetricSteadyNewton(Q,mesh,gas,closed,options),'AxisRoot:ClosedDomain');
bad=options;bad.unknown=1;reject(@()solveAxisymmetricSteadyNewton(Q,mesh,gas,bc,bad),'AxisRoot:Options');
report=struct('status','passed-uniform-open-four-quantity-direct-spatial-Newton-analytic-check', ...
    'matlabEngineProcessID',feature('getpid'),'newtonIterations',solution.acceptedNewtonIterations, ...
    'residual',solution.residualNorm,'boundaryBalanceRelative',solution.boundaryFluxBalanceRelative, ...
    'analyticPrimitiveRelativeError',error,'refusals',refusals,'physicalSolution',false,'experimentalValidation',false);
mkdir(folder);save(fullfile(folder,'steady_root_unit_evidence.mat'),'report','solution','exactP','mesh','gas','bc','-v7');
fid=fopen(fullfile(folder,'steady_root_unit_summary.json'),'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear guard;fprintf('%s\n',jsonencode(report));
    function reject(action,id)
        try,action();catch exception,assert(strcmp(exception.identifier,id));refusals=refusals+1;return;end;error('RootUnit:Refusal','非法输入未拒绝。');
    end
end
