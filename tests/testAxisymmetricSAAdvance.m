function report=testAxisymmetricSAAdvance(folder)
%TESTAXISYMMETRICSAADVANCE 小闭域真五量BE一步、重启和独立账核对。
% 模块一：同一Newton实际更新五量，W不固定，质量仅由原闭域约束保持。
% 模块二：从初末完整残差独立复算BE及5项边界/源账，伪时间不作物理时间。
% 模块三：两步连续与保存重启相同，禁止裁N或以预算终态冒充稳态。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v65_validation','sa_advance');end
if ~exist(folder,'dir'),mkdir(folder);end
[X,R]=ndgrid([0 .01 .02],[.1 .11 .12]);mesh=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));
gas=createTmrAxisymmetricAirGas;rho=.1*ones(mesh.cellCount,1);N=2e-4+1e-3*(mesh.volumeCentroid(:,2)-.1);
Q4=axisymmetricFVConservative([rho,-12*ones(size(rho)),zeros(size(rho)),300*ones(size(rho))],gas);Q=[Q4 rho.*N];
boundary=struct('iMin',struct('type','periodic'),'iMax',struct('type','periodic'), ...
    'jMin',struct('type','wall','temperature',290,'velocity',[-10 0]), ...
    'jMax',struct('type','wall','temperature',310,'velocity',[-12 0]));
scale=[.1 12 12 max(Q(:,4)) max(Q(:,5))];dt=1e-6;
opts=struct('maxSteps',1,'initialPseudoCFL',1e3,'maximumPseudoCFL',1e3,'maximumPseudoStep',dt, ...
    'residualScale',scale,'jacobianRefreshEvery',1,'maxNewtonIterations',18);
one=advanceAxisymmetricSteady(Q,mesh,gas,boundary,opts);
assert(one.acceptedSteps==1&&one.closedMassConstraintActive&&~one.steadyConverged&&one.history.dt(end)==dt);
[rate,info]=axisymmetricRANSResidual(one.Q,mesh,gas,boundary);
be=(one.Q-Q)/dt-rate;normalized=max(abs(be)./scale,[],'all');
stateDefect=max(abs((one.Q-Q-dt*rate)./scale),[],'all');
assert(stateDefect<=one.history.newtonTarget(end)&&isequaln(info,one.finalInfo));
ledger=sum((one.Q-Q).*mesh.volumePerRadian,1)+dt*(info.netBoundaryFlux-info.sourceIntegral);
total=sum(mesh.volumePerRadian)*scale;relative=abs(ledger)./total;
assert(all(relative<=one.options.conservationTolerance)&&numel(ledger)==5&& ...
    abs((sum(one.Q(:,1).*mesh.volumePerRadian)-sum(Q(:,1).*mesh.volumePerRadian))/sum(Q(:,1).*mesh.volumePerRadian))<1e-12);
assert(all(one.Q(:,5)>0)&&any(one.Q(:,5)~=Q(:,5))&&any(one.Q(:,2:4)~=Q(:,2:4),'all'));
save(fullfile(folder,'sa_be_step.mat'),'Q','mesh','gas','boundary','one','ledger','relative','be');
loaded=load(fullfile(folder,'sa_be_step.mat'));restart=advanceAxisymmetricSteady(loaded.one.Q,mesh,gas,boundary,opts);
twoOpts=opts;twoOpts.maxSteps=2;two=advanceAxisymmetricSteady(Q,mesh,gas,boundary,twoOpts);
restartDifference=restart.Q-two.Q;restartRelativeByQuantity=max(abs(restartDifference)./scale,[],1);
save(fullfile(folder,'sa_restart_comparison.mat'),'restart','two','restartDifference','restartRelativeByQuantity');
% 重启重新形成U=Q/scale及质量锚点，会带来浮点往返；保留逐量原差而非要求字节相同。
assert(all(restartRelativeByQuantity<=64*eps)&&restart.acceptedSteps==1&&two.acceptedSteps==2&&~two.steadyConverged, ...
    'AxisSATest:Restart','重启状态差超过64eps的逐量参考尺度。');
report=struct('status','passed-actual-five-variable-BE-step-ledger-and-restart', ...
    'stateBEInfinityNorm',stateDefect,'rateBEInfinityNorm',normalized,'ledgerRaw',ledger,'ledgerRelative',relative, ...
    'actualResidualEvaluations',one.residualEvaluations,'jacobianBuilds',one.jacobianBuilds,'stateWasClipped',false, ...
    'restartDifference',restartDifference,'restartRelativeByQuantity',restartRelativeByQuantity, ...
    'fifthVariableChanged',true,'allFiveQuantitiesCoupled',true,'pseudoTimeIsPhysical',false, ...
    'steadyConverged',false,'physicalSolution',false,'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'sa_advance_evidence.mat'),'report','Q','mesh','gas','boundary','one','restart','two');
fid=fopen(fullfile(folder,'sa_advance_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('SA actual coupled BE/restart passed: %d calls, five-quantity relative ledger %.3g; no steady/physical claim.\n', ...
    one.residualEvaluations,max(relative));
end
