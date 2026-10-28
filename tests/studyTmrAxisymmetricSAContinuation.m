function report=studyTmrAxisymmetricSAContinuation(folder,source,maxSteps)
%STUDYTMRAXISYMMETRICSACONTINUATION 从真实五量检查点继续求稳态，原父证据只读。
% 模块一：核对完整原方程、原尺度、边界、几何和初场SHA，禁止复用旧四量墙热。
% 模块二：只增加工作预算/CFL，原2e-6连续三步、非线性和五量守恒门槛不变。
% 模块三：本段与跨段账分别保留；伪时间不等于飞行时间，未稳态不作物理解。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'),fullfile(root,'internalFlow'));
if nargin<1,folder=fullfile(root,'outputs','v65_validation','tmr_sa_continue1');end
if nargin<2,source=fullfile(root,'outputs','v65_validation','tmr_sa_final','tmr_sa_evidence.mat');end
if nargin<3,maxSteps=50;end
validateattributes(maxSteps,{'numeric'},{'real','finite','scalar','integer','positive'});
sourceSHA=fileSHA(source);parent=load(source,'implicit','mesh','gas','boundary','report');
mesh=parent.mesh;gas=parent.gas;boundary=parent.boundary;old=parent.implicit;
assert(size(old.Q,2)==5&&~old.physicalSolution&&~old.eligibleForSurfaceCoverage&& ...
    ~old.pseudoTimeIsPhysical&&old.physicalTimeAdvanced==0,'AxisSA:Continuation','只接受真实五量研究检查点。');
[before,initialInfo]=axisymmetricRANSResidual(old.Q,mesh,gas,boundary);
assert(isequaln(initialInfo,old.finalInfo)&&max(abs(before)./old.residualScale,[],'all')==old.residualNorm, ...
    'AxisSA:ContinuationReplay','当前方程不能重放父检查点。');
assert(old.options.steadyResidualTolerance==2e-6&&old.options.steadyRequiredSteps==3);
opts=old.options;opts.maxSteps=maxSteps;opts.initialPseudoCFL=8;opts.maximumPseudoCFL=1e6;
opts.maxNewtonIterations=24;opts.jacobianRefreshEvery=1;
assert(~exist(fullfile(folder,'tmr_sa_evidence.mat'),'file'),'AxisSA:Overwrite','不能覆盖父段或已有证据。');
if ~exist(folder,'dir'),mkdir(folder);end
clock=tic;implicit=advanceAxisymmetricSteady(old.Q,mesh,gas,boundary,opts);elapsed=toc(clock);
[rate,info]=axisymmetricRANSResidual(implicit.Q,mesh,gas,boundary);
assert(isequaln(info,implicit.finalInfo)&&max(abs(rate)./implicit.residualScale,[],'all')==implicit.residualNorm&& ...
    strcmp(sourceSHA,fileSHA(source)));
if isfield(parent.report,'cumulative')
    cumulative=parent.report.cumulative;
else
    cumulative=struct('initialConservedIntegral',old.initialConservedIntegral, ...
        'boundaryIntegral',zeros(1,5),'sourceIntegral',zeros(1,5),'acceptedSteps',0, ...
        'rejectedAttempts',0,'solverResidualEvaluations',0,'pseudoTime',0);
    cumulative.boundaryIntegral=old.accumulatedBoundaryIntegral;cumulative.sourceIntegral=old.accumulatedSourceIntegral;
    cumulative.acceptedSteps=old.acceptedSteps;cumulative.rejectedAttempts=old.rejectedAttempts;
    cumulative.solverResidualEvaluations=old.residualEvaluations;cumulative.pseudoTime=old.pseudoTime;
end
cumulative.boundaryIntegral=cumulative.boundaryIntegral+implicit.accumulatedBoundaryIntegral;
cumulative.sourceIntegral=cumulative.sourceIntegral+implicit.accumulatedSourceIntegral;
cumulative.acceptedSteps=cumulative.acceptedSteps+implicit.acceptedSteps;
cumulative.rejectedAttempts=cumulative.rejectedAttempts+implicit.rejectedAttempts;
cumulative.solverResidualEvaluations=cumulative.solverResidualEvaluations+implicit.residualEvaluations;
cumulative.pseudoTime=cumulative.pseudoTime+implicit.pseudoTime;
cumulative.rawFiveQuantityDefect=implicit.finalConservedIntegral-cumulative.initialConservedIntegral+ ...
    cumulative.boundaryIntegral-cumulative.sourceIntegral;
cumulative.fiveQuantityRelativeDefect=abs(cumulative.rawFiveQuantityDefect)./(sum(mesh.volumePerRadian)*implicit.residualScale);
assert(all(cumulative.fiveQuantityRelativeDefect<=opts.conservationTolerance));
report=struct('status','actual-coupled-SA-continuation-not-experimental-validation','modelID',info.sa.modelID, ...
    'parentCheckpoint',source,'parentCheckpointSHA256',sourceSHA,'cumulative',cumulative, ...
    'initialOriginalResidual',old.residualNorm,'finalOriginalResidual',implicit.residualNorm, ...
    'steadyConverged',implicit.steadyConverged,'stopReason',implicit.stopReason,'acceptedSteps',implicit.acceptedSteps, ...
    'rejectedAttempts',implicit.rejectedAttempts,'solverResidualEvaluations',implicit.residualEvaluations, ...
    'diagnosticResidualEvaluations',2,'elapsedSeconds',elapsed,'globalBEFiveQuantityBalanceRelative',implicit.globalBalanceRelative, ...
    'finalNuTildeRange',[min(info.sa.nuTilde),max(info.sa.nuTilde)], ...
    'wallSignedHeatRangeWm2',[min(info.wall.heatFluxIntoWall),max(info.wall.heatFluxIntoWall)], ...
    'pseudoTimeIsPhysical',false,'physicalTimeAdvanced',0,'turbulenceTransportImplemented',true, ...
    'turbulenceValidated',false,'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'experimentalValidation',false,'newProductionSolvedAreaM2',0);
save(fullfile(folder,'tmr_sa_evidence.mat'),'report','implicit','mesh','gas','boundary','sourceSHA','initialInfo');
fid=fopen(fullfile(folder,'tmr_sa_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('TMR SA continuation: %d steps residual %.9g -> %.9g, steady=%d, cumulative ledger %.3g.\n', ...
    implicit.acceptedSteps,old.residualNorm,implicit.residualNorm,implicit.steadyConverged,max(cumulative.fiveQuantityRelativeDefect));
end
function digest=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
