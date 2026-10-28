function report=testRocketAxisymmetricControlledCheckpoint(folder,auditFolder)
%TESTROCKETAXISYMMETRICCONTROLLEDCHECKPOINT 独立逐端点复算真实完整箭体的五量BE链。
% 模块一：从冻结输入重建完整几何，检查源码字节和每个不可变端点SHA。
% 模块二：逐步重算真实初末完整RANS、原Newton缺陷与边界/源积分。
% 模块三：从所有真实端点重建历史、累计五量账和实际工作量，拒绝初猜冒充解。
% 模块四：只做只读审计，不推进原目录、不改变数值门槛或生产资格。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v69_validation','work_arithmetic');end
if nargin<2,auditFolder=fullfile(root,'outputs','v69_validation','controlled_checkpoint_audit');end
saved=load(fullfile(folder,'rocket_axisymmetric_evidence.mat'));
assert(all(isfield(saved,{'caseDefinition','implicit','initialInfo','durable','report'})), ...
    'RocketCheckpointTest:Payload','缺少真实工况/端点/链/报告。');
c=saved.caseDefinition;last=saved.implicit;d=saved.durable;
assert(isfield(c,'caseFactory'),'RocketCheckpointTest:Identity','受控链必须显式登记工况工厂。');
if isequal(c.caseFactory,'createRocketAxisymmetricWallCase')
    kind='rocket-complete-C-domain-SA-controlled-discretization-checkpoints-v1';
    fresh=createRocketAxisymmetricWallCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricCollarCase')
    kind='rocket-complete-C-domain-SA-collar-family-checkpoints-v1';
    fresh=createRocketAxisymmetricCollarCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricPressureCase')
    kind='rocket-complete-C-domain-SA-pressure-source-checkpoints-v1';
    fresh=createRocketAxisymmetricPressureCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricFluxCase')
    kind='rocket-complete-C-domain-SA-interior-flux-checkpoints-v1';
    fresh=createRocketAxisymmetricFluxCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricLocalCase')
    kind='rocket-complete-C-domain-SA-local-Jacobian-checkpoints-v1';
    fresh=createRocketAxisymmetricLocalCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricRootCase')
    kind='rocket-complete-C-domain-SA-spatial-root-confirmed-checkpoints-v1';
    fresh=createRocketAxisymmetricRootCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricSensitivityCase')
    kind='rocket-complete-C-domain-SA-sensitivity-root-confirmed-checkpoints-v1';
    fresh=createRocketAxisymmetricSensitivityCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricWallResolutionCase')
    kind='rocket-complete-C-domain-SA-wall-resolution-root-confirmed-checkpoints-v1';
    fresh=createRocketAxisymmetricWallResolutionCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricContactCase')
    kind='rocket-complete-C-domain-SA-contact-root-confirmed-checkpoints-v1';
    fresh=createRocketAxisymmetricContactCase(c.options);
elseif isequal(c.caseFactory,'createRocketAxisymmetricContactRefinedCase')
    kind='rocket-complete-C-domain-SA-contact-refined-root-confirmed-checkpoints-v1';
    fresh=createRocketAxisymmetricContactRefinedCase(c.options);
else
    error('RocketCheckpointTest:Identity','未知的受控工况工厂不能替代真实链。');
end
assert(strcmp(d.kind,kind)&& ...
    last.acceptedSteps>=1&&numel(d.steps)==last.acceptedSteps&& ...
    ~last.physicalSolution&&~last.eligibleForSurfaceCoverage&& ...
    ~last.pseudoTimeIsPhysical&&last.physicalTimeAdvanced==0, ...
    'RocketCheckpointTest:Identity','只接受已实际推进且未晋升物理解的真实链。');
for k=1:numel(d.sourceBindings)
    binding=d.sourceBindings(k);
    assert(strcmp(fileSHA(fullfile(root,binding.path)),binding.sha256), ...
        'RocketCheckpointTest:Source','端点来源源码已变化：%s。',binding.path);
end
assert(isequaln(fresh,c),'RocketCheckpointTest:Case','冻结工况不能从原输入逐值重建。');
mesh=c.mesh;gas=c.gas;boundary=c.boundary;V=mesh.volumePerRadian;
scale=c.residualScale;totalScale=sum(V)*scale;initial=sum(c.Q.*V,1);
Q=c.Q;[rate,info]=axisymmetricRANSResidual(Q,mesh,gas,boundary);
assert(isequaln(info,saved.initialInfo)&&isequal(last.initialQ,Q));
residuals=max(abs(rate)./scale,[],'all');pseudoTimes=0;balances=0;dtList=0;
newtonTargets=0;newtonDefects=0;stepDefects=zeros(1,5);
sourceIntegral=zeros(1,5);boundaryIntegral=zeros(1,5);pseudoTime=0;
evaluations=0;rejected=0;linearSolves=0;factorizations=0;jacobianBuilds=0;
solverContract=rocketAxisymmetricSolverContract(c);
nextCFL=solverContract.initialPseudoCFL;maxDefectRatio=0;maxStepBalance=0;steps=struct([]);
for k=1:last.acceptedSteps
    entry=d.steps(k);relative=sprintf('steps/step_%06d.mat',k);path=fullfile(folder,relative);
    assert(entry.number==k&&strcmp(entry.relativePath,relative)&&strcmp(fileSHA(path),entry.sha256), ...
        'RocketCheckpointTest:SHA','端点缺失、改字节、断号或被错误重排。');
    payload=load(path,'segment','stepContext');s=payload.segment;ctx=payload.stepContext;
    assert(ctx.number==k&&isequaln(ctx.sourceBindings,d.sourceBindings)&& ...
        ctx.matlabProcessID==entry.matlabProcessID&&isequal(s.initialQ,Q)&& ...
        s.acceptedSteps==1&&isequal(s.residualScale,scale)&& ...
        s.options.initialPseudoCFL==nextCFL&&~s.closedMassConstraintActive&& ...
        ~s.physicalSolution&&~s.eligibleForSurfaceCoverage&&~s.pseudoTimeIsPhysical&&s.physicalTimeAdvanced==0, ...
        'RocketCheckpointTest:Join','真实BE端点未与上一步、源码、实际CFL或开放域合同相接。');
    assert(s.options.steadyResidualTolerance==2e-6&&s.options.steadyRequiredSteps==3&& ...
        s.options.minimumSteadySteps==3&&s.options.conservationTolerance==5e-10&& ...
        s.options.nonlinearAbsoluteTolerance==1e-11&&s.options.nonlinearRelativeTolerance==1e-9&& ...
        s.options.differenceStep==solverContract.differenceStep&&s.options.maximumPseudoCFL==solverContract.maximumPseudoCFL&&strcmp(s.options.differenceScheme,'central')&& ...
        s.options.jacobianRefreshEvery==1&&s.options.maxNewtonIterations==24, ...
        'RocketCheckpointTest:Tolerances','原稳态/BE/守恒验收门槛改变。');
    [endRate,endInfo]=axisymmetricRANSResidual(s.Q,mesh,gas,boundary);
    assert(isequaln(endInfo,s.finalInfo),'RocketCheckpointTest:Replay','完整五量info不能逐值重放。');
    dt=s.history.dt(end);explicit=1/max(info.convectiveFrequency+info.diffusiveFrequency);
    assert(isfinite(dt)&&dt>0&&ctx.explicitStepAtInitialState==explicit&& ...
        isequal(s.accumulatedBoundaryIntegral,dt*endInfo.netBoundaryFlux)&& ...
        isequal(s.accumulatedSourceIntegral,dt*endInfo.sourceIntegral), ...
        'RocketCheckpointTest:Integrals','实际dt或边界/源积分不一致。');
    U0=Q./scale;U1=s.Q./scale;defect=U1-U0-dt*(endRate./scale);
    target=1e-11+1e-9*max(abs(dt*(rate./scale)),[],'all');
    roundoff=512*eps*max(abs(U0)+abs(U1)+abs(dt*(endRate./scale)),[],'all');
    actual=max(abs(defect),[],'all');
    assert(s.history.newtonTarget(end)==target&&actual<=target+roundoff&& ...
        abs(s.history.newtonDefect(end)-actual)<=roundoff, ...
        'RocketCheckpointTest:BE','端点真实BE缺陷没有通过原非线性目标。');
    final=sum(s.Q.*V,1);before=sum(Q.*V,1);
    stepDefect=final-before+s.accumulatedBoundaryIntegral-s.accumulatedSourceIntegral;
    stepBalance=max(abs(stepDefect)./totalScale);
    assert(stepBalance<=5e-10&&isequal(s.initialConservedIntegral,before)&& ...
        isequal(s.finalConservedIntegral,final), ...
        'RocketCheckpointTest:StepLedger','逐步五量账或保存积分不符。');
    sourceIntegral=sourceIntegral+s.accumulatedSourceIntegral;
    boundaryIntegral=boundaryIntegral+s.accumulatedBoundaryIntegral;
    balance=max(abs(final-initial+boundaryIntegral-sourceIntegral)./totalScale);
    assert(balance<=5e-10,'RocketCheckpointTest:ChainLedger','全链五量账不通过。');
    residual=max(abs(endRate)./scale,[],'all');assert(residual==s.residualNorm);
    nextCFL=min(s.options.maximumPseudoCFL,s.options.pseudoCFLGrowth*dt/explicit);
    assert(ctx.nextPseudoCFL==nextCFL);
    pseudoTime=pseudoTime+dt;Q=s.Q;rate=endRate;info=endInfo;
    residuals(end+1)=residual;pseudoTimes(end+1)=pseudoTime;balances(end+1)=balance;dtList(end+1)=dt; %#ok<AGROW>
    newtonTargets(end+1)=target;newtonDefects(end+1)=s.history.newtonDefect(end);stepDefects(end+1,:)=stepDefect; %#ok<AGROW>
    evaluations=evaluations+s.residualEvaluations;rejected=rejected+s.rejectedAttempts;
    linearSolves=linearSolves+s.linearSolves;factorizations=factorizations+s.factorizations;
    jacobianBuilds=jacobianBuilds+s.jacobianBuilds;
    maxDefectRatio=max(maxDefectRatio,actual/(target+roundoff));maxStepBalance=max(maxStepBalance,stepBalance);
    stepReport=struct('number',k,'processID',ctx.matlabProcessID,'residual',residual, ...
        'dt',dt,'beDefect',actual,'newtonTarget',target,'stepBalance',stepBalance,'chainBalance',balance);
    if isempty(steps),steps=stepReport;else,steps(k)=stepReport;end %#ok<AGROW>
end
assert(isequal(Q,last.Q)&&isequaln(info,last.finalInfo)&&isequal(last.history.residual,residuals)&& ...
    isequal(last.history.pseudoTime,pseudoTimes)&&isequal(last.history.dt,dtList)&& ...
    isequal(last.history.balance,balances)&&isequal(last.history.newtonTarget,newtonTargets)&& ...
    isequal(last.history.newtonDefect,newtonDefects)&&isequal(last.history.stepDefect,stepDefects)&& ...
    isequal(last.accumulatedSourceIntegral,sourceIntegral)&&isequal(last.accumulatedBoundaryIntegral,boundaryIntegral)&& ...
    last.globalBalanceRelative==balances(end)&&last.pseudoTime==pseudoTime&&d.nextPseudoCFL==nextCFL&& ...
    last.residualEvaluations==evaluations&&last.rejectedAttempts==rejected&&last.linearSolves==linearSolves&& ...
    last.factorizations==factorizations&&last.jacobianBuilds==jacobianBuilds, ...
    'RocketCheckpointTest:Aggregate','全链历史、成本、CFL或累计五量积分未从真实端点重建。');
expectedSteady=last.acceptedSteps>=3&&all(residuals(end-2:end)<=2e-6);
assert(last.steadyConverged==expectedSteady);
report=struct('status','passed-readonly-controlled-rocket-five-quantity-checkpoint-replay', ...
    'acceptedEndpointsReplayed',last.acceptedSteps,'cellCount',mesh.cellCount,'steps',steps, ...
    'maxBEDefectToOriginalTargetAndRoundoff',maxDefectRatio,'maxStepBalanceRelative',maxStepBalance, ...
    'globalBalanceRelative',balances(end),'originalFinalResidual',residuals(end), ...
    'steadyConverged',last.steadyConverged,'sourceSHAAndCaseRebuildPassed',true, ...
    'sourceEvidenceSHA256',fileSHA(fullfile(folder,'rocket_axisymmetric_evidence.mat')), ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false, ...
    'newProductionSolvedAreaM2',0);
if ~exist(auditFolder,'dir'),mkdir(auditFolder);end
save(fullfile(auditFolder,'checkpoint_replay_evidence.mat'),'report');
fid=fopen(fullfile(auditFolder,'checkpoint_replay_summary.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));clear guard;
fprintf('Rocket checkpoint independent replay: %d real endpoints, residual %.9g, five-quantity ledger %.9g.\n', ...
    last.acceptedSteps,residuals(end),balances(end));
end

function digest=fileSHA(path)
% 模块用途：按实际文件字节绑定源码与端点，不依赖文件名或修改时间。
fid=fopen(path,'rb');assert(fid>=0,'RocketCheckpointTest:File','无法读取%s。',path);
guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
raw=typecast(md.digest(),'uint8');digest=lower(reshape(dec2hex(raw,2)',1,[]));
end
