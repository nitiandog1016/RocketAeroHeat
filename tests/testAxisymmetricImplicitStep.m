function report=testAxisymmetricImplicitStep(folder)
%TESTAXISYMMETRICIMPLICITSTEP 原NS残差的一/两步全局后向Euler、质量与墙热合同。
% 模块一：小型周期环隙使用非均匀rho/ux/ur/T，激活真实共享通量和环向源。
% 模块二：独立复算完整BE方程和四量缺陷，不用Newton收敛冒充稳态收敛。
% 模块三：两步结果与逐步重启互核，检查闭域质量及接受步的边界/源积分。
% 模块四：给定Tw与流体温度分开；伪时间不是物理时间，物理解面积仍为0。
% 仅运行一两步及输入拒绝；不启动环隙长稳态、TMR推进或默认生产UI。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v63_validation','implicit_step');end
bindings=sourceBindings(root);
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
[X,R]=ndgrid([0 .02 .05 .08],[.03 .035 .045 .06]);
mesh=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));
V=mesh.volumePerRadian;nc=mesh.cellCount;xc=mesh.volumeCentroid(:,1);rc=mesh.volumeCentroid(:,2);
phase=2*pi*xc/.08;xi=log(rc/.03)/log(2);
primitive=[1+.03*cos(phase)+.02*sin(pi*xi), ...
    -.4+.08*sin(phase).*(.5+xi), ...
    .12*sin(pi*xi).*(1+.15*cos(phase)), ...
    300+.2*cos(phase)+.1*sin(pi*xi)];
Q=axisymmetricFVConservative(primitive,gas);initialQ=Q;
periodic=struct('type','periodic');
boundary=struct('iMin',periodic,'iMax',periodic, ...
    'jMin',struct('type','wall','temperature',300,'velocity',[-3 0]), ...
    'jMax',struct('type','wall','temperature',302,'velocity',[0 0]));
scale=[1 20 20 800];dt=1e-7;
options=struct('maxSteps',1,'initialPseudoCFL',5,'maximumPseudoCFL',20, ...
    'maximumPseudoStep',dt,'steadyResidualTolerance',0, ...
    'steadyRequiredSteps',3,'minimumSteadySteps',3,'residualScale',scale, ...
    'nonlinearAbsoluteTolerance',1e-11,'nonlinearRelativeTolerance',1e-9, ...
    'conservationTolerance',5e-10,'differenceScheme','central', ...
    'differenceStep',2e-6,'jacobianRefreshEvery',1);
[initialRate,initialInfo]=axisymmetricNavierStokesResidual(Q,mesh,gas,boundary);
assert(nc==9&&max(abs(initialInfo.sourceHoopStress))>1e-5&& ...
    initialInfo.sourceIntegral(3)>0&&all(initialInfo.sourceIntegral([1 2 4])==0), ...
    'AxisSteadyTest:Excitation','初态必须激活非零径向速度/环向应力及柱坐标源。');
assert(max(abs(initialRate)./scale,[],'all')>0);
assert(max(abs(sum(initialRate.*V,1)+initialInfo.netBoundaryFlux-initialInfo.sourceIntegral))<1e-11);

%% 一、真实一步完整BE与固定非线性门槛；步数上限不得被当作稳态
one=advanceAxisymmetricSteady(Q,mesh,gas,boundary,options);
assert(one.acceptedSteps==1&&~one.steadyConverged&& ...
    strcmp(one.stopReason,'maximum-steps-reached')&&one.rejectedAttempts==0, ...
    'AxisSteadyTest:OneStep','应完成恰好一个小BE步，不能放大为稳态或隐瞒回退。');
assert(~isequaln(one.Q,Q)&&isequaln(one.initialQ,Q)&&isequal(one.residualScale,scale));
assert(abs(one.pseudoTime-dt)<64*eps(dt)&&abs(one.history.dt(2)-dt)<64*eps(dt));
oneAudit=auditStep(Q,one.Q,dt,one,2);
auditResult(one,1);
assert(one.history.newtonIterations(2)>=1&&one.linearSolves>=1&&one.residualEvaluations>=3, ...
    'AxisSteadyTest:Work','非均匀一步须有实际Newton/线性求解及完整残差调用。');

%% 二、两步与逐次单步重启一致；同时复算每步与累计四量/闭域质量
twoOptions=options;twoOptions.maxSteps=2;
two=advanceAxisymmetricSteady(Q,mesh,gas,boundary,twoOptions);
restarted=advanceAxisymmetricSteady(one.Q,mesh,gas,boundary,options);
assert(two.acceptedSteps==2&&restarted.acceptedSteps==1&& ...
    ~two.steadyConverged&&strcmp(two.stopReason,'maximum-steps-reached')&& ...
    two.rejectedAttempts==0&&restarted.rejectedAttempts==0);
assert(abs(two.pseudoTime-2*dt)<64*eps(2*dt)&&all(abs(two.history.dt(2:end)-dt)<64*eps(dt)));
restartDifference=max(abs(two.Q-restarted.Q)./scale,[],'all');
assert(restartDifference<5e-10,'AxisSteadyTest:Restart','固定全局dt与scale下两步和逐步重启不一致。');
secondAudit=auditStep(one.Q,two.Q,dt,two,3);
auditStep(one.Q,restarted.Q,dt,restarted,2);
auditResult(two,2);auditResult(restarted,1);
[~,oneInfo]=axisymmetricNavierStokesResidual(one.Q,mesh,gas,boundary);
[~,twoInfo]=axisymmetricNavierStokesResidual(two.Q,mesh,gas,boundary);
expectedBoundary=dt*(oneInfo.netBoundaryFlux+twoInfo.netBoundaryFlux);
expectedSource=dt*(oneInfo.sourceIntegral+twoInfo.sourceIntegral);
assert(max(abs(two.accumulatedBoundaryIntegral-expectedBoundary)./(sum(V)*scale))<1e-12&& ...
    max(abs(two.accumulatedSourceIntegral-expectedSource)./(sum(V)*scale))<1e-12, ...
    'AxisSteadyTest:EndpointIntegral','BE积分须使用每个接受步终态，不能沿用RK2两阶段账。');
mass=initialQ(:,1)'*V;
massDrift=[one.Q(:,1)'*V-mass,two.Q(:,1)'*V-mass]/mass;
assert(max(abs(massDrift))<=options.conservationTolerance, ...
    'AxisSteadyTest:Mass','闭域初态总质量没有保持。');
assert(all(two.history.newtonDefect(2:end)<=two.history.newtonTarget(2:end))&& ...
    all(two.history.residual(2:end)>0)&&~two.steadyConverged, ...
    'AxisSteadyTest:ConvergenceIdentity','BE缺陷小不能代替原始NS稳态残差。');

%% 三、墙热、机械功和给定Tw始终由当前态原残差派生；温度身份不得混淆
wall=two.finalInfo.wall;state=axisymmetricFVState(two.Q,gas);
inner=wall.normal(:,2)<0;outer=~inner;
assert(all(wall.temperatureInput(inner)==300)&&all(wall.temperatureInput(outer)==302)&& ...
    all(wall.velocityInput(inner,1)==-3)&&all(wall.velocityInput(outer,1)==0)&& ...
    all(wall.velocityInput(:,2)==0));
assert(all(isfinite(wall.pressure)&wall.pressure>0)&&all(isfinite(wall.heatFluxIntoWall))&& ...
    any(abs(wall.heatFluxIntoWall)>1e-8)&&isequal(wall.pressure,state(wall.ownerCell,5)));
assert(any(abs(state(wall.ownerCell,4)-wall.temperatureInput)>1e-3), ...
    'AxisSteadyTest:TemperatureIdentity','流体状态温度不能被覆写为给定墙温。');
wallEnergyIdentity=max(abs(wall.totalEnergyFluxOutPerRadian- ...
    (wall.heatFluxIntoWall.*wall.areaPerRadian-wall.mechanicalPowerIntoFluidPerRadian)));
assert(wallEnergyIdentity<1e-11,'AxisSteadyTest:WallEnergy','墙总能量须等于导热外流减机械功输入。');
assert(isequaln(two.finalInfo,twoInfo),'AxisSteadyTest:StaleWall','终态墙热或通量不是当前Q的原始空间残差。');
for key={'wallTemperatureSolved','recoveryTemperature','edgeTemperature'}
    if isfield(two,key{1})
        assert(all(isnan(two.(key{1})),'all'),'AxisSteadyTest:UnknownTemperature','未定义温度不能由流体或边界输入补齐。');
    end
end
assert(isequaln(Q,initialQ),'测试不得改写输入守恒量。');

%% 四、原状态、温域、壁条件与伪时间输入非法时，必须明确拒绝而非裁值
refusals=struct('label',{},'identifier',{});
bad=Q;bad(1,1)=-1;reject('negative-density',@()advanceAxisymmetricSteady(bad,mesh,gas,boundary,options));
bad=complex(Q,ones(size(Q)));reject('complex-state',@()advanceAxisymmetricSteady(bad,mesh,gas,boundary,options));
bad=Q;u=bad(1,2:3)/bad(1,1);bad(1,4)=bad(1,1)*(gas.R/(gas.gamma-1)*199+sum(u.^2)/2);
reject('initial-temperature-out-of-range',@()advanceAxisymmetricSteady(bad,mesh,gas,boundary,options));
badGas=gas;badGas.mu=-1;reject('negative-viscosity',@()advanceAxisymmetricSteady(Q,mesh,badGas,boundary,options));
badBoundary=rmfield(boundary,'iMax');reject('missing-boundary',@()advanceAxisymmetricSteady(Q,mesh,gas,badBoundary,options));
badBoundary=boundary;badBoundary.jMin.temperature=NaN;
reject('nonfinite-prescribed-Tw',@()advanceAxisymmetricSteady(Q,mesh,gas,badBoundary,options));
badBoundary=boundary;badBoundary.jMin.velocity=[-3 .1];
reject('penetrating-wall',@()advanceAxisymmetricSteady(Q,mesh,gas,badBoundary,options));
badBoundary=boundary;badBoundary.jMax=struct('type','axis');
reject('finite-radius-axis',@()advanceAxisymmetricSteady(Q,mesh,gas,badBoundary,options));
badOptions=options;badOptions.unknownOption=true;
reject('unknown-option',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));
badOptions=options;badOptions.maximumPseudoStep=0;
reject('nonpositive-pseudo-step',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));
badOptions=options;badOptions.differenceStep=0;
reject('nonpositive-difference-step',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));
badOptions=options;badOptions.differenceScheme='complex';
reject('unsupported-difference-scheme',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));
badOptions=options;badOptions.residualScale=[1 0 1 1];
reject('nonpositive-residual-scale',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));
badOptions=options;badOptions.maxNewtonIterations=0;
reject('nonpositive-Newton-work',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));
badOptions=options;badOptions.steadyResidualTolerance=NaN;
reject('nonfinite-steady-tolerance',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));
badOptions=options;badOptions.massConstraint='pressure';
reject('unsupported-mass-constraint',@()advanceAxisymmetricSteady(Q,mesh,gas,boundary,badOptions));

%% 五、冻结实际测试证据；未进行长稳态、TMR或生产区域验收
assert(isequal(bindings,sourceBindings(root)),'AxisSteadyTest:Source','测试期间源码改变，不能绑定同一通过结果。');
report=struct('status','passed-one-two-global-BE-steps-and-closed-mass-audit', ...
    'cellCount',nc,'globalPseudoStepS',dt,'oneStep',oneAudit,'secondStep',secondAudit, ...
    'restartMaximumScaledDifference',restartDifference,'closedMassRelativeDrift',massDrift, ...
    'wallEnergyIdentityMaximumPerRadianW',wallEnergyIdentity, ...
    'residualEvaluations',[one.residualEvaluations two.residualEvaluations restarted.residualEvaluations], ...
    'linearSolves',[one.linearSolves two.linearSolves restarted.linearSolves], ...
    'refusals',refusals,'refusalCount',numel(refusals),'sourceBindings',bindings, ...
    'backwardEulerIsNotSteadyConvergence',true,'pseudoTimeIsPhysical',false, ...
    'physicalTimeAdvanced',0,'wallTemperatureIsPrescribed',true,'longSteadyRunPerformed',false, ...
    'tmrRunPerformed',false,'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'experimentalValidation',false,'newSolvedSurfaceAreaM2',0);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'implicit_step_evidence.mat'),'report','mesh','gas','boundary', ...
    'initialQ','primitive','initialRate','initialInfo','options','one','two','restarted');
fid=fopen(fullfile(folder,'implicit_step_summary.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report,'ConvertInfAndNaN',true));
fprintf('Implicit BE step contract passed: %d cells, step=%.9g, mass drift %.3g, %d refusals; no steady/TMR/physical coverage.\n', ...
    nc,dt,max(abs(massDrift)),numel(refusals));

    function item=auditStep(oldQ,newQ,currentDt,current,index)
        [oldRate,~]=axisymmetricNavierStokesResidual(oldQ,mesh,gas,boundary);
        [newRate,newInfo]=axisymmetricNavierStokesResidual(newQ,mesh,gas,boundary);
        G=newQ-oldQ-currentDt*newRate;
        defectNorm=max(abs(G)./scale,[],'all');
        target=options.nonlinearAbsoluteTolerance+ ...
            options.nonlinearRelativeTolerance*currentDt*max(abs(oldRate)./scale,[],'all');
        roundoff=512*eps(max(1,max(abs(newQ)./scale,[],'all')));
        assert(defectNorm<=target+roundoff&& ...
            abs(current.history.newtonTarget(index)-target)<=roundoff&& ...
            abs(current.history.newtonDefect(index)-defectNorm)<=roundoff, ...
            'AxisSteadyTest:BEEquation','完整BE方程或固定Newton目标不匹配原始残差。');
        integral=sum((newQ-oldQ).*V,1)+currentDt*(newInfo.netBoundaryFlux-newInfo.sourceIntegral);
        weightedG=sum(G.*V,1);totalScale=sum(V)*scale;
        assert(max(abs(integral-weightedG)./totalScale)<1e-12&& ...
            max(abs(current.history.stepDefect(index,:)-integral)./totalScale)<1e-12, ...
            'AxisSteadyTest:StepDefect','四量缺陷须等于原BE缺陷的真实体积积分。');
        assert(abs(current.history.residual(index)-max(abs(newRate)./scale,[],'all'))<1e-9, ...
            'AxisSteadyTest:OriginalResidual','稳态历史须为原NS残差，不能用BE缺陷或预条件残差替代。');
        item=struct('scaledBEdefect',defectNorm,'nonlinearTarget',target, ...
            'integratedStepDefect',integral,'volumeIntegralOfBEdefect',weightedG, ...
            'originalNSResidual',max(abs(newRate)./scale,[],'all'));
    end
    function auditResult(current,steps)
        assert(current.acceptedSteps==steps&&numel(current.history.dt)==steps+1&& ...
            numel(current.history.pseudoTime)==steps+1&&isequal(size(current.history.stepDefect),[steps+1 4]));
        assert(current.history.dt(1)==0&&current.history.pseudoTime(1)==0&& ...
            all(current.history.stepDefect(1,:)==0)&&current.history.newtonDefect(1)==0&& ...
            current.history.newtonTarget(1)==0&&current.history.balance(1)==0&&current.history.newtonIterations(1)==0);
        assert(~current.pseudoTimeIsPhysical&&current.physicalTimeAdvanced==0&&~isfield(current,'time')&& ...
            ~current.physicalSolution&&~current.eligibleForSurfaceCoverage&& ...
            strcmp(current.wallTemperatureValueKind,'prescribed'));
        assert(all(isfinite(current.Q),'all')&&all(current.history.dt(2:end)>0)&& ...
            abs(current.pseudoTime-sum(current.history.dt))<64*eps(max(current.pseudoTime,realmin)));
        assert(isequal(current.geometryBinding.X,mesh.X)&&isequal(current.geometryBinding.R,mesh.R)&& ...
            isequal(current.geometryBinding.cells,mesh.cells)&&isequal(current.geometryBinding.volumePerRadian,V));
        assert(isequal(current.gas,axisymmetricFVGas(gas)));
        change=sum(current.Q.*V,1)-sum(current.initialQ.*V,1);
        globalDefect=change+current.accumulatedBoundaryIntegral-current.accumulatedSourceIntegral;
        relative=max(abs(globalDefect)./(sum(V)*scale));
        assert(relative<=options.conservationTolerance&&abs(relative-current.globalBalanceRelative)<1e-13&& ...
            max(abs(sum(current.history.stepDefect,1)-globalDefect)./(sum(V)*scale))<1e-12);
        currentMass=sum(current.Q(:,1).*V);oldMass=sum(current.initialQ(:,1).*V);
        assert(abs(currentMass-oldMass)/oldMass<=options.conservationTolerance);
    end
    function reject(label,action)
        try
            action();
        catch exception
            assert(startsWith(exception.identifier,{'AxisSteady:','AxisFV:','MATLAB:expected', ...
                'MATLAB:invalidType','MATLAB:validators:','RocketAeroHeat:NSDraft:'}), ...
                'AxisSteadyTest:UnexpectedRefusal','拒绝%s时出现未分类代码错误：%s。',label,exception.identifier);
            refusals(end+1)=struct('label',label,'identifier',exception.identifier); %#ok<AGROW>
            return;
        end
        error('AxisSteadyTest:MissingRefusal','非法输入%s未被拒绝。',label);
    end
end

function bindings=sourceBindings(root)
% 绑定当前区域模块和本测试的完整源码；本测试不把历史运行当当前源码证据。
files=dir(fullfile(root,'regionalFlow','axisymmetric','*.m'));
paths=cellfun(@(name)['regionalFlow/axisymmetric/' name],{files.name},'UniformOutput',false);
paths{end+1}='tests/testAxisymmetricImplicitStep.m';paths=sort(paths);bindings=struct([]);
for j=1:numel(paths)
    fid=fopen(fullfile(root,paths{j}),'rb');assert(fid>=0,'AxisSteadyTest:Source','缺少来源文件。');
    cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
    digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
    item=struct('path',paths{j},'sha256',lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[])));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
end
end
