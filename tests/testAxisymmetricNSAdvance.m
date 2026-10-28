function report=testAxisymmetricNSAdvance(folder)
%TESTAXISYMMETRICNSADVANCE 环隙负X轴向Couette的真实推进、稳态和壁热校核。
% 模块一：从静止等温初态实际推进，不从解析终态伪造稳态历史。
% 模块二：比较独立积分的轴对称FV保守均值、完整速度/温度及墙热剪切。
% 模块三：检查内壁机械功、两墙导热、四量边界/几何源积分的闭合。
% 模块四：非零径向启动激活环向黏性源；区分工作上限、到时和稳态。
% 本例沿用人工常物性气体R=1、mu=.5，明确不是空气或TMR实验验证；
% 无轴环隙也不能认证r=0、任意弯曲网格、湍流或真实箭体全面积求解。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v62_validation','advance');end
if ~exist(folder,'dir'),mkdir(folder);end
bindings=sourceBindings(root);
gas=struct('gamma',1.4,'R',1,'mu',.5,'Pr',.72,'temperatureRange',[200 500]);
ri=.03;ro=.06;axialLength=.06;innerU=-3;outerU=0;Tw=300;initialRho=1;
periodic=struct('type','periodic');
boundary=struct('iMin',periodic,'iMax',periodic, ...
    'jMin',struct('type','wall','temperature',Tw,'velocity',[innerU 0]), ...
    'jMax',struct('type','wall','temperature',Tw,'velocity',[outerU 0]));
options=struct('maxSteps',100000,'convectiveCFL',.35,'diffusiveCFL',.35, ...
    'steadyResidualTolerance',2e-6,'steadyRequiredSteps',10,'minimumSteadySteps',10, ...
    'conservationTolerance',1e-10);
grids=[8 16];cases=struct([]);solutions=cell(size(grids));references=cell(size(grids));

%% 一、两档真实环隙从静止初态收敛，比较独立连续解及真正FV均值
for k=1:numel(grids)
    nr=grids(k);[X,R]=ndgrid([0 axialLength],linspace(ri,ro,nr+1));
    mesh=buildAxisymmetricFVGeometry(X,R,struct('periodicI',true));
    V=mesh.volumePerRadian;nc=mesh.cellCount;
    assert(abs(sum(V)-axialLength*(ro^2-ri^2)/2)<1e-15,'环隙每弧度体积不符。');
    initialP=[initialRho*ones(nc,1),zeros(nc,2),Tw*ones(nc,1)];
    initialQ=axisymmetricFVConservative(initialP,gas);
    mass=sum(initialQ(:,1).*V);
    exact=analyticReference(mesh,gas,ri,ro,axialLength,innerU,outerU,Tw,mass);
    timer=tic;
    result=advanceAxisymmetricNavierStokes(initialQ,mesh,gas,boundary,options);
    elapsed=toc(timer);
    save(fullfile(folder,sprintf('annular_n%02d_attempt.mat',nr)),'result','mesh','exact','elapsed');
    fprintf('Annular nr=%d: steady=%d steps=%d residual=%.12g seconds=%.3g time=%.9g.\n', ...
        nr,result.steadyConverged,result.acceptedSteps,result.residualNorm,elapsed,result.time);
    assert(result.steadyConverged&&strcmp(result.stopReason,'steady-residual-converged'), ...
        '环隙%d格只达到工作上限或目标时间，没有真实稳态。',nr);
    assert(result.acceptedSteps>=options.minimumSteadySteps&&~isequal(result.Q,initialQ));
    assert(result.residualNorm<=options.steadyResidualTolerance);
    assert(all(result.history.residual(end-options.steadyRequiredSteps+1:end)<=options.steadyResidualTolerance));
    assert(result.history.residual(1)>1e4*result.residualNorm,'没有可核对的残差下降。');
    assert(~result.physicalSolution&&~result.eligibleForSurfaceCoverage);
    assert(strcmp(result.wallTemperatureValueKind,'prescribed'));
    assert(result.globalBalanceRelative<=options.conservationTolerance);
    assert(all(isfinite(result.Q),'all')&&all(result.history.dt(2:end)>0));
    P=axisymmetricFVState(result.Q,gas);w=result.finalInfo.wall;
    inner=w.normal(:,2)<-.9;outer=w.normal(:,2)>.9;
    assert(any(inner)&&any(outer)&&all(inner|outer),'周期环隙只应有内外两组墙。');
    assert(all(w.temperatureInput==Tw));
    assert(all(w.velocityInput(inner,1)==innerU)&&all(w.velocityInput(outer,1)==outerU));
    assert(all(w.velocityInput(:,2)==0));

    % Q是真体积均值；EOS(Q均值)不等于质心温度或真体积平均温度。
    scale=[initialRho,initialRho*abs(innerU-outerU),initialRho*abs(innerU-outerU), ...
        initialRho*gas.R*Tw/(gas.gamma-1)];
    difference=(result.Q-exact.QMean)./scale;
    qMeanRelativeRMS=sqrt(sum(sum(difference.^2,2).*V)/sum(V));
    velocityRelativeRMS=sqrt(sum((P(:,2)-exact.pointVelocity).^2.*V)/sum(V))/abs(innerU-outerU);
    temperatureRiseRelativeRMS=sqrt(sum((P(:,4)-exact.pointTemperature).^2.*V)/sum(V))/exact.maximumTemperatureRise;
    maximumRadialSpeed=max(abs(P(:,3)));
    pressureNonuniformRelative=(max(P(:,5))-min(P(:,5)))/mean(P(:,5));
    pressureReferenceErrorPa=max(abs(P(:,5)-exact.uniformPressure));
    assert(all(P(:,2)<0)&&all(P(:,2)>innerU),'负X壁驱动方向或完整速度剖面不符。');
    assert(max(P(:,4))>Tw+.25*exact.maximumTemperatureRise,'未形成实际黏性发热温升。');
    assert(maximumRadialSpeed<2e-3&&pressureNonuniformRelative<5e-4, ...
        '轴向Couette终态仍有显著径向流或径向压力梯度。');
    assert(pressureReferenceErrorPa<.05,'闭盒质量确定的最终压力与连续参考差异过大。');
    assert(qMeanRelativeRMS<.06&&velocityRelativeRMS<.06&&temperatureRiseRelativeRMS<.30, ...
        '完整剖面或FV均值尚未达到首阶段数值校核门槛。');

    % 不同半径墙热密度不能设成同值；曲率与机械功必须分别检查。
    expectedHeat=zeros(size(w.heatFluxIntoWall));expectedShear=expectedHeat;
    expectedHeat(inner)=exact.innerHeatFlux;expectedHeat(outer)=exact.outerHeatFlux;
    expectedShear(inner)=exact.innerViscousTractionX;expectedShear(outer)=exact.outerViscousTractionX;
    wallHeatRelativeRMS=sqrt(mean((w.heatFluxIntoWall./expectedHeat-1).^2));
    wallShearRelativeRMS=sqrt(mean((w.viscousTraction(:,1)./expectedShear-1).^2));
    assert(all(w.heatFluxIntoWall>0),'同温两壁应接收黏性发热，不能用绝对值掩盖符号。');
    assert(all(w.viscousTraction(inner,1)<0)&&all(w.viscousTraction(outer,1)>0));
    assert(wallHeatRelativeRMS<.25&&wallShearRelativeRMS<.15, ...
        '环隙墙热或剪切不符合曲率解析参考。');
    totalWallHeat=sum(w.heatFluxIntoWall.*w.areaPerRadian);
    movingWallPower=sum(w.mechanicalPowerIntoFluidPerRadian);
    wallEnergyOut=sum(w.totalEnergyFluxOutPerRadian);
    powerClosure=abs(totalWallHeat-movingWallPower)/max(abs(movingWallPower),realmin);
    assert(movingWallPower>0&&powerClosure<1e-5,'内壁机械功未转为两墙总热。');
    assert(abs(wallEnergyOut-(totalWallHeat-movingWallPower))<1e-11, ...
        '墙能量通量必须逐值等于导热外流减机械功输入。');
    assert(abs(totalWallHeat/exact.totalWallHeatPower-1)<.15);
    assert(abs(movingWallPower/exact.totalWallHeatPower-1)<.15);
    checkIntegratedBalance(result,V);
    item=struct('radialCells',nr,'axialCells',1,'acceptedSteps',result.acceptedSteps, ...
        'seconds',elapsed,'time',result.time,'residualNorm',result.residualNorm, ...
        'initialResidualNorm',result.history.residual(1),'rejectedAttempts',result.rejectedAttempts, ...
        'QMeanRelativeRMS',qMeanRelativeRMS,'velocityRelativeRMS',velocityRelativeRMS, ...
        'temperatureRiseRelativeRMS',temperatureRiseRelativeRMS, ...
        'maximumRadialSpeedMps',maximumRadialSpeed,'pressureNonuniformRelative',pressureNonuniformRelative, ...
        'pressureReferenceErrorPa',pressureReferenceErrorPa,'expectedFinalPressurePa',exact.uniformPressure, ...
        'wallHeatRelativeRMS',wallHeatRelativeRMS,'wallShearRelativeRMS',wallShearRelativeRMS, ...
        'totalWallHeatPerRadianW',totalWallHeat,'movingWallPowerPerRadianW',movingWallPower, ...
        'expectedHeatPowerPerRadianW',exact.totalWallHeatPower,'steadyPowerRelativeClosure',powerClosure, ...
        'globalBalanceRelative',result.globalBalanceRelative);
    if isempty(cases),cases=item;else,cases(end+1)=item;end %#ok<AGROW>
    solutions{k}=struct('mesh',mesh,'result',result);references{k}=exact;
end
% 两档趋势检验不等于声称二阶；发现非单调须诊断，不能事后隐藏。
assert(cases(2).QMeanRelativeRMS<.9*cases(1).QMeanRelativeRMS, ...
    '真实FV均值误差未随8到16格减小。');
assert(cases(2).wallHeatRelativeRMS<.9*cases(1).wallHeatRelativeRMS, ...
    '真实墙热误差未随8到16格减小。');
assert(cases(2).temperatureRiseRelativeRMS<cases(1).temperatureRiseRelativeRMS, ...
    '黏性温升的连续解误差未随网格减小。');

%% 二、同一封闭环隙的径向启动：激活环向黏性源并实推进边界/源项账
mesh=solutions{1}.mesh;V=mesh.volumePerRadian;nc=mesh.cellCount;
radius=mesh.volumeCentroid(:,2);logRatio=log(ro/ri);amplitude=.2;
radialVelocity=amplitude*sin(pi*log(radius/ri)/logRatio);
initialP=[initialRho*ones(nc,1),zeros(nc,1),radialVelocity,Tw*ones(nc,1)];
radialInitialQ=axisymmetricFVConservative(initialP,gas);
[~,initialInfo]=axisymmetricNavierStokesResidual(radialInitialQ,mesh,gas,boundary);
% 独立连续积分：两壁vr=0，integral(tau_theta_theta)dA
% = (8/3)*mu*Lx*amplitude*ln(ro/ri)/pi。不能用p大项掩盖漏掉的hoop项。
exactHoopIntegral=(8/3)*gas.mu*axialLength*amplitude*logRatio/pi;
computedHoopIntegral=sum(initialInfo.sourceHoopStress.*mesh.planarArea);
hoopIntegralRelativeError=abs(computedHoopIntegral/exactHoopIntegral-1);
assert(hoopIntegralRelativeError<.20,'真实径向启动没有正确激活环向黏性源。');
pressureIntegral=initialRho*gas.R*Tw*sum(mesh.planarArea);
assert(abs((pressureIntegral-initialInfo.sourceIntegral(3))-computedHoopIntegral)<1e-11);
shortOptions=struct('targetTime',1e-5,'maxSteps',200,'convectiveCFL',.35, ...
    'diffusiveCFL',.35,'steadyResidualTolerance',0,'conservationTolerance',1e-10);
short=advanceAxisymmetricNavierStokes(radialInitialQ,mesh,gas,boundary,shortOptions);
assert(short.completedRequestedTime&&~short.steadyConverged&&strcmp(short.stopReason,'target-time-reached'));
assert(short.acceptedSteps>=2&&~isequal(short.Q,radialInitialQ));
assert(short.time==shortOptions.targetTime&&short.accumulatedSourceIntegral(3)>0);
assert(all(short.accumulatedSourceIntegral([1 2 4])==0));
checkIntegratedBalance(short,V);
assert(abs(short.finalConservedIntegral(1)-short.initialConservedIntegral(1)) ...
    <1e-12*short.initialConservedIntegral(1),'不可穿透周期环隙发生质量流失。');

%% 三、工作上限不能当稳态；完成物理时间也不能自动标记稳态
one=advanceAxisymmetricNavierStokes(radialInitialQ,mesh,gas,boundary, ...
    struct('maxSteps',1,'steadyResidualTolerance',0,'steadyRequiredSteps',8,'minimumSteadySteps',8));
assert(one.acceptedSteps==1&&~one.steadyConverged&&~one.completedRequestedTime);
assert(strcmp(one.stopReason,'maximum-steps-reached')&&one.time>0);
assert(~isequal(one.Q,radialInitialQ)&&~one.physicalSolution&&~one.eligibleForSurfaceCoverage);
checkIntegratedBalance(one,V);

%% 四、冻结完整测试证据；明确解析校核与实际空气/实验验证的区别
assert(isequal(bindings,sourceBindings(root)),'本测试运行期间模块源码改变，不能绑定同一通过结果。');
report=struct('status','passed-real-annular-NS-time-advance-and-analytic-steady-check', ...
    'cases',cases,'testGas',gas,'testGasIsAir',false,'axisIncluded',false, ...
    'radialStartupHoopIntegralExactPerRadianN',exactHoopIntegral, ...
    'radialStartupHoopIntegralComputedPerRadianN',computedHoopIntegral, ...
    'radialStartupHoopIntegralRelativeError',hoopIntegralRelativeError, ...
    'radialStartupAcceptedSteps',short.acceptedSteps,'radialStartupFinalTime',short.time, ...
    'radialStartupGlobalBalanceRelative',short.globalBalanceRelative, ...
    'maximumStepsIsNotConvergence',true,'targetTimeIsNotSteadyConvergence',true, ...
    'sourceBindings',bindings,'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'experimentalValidation',false,'newSolvedSurfaceAreaM2',0, ...
    'scope','人工常物性连续NS的无源环隙真实推进及解析稳态；非空气、非TMR、无轴点或真实箭体资格');
save(fullfile(folder,'axisymmetric_advance_evidence.mat'),'report','solutions','references','short','one', ...
    'radialInitialQ','initialInfo','boundary','gas','options');
fid=fopen(fullfile(folder,'axisymmetric_advance_summary.json'),'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s',jsonencode(report));
fprintf('Axisymmetric annular NS actual steady passed: Nr8/16=%d/%d steps, wall q errors %.6g/%.6g, radial transient %d steps; no physical coverage.\n', ...
    cases(1).acceptedSteps,cases(2).acceptedSteps,cases(1).wallHeatRelativeRMS,cases(2).wallHeatRelativeRMS,short.acceptedSteps);
end

function exact=analyticReference(mesh,gas,ri,ro,Lx,Ui,Uo,Tw,mass)
% 独立解析参考：严格按r权积分保守均值，不借用求解器梯度或通量。
L=log(ro/ri);du=Uo-Ui;cv=gas.R/(gas.gamma-1);cp=gas.gamma*cv;
k=gas.mu*cp/gas.Pr;xi=@(r)log(r/ri)/L;
u=@(r)Ui+du*xi(r);T=@(r)Tw+gas.mu*du^2/(2*k)*xi(r).*(1-xi(r));
p=mass*gas.R/(Lx*integral(@(r)r./T(r),ri,ro,'AbsTol',1e-14,'RelTol',1e-12));
rho=@(r)p./(gas.R*T(r));Q=zeros(mesh.cellCount,4);Tmean=zeros(mesh.cellCount,1);
for c=1:mesh.cellCount
    nodes=mesh.vertices(mesh.cells(c,:),2);lo=min(nodes);hi=max(nodes);radialVolume=(hi^2-lo^2)/2;
    average=@(fun)integral(@(r)r.*fun(r),lo,hi,'AbsTol',1e-13,'RelTol',1e-11)/radialVolume;
    Q(c,:)=[average(rho),average(@(r)rho(r).*u(r)),0, ...
        p/(gas.gamma-1)+.5*average(@(r)rho(r).*u(r).^2)];
    Tmean(c)=average(T);
end
exact=struct('QMean',Q,'trueVolumeMeanTemperature',Tmean,'uniformPressure',p, ...
    'pointVelocity',u(mesh.volumeCentroid(:,2)),'pointTemperature',T(mesh.volumeCentroid(:,2)), ...
    'maximumTemperatureRise',gas.mu*du^2/(8*k), ...
    'innerHeatFlux',gas.mu*du^2/(2*L*ri),'outerHeatFlux',gas.mu*du^2/(2*L*ro), ...
    'innerViscousTractionX',-gas.mu*du/(L*ri),'outerViscousTractionX',gas.mu*du/(L*ro), ...
    'totalWallHeatPower',Lx*gas.mu*du^2/L);
assert(abs(sum(Q(:,1).*mesh.volumePerRadian)-mass)<1e-13*mass);
end

function checkIntegratedBalance(result,V)
% 径向动量含几何源；不能照搬笛卡尔闭盒四量全无源的余额。
change=result.finalConservedIntegral-result.initialConservedIntegral;
defect=change+result.accumulatedBoundaryIntegral-result.accumulatedSourceIntegral;
scale=sum(V)*reshape(result.residualScale,1,4);
relative=max(abs(defect)./scale);
assert(relative<1e-10&&abs(relative-result.globalBalanceRelative)<1e-14);
assert(all(isfinite(result.accumulatedSourceIntegral))&&all(isfinite(result.accumulatedBoundaryIntegral)));
end

function bindings=sourceBindings(root)
% 完整模块与本测试逐文件绑定，历史运行不冒称预锁源码。
names={'advanceAxisymmetricNavierStokes.m','axisymmetricNavierStokesResidual.m', ...
    'axisymmetricFVGas.m','axisymmetricFVState.m','axisymmetricFVConservative.m', ...
    'buildAxisymmetricFVGeometry.m'};
paths=cellfun(@(s)['regionalFlow/axisymmetric/' s],names,'UniformOutput',false);
paths{end+1}='tests/testAxisymmetricNSAdvance.m';paths=sort(paths);bindings=struct([]);
for j=1:numel(paths)
    path=paths{j};fid=fopen(fullfile(root,path),'rb');assert(fid>=0);
    cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup;
    digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
    item=struct('path',path,'sha256',lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[])));
    if isempty(bindings),bindings=item;else,bindings(end+1)=item;end %#ok<AGROW>
end
end
