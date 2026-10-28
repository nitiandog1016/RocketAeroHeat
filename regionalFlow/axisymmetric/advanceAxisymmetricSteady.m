function result=advanceAxisymmetricSteady(Q,mesh,gas,boundary,options)
%ADVANCEAXISYMMETRICSTEADY 同一完整NS残差的全局后向Euler伪时间稳态推进。
% 模块一：四守恒量耦合Newton；稀疏有色差分包含声学、黏性、导热及几何源。
% 模块二：每个全局伪步真正解BE方程；正性/温域回溯，不裁剪状态或重标压力。
% 模块三：闭域质量由约束锚定初态，完整原残差仍独立验收；开放域不钉总质量。
% 模块四：逐步和累计核对BE边界/源积分及非线性缺陷，伪时间不是物理时间。
% 本模块不改变v62空间离散，不具备真实箭体/湍流/实验或生产覆盖资格。
if nargin<5,options=struct;end
defaults=struct('maxSteps',100,'initialPseudoCFL',5,'maximumPseudoCFL',1e6, ...
    'pseudoCFLGrowth',2,'maximumPseudoStep',Inf,'minimumPseudoStep',1e-16, ...
    'maxStepRetries',8,'maxNewtonIterations',12,'maxLineSearch',14, ...
    'nonlinearAbsoluteTolerance',1e-11,'nonlinearRelativeTolerance',1e-9, ...
    'conservationTolerance',5e-10,'steadyResidualTolerance',2e-6, ...
    'steadyRequiredSteps',3,'minimumSteadySteps',3,'residualScale',[], ...
    'differenceStep',2e-6,'differenceScheme','central','jacobianRefreshEvery',3, ...
    'massConstraint','auto');
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'AxisSteady:Options','未知稳态选项。');
for name=fieldnames(defaults)',if ~isfield(options,name{1}),options.(name{1})=defaults.(name{1});end,end
for name={'maxSteps','maxNewtonIterations','maxLineSearch','steadyRequiredSteps','minimumSteadySteps','jacobianRefreshEvery'}
    validateattributes(options.(name{1}),{'numeric'},{'real','finite','scalar','positive','integer'});
end
validateattributes(options.maxStepRetries,{'numeric'},{'real','finite','scalar','nonnegative','integer'});
for name={'initialPseudoCFL','maximumPseudoCFL','pseudoCFLGrowth','minimumPseudoStep', ...
        'nonlinearAbsoluteTolerance','nonlinearRelativeTolerance','conservationTolerance','differenceStep'}
    validateattributes(options.(name{1}),{'numeric'},{'real','finite','scalar','positive'});
end
validateattributes(options.maximumPseudoStep,{'numeric'},{'real','scalar','positive','nonnan'});
validateattributes(options.steadyResidualTolerance,{'numeric'},{'real','finite','scalar','nonnegative'});
assert(options.maximumPseudoCFL>=options.initialPseudoCFL&&options.pseudoCFLGrowth>=1&& ...
    ismember(options.differenceScheme,{'central','forward'})&&ismember(options.massConstraint,{'auto','off'}), ...
    'AxisSteady:Options','CFL范围、差分方式或质量约束合同无效。');
nq=size(Q,2);assert(ismember(nq,[4 5]),'AxisSteady:State','需要四量NS或五量SA完整状态。');
gas=axisymmetricFVGas(gas);P=axisymmetricFVState(Q(:,1:4),gas);
V=mesh.volumePerRadian;nc=mesh.cellCount;
assert(size(Q,1)==nc,'AxisSteady:Geometry','状态与控制体不匹配。');
if isempty(options.residualScale)
    rho=max(P(:,1));speed=max(max(P(:,6)),max(vecnorm(P(:,2:3),2,2)));
    scale=[rho,rho*speed,rho*speed,max(abs(Q(:,4)))];
    if nq==5
        % W=rho*N与mu同量纲；全零初场也给固定正参考，不在Jacobian内微分全局max。
        laminar=axisymmetricFVTransport(P(:,4),gas);scale(5)=max([max(abs(Q(:,5))),max(laminar.mu)]);
    end
else
    validateattributes(options.residualScale,{'numeric'},{'real','finite','positive','numel',nq});
    scale=reshape(options.residualScale,1,nq);
end
initialQ=Q;initial=sum(Q.*V,1);totalScale=sum(V)*scale;
U=Q./scale;massWeights=V/sum(V);massInitial=sum(massWeights.*U(:,1));
assert(isstruct(boundary)&&isscalar(boundary)&&all(isfield(boundary,{'iMin','iMax','jMin','jMax'})), ...
    'AxisSteady:Boundary','须明确给出四条逻辑边界。');
closed=strcmp(options.massConstraint,'auto');
for name={'iMin','iMax','jMin','jMax'}
    entry=boundary.(name{1});
    assert(isstruct(entry)&&isscalar(entry)&&isfield(entry,'type'),'AxisSteady:Boundary','每条边须给出单一type。');
    closed=closed&&ismember(boundary.(name{1}).type,{'wall','axis','periodic'});
end
constraint=sparse((1:nc)',ones(nc,1),massWeights,nq*nc,1);
[rate,info]=evaluate(Q);residual=normRate(rate);pseudoTime=0;accepted=0;rejected=0;
% 计数从首个完整残差开始，Jacobian返回的所有试探调用另行累计。
evaluations=1;linearSolves=0;factorizations=0;jacobians=0;pattern=[];
history=struct('pseudoTime',0,'dt',0,'residual',residual,'balance',0, ...
    'newtonIterations',0,'newtonDefect',0,'newtonTarget',0,'stepDefect',zeros(1,nq), ...
    'massRelativeDrift',0,'lineSearchTrials',0,'residualEvaluations',1,'linearDefect',0, ...
    'augmentedLinearDefect',0,'massConstraintDefect',0);
net=zeros(1,nq);src=net;steady=false;steadyCount=0;reason='maximum-steps-reached';
rejectLog=struct('step',{},'dt',{},'reason',{},'nonlinearDefect',{});
pseudoCFL=options.initialPseudoCFL;
for step=1:options.maxSteps
    % 原始稳定频率只用于选择伪步，不进入Jacobian，避免全局虚假耦合。
    explicit=1/max(info.convectiveFrequency+info.diffusiveFrequency);
    dt=min(options.maximumPseudoStep,pseudoCFL*explicit);
    assert(isfinite(dt)&&dt>=options.minimumPseudoStep,'AxisSteady:PseudoStep','伪步已无有效正值。');
    oldU=U;oldQ=Q;success=false;
    for attempt=0:options.maxStepRetries
        assert(dt>=options.minimumPseudoStep,'AxisSteady:MinimumPseudoStep','回退达到最小伪步，未裁状态。');
        trialU=oldU;trialRate=rate;trialInfo=info;
        g=-dt*(rate./scale);target=options.nonlinearAbsoluteTolerance+ ...
            options.nonlinearRelativeTolerance*max(abs(g),[],'all');
        cache=[];inner=0;lineTrials=0;linearDefect=0;augmentedDefect=0;massDefect=0;lastFailure='maximum-newton-iterations';
        for iteration=1:options.maxNewtonIterations
            g=trialU-oldU-dt*(trialRate./scale);gnorm=max(abs(g),[],'all');
            defect=sum((trialU-oldU).*scale.*V,1)+dt*(trialInfo.netBoundaryFlux-trialInfo.sourceIntegral);
            candidateBalance=max(abs(sum(trialU.*scale.*V,1)-initial+net+dt*trialInfo.netBoundaryFlux-src-dt*trialInfo.sourceIntegral)./totalScale);
            stepBalance=max(abs(defect)./totalScale);
            if gnorm<=target&&candidateBalance<=options.conservationTolerance&&stepBalance<=options.conservationTolerance
                success=true;break;
            end
            if isempty(cache)||mod(iteration-1,options.jacobianRefreshEvery)==0
                jacOptions=struct('differenceStep',options.differenceStep,'differenceScheme',options.differenceScheme);
                [J,jreport,pattern]=axisymmetricFVJacobian(trialU.*scale,mesh,gas,boundary,scale,jacOptions,pattern);
                evaluations=evaluations+jreport.residualEvaluations;jacobians=jacobians+1;
                if ~jreport.success,lastFailure='jacobian-perturbation-failed';break;end
                A=speye(nq*nc)-dt*J;
                if closed,system=[A,constraint;constraint',sparse(1,1)];else,system=A;end
                cache=decomposition(system,'lu');factorizations=factorizations+1;
            end
            rhs=-g(:);
            if closed,rhs=[rhs;massInitial-sum(massWeights.*trialU(:,1))];end
            increment=cache\rhs;linearSolves=linearSolves+1;
            du=reshape(increment(1:nq*nc),nc,nq);
            linearDefect=max(abs(A*du(:)+g(:)))/max(1,max(abs(g),[],'all'));
            augmentedDefect=max(abs(system*increment-rhs))/max(1,max(abs(rhs)));
            if closed,massDefect=abs(constraint'*du(:)-rhs(end));end
            if ~all(isfinite(du),'all'),lastFailure='nonfinite-linear-update';break;end
            found=false;alpha=1;
            for search=0:options.maxLineSearch
                lineTrials=lineTrials+1;
                try
                    candidate=trialU+alpha*du;
                    [newRate,newInfo]=evaluate(candidate.*scale);evaluations=evaluations+1;
                    newg=candidate-oldU-dt*(newRate./scale);newnorm=max(abs(newg),[],'all');
                    if newnorm<=(1-1e-4*alpha)*gnorm||newnorm<=target
                        found=true;break;
                    end
                catch exception
                    % 失败调用也属于实际工作量；输入/代码错误不能吞掉。
                    evaluations=evaluations+1;
                    if ~stateFailure(exception),rethrow(exception);end
                end
                alpha=alpha/2;
            end
            if ~found
                if mod(iteration-1,options.jacobianRefreshEvery)~=0
                    cache=[];continue;
                end
                lastFailure='line-search-failed';break;
            end
            trialU=candidate;trialRate=newRate;trialInfo=newInfo;inner=inner+1;
        end
        % 最后一次Newton更新之后仍需验收，不能少算一次而误判工作上限。
        g=trialU-oldU-dt*(trialRate./scale);gnorm=max(abs(g),[],'all');
        defect=sum((trialU-oldU).*scale.*V,1)+dt*(trialInfo.netBoundaryFlux-trialInfo.sourceIntegral);
        candidateBalance=max(abs(sum(trialU.*scale.*V,1)-initial+net+dt*trialInfo.netBoundaryFlux-src-dt*trialInfo.sourceIntegral)./totalScale);
        stepBalance=max(abs(defect)./totalScale);
        success=gnorm<=target&&candidateBalance<=options.conservationTolerance&&stepBalance<=options.conservationTolerance;
        if success,break;end
        rejected=rejected+1;rejectLog(end+1)=struct('step',step,'dt',dt,'reason',lastFailure,'nonlinearDefect',gnorm); %#ok<AGROW>
        dt=dt/4;
    end
    if ~success,reason='nonlinear-step-rejected';break;end
    U=trialU;Q=U.*scale;rate=trialRate;info=trialInfo;
    net=net+dt*info.netBoundaryFlux;src=src+dt*info.sourceIntegral;
    pseudoTime=pseudoTime+dt;accepted=accepted+1;residual=normRate(rate);
    balance=max(abs(sum(Q.*V,1)-initial+net-src)./totalScale);
    drift=(sum(Q(:,1).*V)-initial(1))/initial(1);
    assert(~closed||abs(drift)<=options.conservationTolerance,'AxisSteady:Mass','闭域初态质量未保持，没有重标rho/p。');
    history.pseudoTime(end+1)=pseudoTime;history.dt(end+1)=dt;history.residual(end+1)=residual;
    history.balance(end+1)=balance;history.newtonIterations(end+1)=inner;
    history.newtonDefect(end+1)=gnorm;history.newtonTarget(end+1)=target;
    history.stepDefect(end+1,:)=defect;history.massRelativeDrift(end+1)=drift;
    history.lineSearchTrials(end+1)=lineTrials;history.residualEvaluations(end+1)=evaluations;
    history.linearDefect(end+1)=linearDefect;
    history.augmentedLinearDefect(end+1)=augmentedDefect;history.massConstraintDefect(end+1)=massDefect;
    if residual<=options.steadyResidualTolerance,steadyCount=steadyCount+1;else,steadyCount=0;end
    if accepted>=options.minimumSteadySteps&&steadyCount>=options.steadyRequiredSteps
        steady=true;reason='steady-residual-converged';break;
    end
    % 回退后的实际dt决定下一次CFL，避免每步重复相同失败的大伪步。
    pseudoCFL=min(options.maximumPseudoCFL,options.pseudoCFLGrowth*dt/explicit);
    fprintf('AxisSteady step=%d dt=%.5g residual=%.9g Newton=%d calls=%d balance=%.3g\n', ...
        accepted,dt,residual,inner,evaluations,balance);
end
result=struct('Q',Q,'initialQ',initialQ,'pseudoTime',pseudoTime, ...
    'acceptedSteps',accepted,'rejectedAttempts',rejected,'steadyConverged',steady, ...
    'stopReason',reason,'residualNorm',residual,'residualScale',scale,'history',history, ...
    'rejectLog',rejectLog,'initialConservedIntegral',initial,'finalConservedIntegral',sum(Q.*V,1), ...
    'accumulatedBoundaryIntegral',net,'accumulatedSourceIntegral',src, ...
    'globalBalanceRelative',max(abs(sum(Q.*V,1)-initial+net-src)./totalScale), ...
    'finalInfo',info,'gas',gas,'options',options,'residualEvaluations',evaluations, ...
    'linearSolves',linearSolves,'factorizations',factorizations,'jacobianBuilds',jacobians, ...
    'closedMassConstraintActive',closed,'geometryBinding',struct('X',mesh.X,'R',mesh.R, ...
    'cells',mesh.cells,'volumePerRadian',V),'pseudoTimeIsPhysical',false,'physicalTimeAdvanced',0, ...
    'wallTemperatureValueKind','prescribed','physicalSolution',false,'eligibleForSurfaceCoverage',false);
    function [r,i]=evaluate(state)
        if nq==5,[r,i]=axisymmetricRANSResidual(state,mesh,gas,boundary);
        else,[r,i]=axisymmetricNavierStokesResidual(state,mesh,gas,boundary);end
    end
    function value=normRate(r),value=max(abs(r)./scale,[],'all');end
    function yes=stateFailure(e)
        yes=ismember(e.identifier,{'AxisFV:State','AxisFV:Temperature','AxisFV:HLLC','AxisFV:HLLE', ...
            'AxisSA:State','AxisSA:Arithmetic','AxisFV:SAModelState', ...
            'AxisFV:AxisTraceState','AxisSA:AxisTraceState', ...
            'AxisFV:TransportTemperature','AxisFV:TransportState','AxisFV:WallThermalState', ...
            'AxisFV:PressureSource','AxisFV:OutflowMach', ...
            'RocketAeroHeat:NSDraft:CharacteristicSonic','RocketAeroHeat:NSDraft:CharacteristicRegime', ...
            'RocketAeroHeat:NSDraft:CharacteristicAmbiguous','RocketAeroHeat:NSDraft:CharacteristicState', ...
            'RocketAeroHeat:NSDraft:CharacteristicSound'});
    end
end
