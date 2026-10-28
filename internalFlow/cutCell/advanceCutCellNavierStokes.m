function result=advanceCutCellNavierStokes(Q,cut,gas,boundary,options)
%ADVANCECUTCELLNAVIERSTOKES 真实变体积曲面控制体的全局SSP-RK2实验推进器。
% 模块用途：全体流体组共享物理dt，按原体积积分五量及两阶段边界账。
% 实际声学面通量尺度和LS/完整位移修正扩散模板的绝对行和用于选步；
% 它们是刚度估计，不是非正交离散、特征边界或非线性NS稳定性证明。
% 每个阶段检查状态和新步长估计；失败撤销整步，不裁密度/内能/温度。
% 达到targetTime或工作步数不等于稳态。steadyResidualTolerance需明确。
% 本模块只推进固定工况流体场，不是轨迹/累计壁热/壁体瞬态响应扩展。
if nargin<5,options=struct;end
o=defaults(options);v=cut.cells.volume(:);n=numel(v);fluid=v>0;
assert(isequal(size(Q),[n 5])&&isreal(Q),'CutAdvance:Shape','Q须为真实控制体数×5。');
assert(all(isfinite(v)&v>=0)&&any(fluid),'CutAdvance:Volume','须有正真实流体体积。');
assert(all(isnan(Q(~fluid,:)),'all'),'CutAdvance:Solid','固体保守量必须NaN，不自动覆盖非法输入。');
assert(all(isfield(gas,{'gamma','R','mu','Pr','temperatureRange'})),'CutAdvance:Gas','显式常物性不完整。');
validateattributes([gas.gamma gas.R gas.mu gas.Pr],{'numeric'},{'real','finite','positive'});
assert(gas.gamma>1 && numel(gas.temperatureRange)==2 && all(isfinite(gas.temperatureRange)) && ...
    gas.temperatureRange(1)>0 && gas.temperatureRange(2)>gas.temperatureRange(1),'CutAdvance:Gas','物性域无效。');
cv=gas.R/(gas.gamma-1);cp=gas.gamma*cv;k=gas.mu*cp/gas.Pr;
if isfield(gas,'k'),assert(abs(gas.k/k-1)<1e-12,'CutAdvance:Gas','k/Pr合同不一致。');end
P0=state(Q);initialState=Q;
if isempty(o.residualScale)
    density=max(P0(fluid,1));speed=max([sqrt(gas.gamma*gas.R*P0(fluid,5));vecnorm(P0(fluid,2:4),2,2)]);
    scale=[density,repmat(density*speed,1,3),max(Q(fluid,5))];
else,scale=o.residualScale(:)';validateattributes(scale,{'numeric'},{'real','finite','positive','numel',5});end
integralScale=scale*sum(v);initialTotal=total(Q);accumulated=zeros(1,5);
time=o.initialTime;[R,info]=evaluate(Q);initialSpatialInfo=info;initialNorm=residualNorm(R);currentNorm=initialNorm;
stencil=[];stencilBuilds=0;[dtBound,initialStep]=stepBound(Q);
history=struct('time',time,'acceptedDt',0,'proposedDt',0,'retries',0, ...
    'acousticFrequency',initialStep.acousticFrequency,'diffusiveFrequency',initialStep.diffusiveFrequency, ...
    'residualNorm',currentNorm,'conservationError',0);
rejects=struct('step',{},'attempt',{},'dt',{},'identifier',{},'message',{});
trace=cell(0,1);accepted=0;rejected=0;steady=false;steadyCount=0;targetReached=false;
reason='maximum-steps-reached';
for step=1:o.maxSteps
    if reachedTarget(time),targetReached=true;reason='target-time-reached';break;end
    [dtBound,currentStep]=stepBound(Q);proposed=min(dtBound,o.targetTime-time);trial=proposed;
    before=Q;beforeR=R;beforeInfo=info;done=false;
    for attempt=0:o.maxRetries
        assert(isfinite(trial)&&trial>=o.minimumTimeStep,'CutAdvance:MinimumTimeStep', ...
            '共享dt低于显式下限；没有裁状态或删除小单元。');
        try
            stage=combine(Q,1,Q,0,R,trial);state(stage);[R1,info1]=evaluate(stage);
            [stageBound,~]=stepBound(stage);
            assert(trial<=stageBound*(1+1e-12),'CutAdvance:StageTimeStep','第一阶段状态改变使共享步长超过新的刚度估计。');
            candidate=combine(Q,.5,stage,.5,R1,.5*trial);state(candidate);[Rnew,newInfo]=evaluate(candidate);
            boundaryStep=.5*trial*(info.netBoundaryFlux+info1.netBoundaryFlux);
            nextAccumulated=accumulated+boundaryStep;
            balance=max(abs(total(candidate)-initialTotal+nextAccumulated)./integralScale);
            assert(balance<=o.conservationTolerance,'CutAdvance:Conservation','累计五量改变不符两阶段边界积分。');
            done=true;break;
        catch err
            % 模块用途：只有候选物理态/波系错误可减步；几何、合同和程序错误原样停止。
            recoverable={'CutAdvance:State','CutAdvance:InternalEnergy','CutAdvance:Temperature','CutAdvance:StageTimeStep', ...
                'CutNS:State','CutNS:InternalEnergy','CutNS:Temperature','CutNS:HLLC','CutNS:HLLCStarPressure', ...
                'CutNS:BoundaryMach','RocketAeroHeat:NSDraft:CharacteristicState', ...
                'RocketAeroHeat:NSDraft:CharacteristicSound','RocketAeroHeat:NSDraft:CharacteristicRegime', ...
                'RocketAeroHeat:NSDraft:CharacteristicSonic','RocketAeroHeat:NSDraft:CharacteristicAmbiguous'};
            if ~any(strcmp(err.identifier,recoverable)),rethrow(err);end
            rejected=rejected+1;rejects(end+1)=struct('step',step,'attempt',attempt+1,'dt',trial, ...
                'identifier',err.identifier,'message',err.message); %#ok<AGROW>
            if attempt==o.maxRetries
                failure=MException('CutAdvance:RetryLimit','第%d步候选失败；原状态/时间/边界积分未提交。',step);
                throw(addCause(failure,err));
            end
            trial=trial/2;
        end
    end
    assert(done,'CutAdvance:RetryLimit','完整步未接受。');
    Q=candidate;R=Rnew;info=newInfo;accumulated=nextAccumulated;time=time+trial;accepted=accepted+1;
    currentNorm=residualNorm(R);
    if o.logEvery>0 && mod(accepted,o.logEvery)==0
        fprintf('CUT_ADVANCE_ACCEPT step=%d dt=%.9g t=%.9g residual=%.9g retries=%d balance=%.3g\n', ...
            accepted,trial,time,currentNorm,attempt,balance);
    end
    history.time(end+1)=time;history.acceptedDt(end+1)=trial;history.proposedDt(end+1)=proposed;
    history.retries(end+1)=attempt;history.acousticFrequency(end+1)=currentStep.acousticFrequency;
    history.diffusiveFrequency(end+1)=currentStep.diffusiveFrequency;
    history.residualNorm(end+1)=currentNorm;history.conservationError(end+1)=balance;
    if o.saveStageTrace
        trace{end+1}=struct('initialQ',before,'initialRate',beforeR,'stageQ',stage,'stageRate',R1, ...
            'finalQ',Q,'finalRate',R,'dt',trial,'initialBoundaryFlux',beforeInfo.netBoundaryFlux, ...
            'stageBoundaryFlux',info1.netBoundaryFlux,'integratedBoundary',boundaryStep); %#ok<AGROW>
    end
    if o.steadyEnabled && currentNorm<=o.steadyResidualTolerance,steadyCount=steadyCount+1;else,steadyCount=0;end
    if o.steadyEnabled && accepted>=o.minimumSteadySteps && steadyCount>=o.steadyRequiredSteps
        steady=true;reason='steady-residual-converged';break;
    end
    if reachedTarget(time),targetReached=true;reason='target-time-reached';break;end
end
finalTotal=total(Q);balanceVector=finalTotal-initialTotal+accumulated;
result=struct('Q',Q,'initialQ',initialState,'time',time,'initialTime',o.initialTime, ...
    'stopReason',reason,'completedRequestedTime',targetReached,'steadyConverged',steady, ...
    'targetTimeDifference',o.targetTime-time,'targetTimeRoundoffTolerance',timeRoundoff(time), ...
    'acceptedSteps',accepted,'rejectedAttempts',rejected,'rejections',rejects,'history',history, ...
    'residual',R,'residualNorm',currentNorm,'initialResidualNorm',initialNorm,'residualScale',scale, ...
    'initialConservedIntegral',initialTotal,'finalConservedIntegral',finalTotal, ...
    'accumulatedBoundaryIntegral',accumulated,'globalBalance',balanceVector, ...
    'globalBalanceRelative',max(abs(balanceVector)./integralScale),'initialSpatialInfo',initialSpatialInfo,'finalSpatialInfo',info, ...
    'stageTrace',{trace},'stencilBuilds',stencilBuilds,'initialStepEstimate',initialStep, ...
    'wallTemperatureValueKind','prescribed','wallPressureValueKind','derived-final-fluid-state', ...
    'wallHeatFluxValueKind','derived-final-fluid-state','wallFieldsSteadyConverged',steady, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'timeIntegratorImplemented',true,'fullNavierStokesStabilityProven',false, ...
    'positivityClippingCount',0,'options',o, ...
    'scope','变体积共享dt SSP-RK2实验推进；没有近壁/网格/外域/物性/湍流或实验精度验证');
result.options.residualFunction=func2str(o.residualFunction);

    function P=state(q)
        z=q(fluid,:);assert(isreal(z)&&all(isfinite(z),'all')&&all(z(:,1)>0),'CutAdvance:State','流体rho/保守量非法。');
        rho=z(:,1);u=z(:,2:4)./rho;ei=z(:,5)-.5*rho.*sum(u.^2,2);
        assert(all(isfinite(ei)&ei>0),'CutAdvance:InternalEnergy','候选内能非正，不裁能量。');
        T=ei./(rho*cv);assert(all(isfinite(T)&T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)), ...
            'CutAdvance:Temperature','候选温度超出明确物性域，不裁温。');
        P=nan(n,5);P(fluid,:)=[rho,u,T];
    end
    function [r,i]=evaluate(q)
        [r,i]=o.residualFunction(q,cut,gas,boundary);
        assert(isequal(size(r),[n 5])&&all(isfinite(r(fluid,:)),'all')&&all(isnan(r(~fluid,:)),'all'), ...
            'CutAdvance:Operator','空间残差尺寸、有限值或固体合同失效。');
        assert(isfield(i,'netBoundaryFlux')&&numel(i.netBoundaryFlux)==5&&all(isfinite(i.netBoundaryFlux)), ...
            'CutAdvance:Operator','算子没有有限五量边界通量。');
        i.netBoundaryFlux=reshape(i.netBoundaryFlux,1,5);
    end
    function [dt,s]=stepBound(q)
        p=state(q);b=cutCellBoundaryStencilState(q,cut,gas,boundary);
        if isempty(stencil)||~isequal(stencil.incoming,b.incoming)||~isequal(stencil.wall,b.wall)
            if o.logEvery>0,fprintf('CUT_ADVANCE_STENCIL_BUILD %d\n',stencilBuilds+1);end
            stencil=buildCutCellDiffusionStencil(cut,b);stencilBuilds=stencilBuilds+1;
        end
        frequency=zeros(n,1);f=cut.openFaces;w=cut.walls;a=sqrt(gas.gamma*gas.R*p(:,5));
        for kk=find(f.area>0)'
            m=f.cellMinus(kk);pp=f.cellPlus(kk);normal=f.normal(kk,:);
            ids=[m pp];ids=ids(ids>0);wave=max(abs(p(ids,2:4)*normal')+a(ids));
            if m==0||pp==0
                given=b.primitive(kk,:);wave=max(wave,abs(dot(given(2:4),normal))+sqrt(gas.gamma*gas.R*given(5)));
            end
            for id=ids,frequency(id)=frequency(id)+f.area(kk)*wave/v(id);end
        end
        for kk=find(w.ownerCell>0)'
            id=w.ownerCell(kk);wave=abs(dot(p(id,2:4),w.normal(kk,:)))+a(id);
            frequency(id)=frequency(id)+w.area(kk)*wave/v(id);
        end
        thermal=k./(p(fluid,1)*cv).*stencil.heatAbsoluteRowSum(fluid);
        stress=gas.mu./repmat(p(:,1),3,1).*stencil.stressAbsoluteRowSum;
        diffusion=max([thermal;stress(repmat(fluid,3,1))]);acoustic=max(frequency(fluid));
        dt=min(o.maximumTimeStep,1/(acoustic/o.acousticCFL+diffusion/o.diffusiveCFL));
        assert(isfinite(dt)&&dt>0,'CutAdvance:TimeStep','没有可表示的正共享dt。');
        s=struct('acousticFrequency',acoustic,'diffusiveFrequency',diffusion, ...
            'method','真实共享面声速 + 冻结边界LS修正扩散矩阵绝对行和', ...
            'boundaryStateSensitivityIncluded',false,'fullNavierStokesStabilityProven',false);
    end
    function q=combine(a,wa,b,wb,r,dt)
        q=nan(n,5);q(fluid,:)=wa*a(fluid,:)+wb*b(fluid,:)+dt*r(fluid,:);
    end
    function z=total(q),z=sum(q(fluid,:).*v(fluid),1);end
    function z=residualNorm(r),z=max(abs(r(fluid,:))./scale,[],'all');end
    function z=timeRoundoff(t)
        if isfinite(o.targetTime),z=8*eps(max(abs([t,o.targetTime])));else,z=0;end
    end
    function z=reachedTarget(t)
        % 只接受可量化的浮点终时差；实际积分时间保留，不补走或篡改边界账。
        z=isfinite(o.targetTime) && o.targetTime-t<=timeRoundoff(t);
    end
end
function o=defaults(o)
d=struct('residualFunction',@cutCellNavierStokesResidual,'initialTime',0,'targetTime',Inf, ...
    'maxSteps',1000,'maxRetries',12,'minimumTimeStep',1e-14,'maximumTimeStep',Inf, ...
    'acousticCFL',.35,'diffusiveCFL',.35,'conservationTolerance',1e-9,'residualScale',[], ...
    'steadyResidualTolerance',NaN,'steadyRequiredSteps',5,'minimumSteadySteps',5,'saveStageTrace',false,'logEvery',0);
for key=fieldnames(d)',if ~isfield(o,key{1}),o.(key{1})=d.(key{1});end,end
assert(isa(o.residualFunction,'function_handle'),'CutAdvance:Options','须明确函数句柄。');
for key={'acousticCFL','diffusiveCFL'},validateattributes(o.(key{1}),{'numeric'},{'real','finite','scalar','positive','<=',1});end
for key={'maxSteps','maxRetries','steadyRequiredSteps','minimumSteadySteps','logEvery'}
    validateattributes(o.(key{1}),{'numeric'},{'real','finite','scalar','integer','nonnegative'});
end
assert(o.maxSteps>0&&o.steadyRequiredSteps>0,'CutAdvance:Options','工作步数及稳态持续步数须正。');
validateattributes(o.initialTime,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(o.minimumTimeStep,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(o.conservationTolerance,{'numeric'},{'real','finite','scalar','positive'});
assert(isscalar(o.maximumTimeStep)&&isreal(o.maximumTimeStep)&&o.maximumTimeStep>0,'CutAdvance:Options','最大dt非法。');
assert(isscalar(o.steadyResidualTolerance)&&isreal(o.steadyResidualTolerance),'CutAdvance:Options','稳态门限须实标量。');
o.steadyEnabled=isfinite(o.steadyResidualTolerance);
assert(~o.steadyEnabled||o.steadyResidualTolerance>0,'CutAdvance:Options','稳态门限须正。');
assert(isscalar(o.targetTime)&&isreal(o.targetTime)&&o.targetTime>=o.initialTime&& ...
    (isfinite(o.targetTime)||o.steadyEnabled),'CutAdvance:Target','须明确目标时间或稳态残差门限，步数上限不是收敛。');
assert(isscalar(o.saveStageTrace)&&islogical(o.saveStageTrace),'CutAdvance:Options','轨迹保存须逻辑标量。');
end
