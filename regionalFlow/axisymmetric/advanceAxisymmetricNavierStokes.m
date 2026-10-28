function result=advanceAxisymmetricNavierStokes(Q,mesh,gas,boundary,options)
%ADVANCEAXISYMMETRICNAVIERSTOKES 轴对称贴体守恒量的真实全局步长SSP-RK2。
% 模块一：声学/对流、黏性/导热、圆柱几何扩散频率共同限制dt。
% 模块二：两阶段均保留压力/环向应力源；完整步失败才减半，拒绝裁值。
% 模块三：累计边界与源积分核对四量变化；径向动量不是无源Cartesian量。
% 模块四：目标时间/工作上限/真实连续残差稳态分别报告，未收敛不可假通过。
% 全部账为每弧度；壁温仍输入。当前尚无实际箭体/实验覆盖资格。
if nargin<5,options=struct;end
defaults=struct('maxSteps',2000,'targetTime',Inf,'convectiveCFL',.35,'diffusiveCFL',.35, ...
    'maximumTimeStep',Inf,'minimumTimeStep',1e-16,'maxRetries',10, ...
    'steadyResidualTolerance',0,'steadyRequiredSteps',8,'minimumSteadySteps',8, ...
    'conservationTolerance',1e-10,'residualScale',[]);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'AxisFV:Options','未知推进选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
for key={'maxSteps','steadyRequiredSteps','minimumSteadySteps'}
    validateattributes(options.(key{1}),{'numeric'},{'real','finite','scalar','integer','positive'});
end
validateattributes(options.maxRetries,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
for key={'convectiveCFL','diffusiveCFL','minimumTimeStep','conservationTolerance'}
    validateattributes(options.(key{1}),{'numeric'},{'real','finite','scalar','positive'});
end
for key={'targetTime','maximumTimeStep'}
    validateattributes(options.(key{1}),{'numeric'},{'real','scalar','positive','nonnan'});
end
validateattributes(options.steadyResidualTolerance,{'numeric'},{'real','finite','scalar','nonnegative'});
gas=axisymmetricFVGas(gas);P=axisymmetricFVState(Q,gas);V=mesh.volumePerRadian;
if isempty(options.residualScale)
    rho=max(P(:,1));speed=max(max(P(:,6)),max(vecnorm(P(:,2:3),2,2)));
    scale=[rho,rho*speed,rho*speed,max(abs(Q(:,4)))];
else,scale=options.residualScale;validateattributes(scale,{'numeric'},{'real','finite','positive','numel',4});scale=reshape(scale,1,4);end
initialQ=Q;initial=sum(Q.*V,1);totalScale=sum(V)*scale;
net=zeros(1,4);src=net;time=0;accepted=0;rejected=0;steadyCount=0;steady=false;
evaluations=0;[rate,info]=evaluate(Q);
residual=normRate(rate);history=struct('time',0,'dt',0,'residual',residual,'balance',0,'retries',0, ...
    'residualEvaluations',evaluations);
reason='maximum-steps-reached';reached=false;rejectLog=struct('step',{},'dt',{},'identifier',{},'message',{});
for step=1:options.maxSteps
    if time>=options.targetTime,reached=true;reason='target-time-reached';break;end
    dt=1/max(info.convectiveFrequency/options.convectiveCFL+info.diffusiveFrequency/options.diffusiveCFL);
    dt=min([dt,options.maximumTimeStep,options.targetTime-time]);
    assert(isfinite(dt)&&dt>0,'AxisFV:TimeStep','没有有限正步长。');
    for attempt=0:options.maxRetries
        assert(dt>=options.minimumTimeStep,'AxisFV:MinimumTimeStep','整步回退达到最小dt，没有裁状态。');
        try
            Q1=Q+dt*rate;axisymmetricFVState(Q1,gas);
            [r1,i1]=evaluate(Q1);
            trial=.5*Q+.5*(Q1+dt*r1);axisymmetricFVState(trial,gas);
            [rnew,inew]=evaluate(trial);
            bnew=net+.5*dt*(info.netBoundaryFlux+i1.netBoundaryFlux);
            snew=src+.5*dt*(info.sourceIntegral+i1.sourceIntegral);
            defect=sum(trial.*V,1)-initial+bnew-snew;balance=max(abs(defect)./totalScale);
            assert(isfinite(balance)&&balance<=options.conservationTolerance,'AxisFV:Conservation','四量改变与RK2边界/源积分不一致。');
            break;
        catch exception
            if ~ismember(exception.identifier,{'AxisFV:State','AxisFV:Temperature','AxisFV:HLLC','AxisFV:HLLE', ...
                    'AxisFV:AxisTraceState', ...
                    'AxisFV:TransportTemperature','AxisFV:TransportState','AxisFV:PressureSource','AxisFV:OutflowMach'})
                rethrow(exception);
            end
            rejected=rejected+1;rejectLog(end+1)=struct('step',step,'dt',dt,'identifier',exception.identifier,'message',exception.message); %#ok<AGROW>
            if attempt==options.maxRetries,rethrow(exception);end
            dt=dt/2;
        end
    end
    Q=trial;rate=rnew;info=inew;net=bnew;src=snew;time=time+dt;accepted=accepted+1;
    residual=normRate(rate);history.time(end+1)=time;history.dt(end+1)=dt;history.residual(end+1)=residual;
    history.balance(end+1)=balance;history.retries(end+1)=attempt;
    history.residualEvaluations(end+1)=evaluations;
    if residual<=options.steadyResidualTolerance,steadyCount=steadyCount+1;else,steadyCount=0;end
    if accepted>=options.minimumSteadySteps&&steadyCount>=options.steadyRequiredSteps
        steady=true;reason='steady-residual-converged';break;
    end
end
if time>=options.targetTime,reached=true;if ~steady,reason='target-time-reached';end,end
result=struct('Q',Q,'initialQ',initialQ,'time',time,'acceptedSteps',accepted,'rejectedAttempts',rejected, ...
    'steadyConverged',steady,'completedRequestedTime',reached,'stopReason',reason,'residualNorm',residual, ...
    'history',history,'rejectLog',rejectLog,'residualEvaluations',evaluations,'initialConservedIntegral',initial, ...
    'finalConservedIntegral',sum(Q.*V,1),'accumulatedBoundaryIntegral',net,'accumulatedSourceIntegral',src, ...
    'globalBalanceRelative',max(abs(sum(Q.*V,1)-initial+net-src)./totalScale), ...
    'finalInfo',info,'options',options,'residualScale',scale, ...
    'gas',gas,'geometryBinding',struct('X',mesh.X,'R',mesh.R,'cells',mesh.cells,'volumePerRadian',V), ...
    'wallTemperatureValueKind','prescribed','physicalSolution',false,'eligibleForSurfaceCoverage',false);
    function n=normRate(r),n=max(abs(r)./scale,[],'all');end
    function [r,i]=evaluate(state)
        % 只数完整空间算子，EOS检查不计；失败的真实调用仍属于工作量。
        evaluations=evaluations+1;
        [r,i]=axisymmetricNavierStokesResidual(state,mesh,gas,boundary);
    end
end
