function result=advanceCompressibleNavierStokes3D(Q,grid,solid,gas,boundary,options)
%ADVANCECOMPRESSIBLENAVIERSTOKES3D 三维黏性/导热保守量的全局步长SSP-RK2推进。
% 模块作用：调用空间残差内核，按声学/对流和黏性/热扩散频率选择
% 共享dt。整步两阶段失败时原状态不变、dt减半；有限次后明确报错，
% 不裁密度、内能、温度，也不以HLLC非法状态假造收敛。
% Q为Nx*Ny*Nz*5，SI，solid内部保持NaN。定比热气体/物性及边界
% 限制继承所选空间内核。options.residualFunction可显式提供算子。
% targetTime达到只说明完成该段物理时间；真实稳态还须归一空间残差
% 连续steadyRequiredSteps步不超过steadyResidualTolerance [1/s]。
% maxSteps仅是工作上限，达到不能标converged。两阶段净边界通量积分
% 和保守量总改变单独核对。最终墙Tw标输入，q/p由最终流体态派生；
% 无网格/边界收敛或实验验证，不自动计入表面物理解面积。
if nargin<6,options=struct;end
options=defaults(options);h=double(grid.spacing(:)');
validateattributes(h,{'double'},{'finite','positive','numel',3});
shape=[size(Q,1),size(Q,2),size(Q,3)];n=prod(shape);volume=prod(h);
assert(size(Q,4)==5 && ndims(Q)<=4,'RocketAeroHeat:NSAdvance:Shape','Q须为Nx*Ny*Nz*5。');
assert(islogical(solid) && isequal([size(solid,1),size(solid,2),size(solid,3)],shape), ...
    'RocketAeroHeat:NSAdvance:Mask','solid须与流体网格一致。');
fluid=~solid(:);assert(any(fluid),'RocketAeroHeat:NSAdvance:NoFluid','没有流体控制体。');
flat=reshape(double(Q),n,5);flat(~fluid,:)=NaN;Q=reshape(flat,[shape,5]);
[primitive0,cp,cv,k]=state(Q);
if isempty(options.residualScale)
    rhoScale=max(primitive0(:,1));speedScale=max(sqrt(gas.gamma*gas.R*primitive0(:,5)));
    velocityScale=max(vecnorm(primitive0(:,2:4),2,2));momentumScale=rhoScale*max(speedScale,velocityScale);
    scale=[rhoScale,momentumScale,momentumScale,momentumScale,max(abs(flat(fluid,5)))];
else
    scale=double(options.residualScale(:)');
    validateattributes(scale,{'double'},{'real','finite','positive','numel',5});
end
totalScale=scale*(nnz(fluid)*volume);initialTotal=total(Q);
initialState=Q;time=options.initialTime;initialTime=time;accumulatedBoundary=zeros(1,5);
operator=options.residualFunction;
[R,currentInfo]=evaluate(Q);currentNorm=normResidual(R);initialNorm=currentNorm;
history=struct('time',time,'acceptedDt',0,'proposedDt',0,'retries',0, ...
    'convectiveFrequency',frequency(primitive0),'diffusiveFrequency',diffusive(primitive0), ...
    'residualNorm',currentNorm,'conservationError',0);
rejectLog=struct('step',{},'attempt',{},'dt',{},'identifier',{},'message',{});
steadyCount=0;steady=false;targetReached=false;stopReason='maximum-steps-reached';
accepted=0;rejections=0;
for step=1:options.maxSteps
    if reachedTarget(time),targetReached=true;stopReason='target-time-reached';break;end
    [p,~,~,~]=state(Q);convective=frequency(p);diffusion=diffusive(p);
    % 合并频率比单独取min更保守，保持所有流体单元使用相同物理dt。
    proposed=1/(convective/options.convectiveCFL+diffusion/options.diffusiveCFL);
    proposed=min(proposed,options.maximumTimeStep);
    if isfinite(options.targetTime),proposed=min(proposed,options.targetTime-time);end
    assert(isfinite(proposed) && proposed>0,'RocketAeroHeat:NSAdvance:TimeStep','没有有限正时间步。');
    trial=proposed;done=false;
    for attempt=0:options.maxRetries
        if trial<options.minimumTimeStep
            error('RocketAeroHeat:NSAdvance:MinimumTimeStep', ...
                '第%d步dt=%.17g低于minimumTimeStep=%.17g；没有裁物理状态。',step,trial,options.minimumTimeStep);
        end
        try
            stage1=combine(Q,1,Q,0,R,trial);state(stage1);
            [R1,info1]=evaluate(stage1);
            candidate=combine(Q,.5,stage1,.5,R1,.5*trial);state(candidate);
            % 接受前还评估终态；非法HLLC终态不能进入下一步或墙面结果。
            [Rnew,infoNew]=evaluate(candidate);newNorm=normResidual(Rnew);
            boundaryStep=.5*trial*(currentInfo.netBoundaryFlux+info1.netBoundaryFlux);
            candidateBoundary=accumulatedBoundary+boundaryStep;
            errorVector=total(candidate)-initialTotal+candidateBoundary;
            balance=max(abs(errorVector)./totalScale);
            assert(balance<=options.conservationTolerance, ...
                'RocketAeroHeat:NSAdvance:Conservation','累计五守恒量改变与RK2边界积分不一致。');
            done=true;break;
        catch exception
            % 模块作用：只有候选流体/波系/特征状态失效才整步减dt。
            % 程序、算子接口或守恒合同错误原样抛出，不能重复试13次掩盖它。
            recoverable={'RocketAeroHeat:NSAdvance:State','RocketAeroHeat:NSAdvance:InternalEnergy', ...
                'RocketAeroHeat:NSAdvance:Temperature','RocketAeroHeat:NSDraft:State', ...
                'RocketAeroHeat:NSDraft:InternalEnergy','RocketAeroHeat:NSDraft:Temperature', ...
                'RocketAeroHeat:NSDraft:HLLC','RocketAeroHeat:NSDraft:CharacteristicState', ...
                'RocketAeroHeat:NSDraft:CharacteristicSound','RocketAeroHeat:NSDraft:CharacteristicRegime', ...
                'RocketAeroHeat:NSDraft:CharacteristicSonic','RocketAeroHeat:NSDraft:CharacteristicAmbiguous', ...
                'RocketAeroHeat:NSDraft:OutflowMach'};
            if ~any(startsWith(exception.identifier,recoverable)),rethrow(exception);end
            rejections=rejections+1;
            rejectLog(end+1)=struct('step',step,'attempt',attempt+1,'dt',trial, ...
                'identifier',exception.identifier,'message',exception.message); %#ok<AGROW>
            if attempt==options.maxRetries
                failure=MException('RocketAeroHeat:NSAdvance:RetryLimit', ...
                    '第%d步在%d次候选失败后仍无可接受NS状态；最后dt=%.17g，原步未提交。', ...
                    step,attempt+1,trial);failure=addCause(failure,exception);throw(failure);
            end
            trial=trial/2;
        end
    end
    assert(done,'RocketAeroHeat:NSAdvance:RetryLimit','整步没有被接受。');
    % 只有完整RK2和终态通过才同时提交状态、时间及边界积分。
    Q=candidate;R=Rnew;currentInfo=infoNew;currentNorm=newNorm;
    accumulatedBoundary=candidateBoundary;time=time+trial;accepted=accepted+1;
    history.time(end+1)=time;history.acceptedDt(end+1)=trial;history.proposedDt(end+1)=proposed;
    history.retries(end+1)=attempt;history.convectiveFrequency(end+1)=convective;
    history.diffusiveFrequency(end+1)=diffusion;history.residualNorm(end+1)=currentNorm;
    history.conservationError(end+1)=balance;
    if options.steadyEnabled && currentNorm<=options.steadyResidualTolerance
        steadyCount=steadyCount+1;
    else
        steadyCount=0;
    end
    if options.steadyEnabled && accepted>=options.minimumSteadySteps && steadyCount>=options.steadyRequiredSteps
        steady=true;stopReason='steady-residual-converged';break;
    end
    if reachedTarget(time),targetReached=true;stopReason='target-time-reached';break;end
end
finalTotal=total(Q);globalBalance=finalTotal-initialTotal+accumulatedBoundary;
result=struct('Q',Q,'initialState',initialState,'time',time,'initialTime',initialTime, ...
    'stopReason',stopReason,'completedRequestedTime',targetReached,'steadyConverged',steady, ...
    'converged',steady,'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'acceptedSteps',accepted,'rejectedAttempts',rejections,'rejections',rejectLog,'history',history, ...
    'residual',R,'residualNorm',currentNorm,'initialResidualNorm',initialNorm,'residualScale',scale, ...
    'residualNormUnits','1/s；保守量按固定初始尺度归一后的最大空间残差', ...
    'initialConservedIntegral',initialTotal,'finalConservedIntegral',finalTotal, ...
    'conservedIntegralChange',finalTotal-initialTotal,'accumulatedBoundaryIntegral',accumulatedBoundary, ...
    'globalBalance',globalBalance,'globalBalanceRelative',max(abs(globalBalance)./totalScale), ...
    'conservedIntegralUnits',{{'kg','kg*m/s','kg*m/s','kg*m/s','J'}}, ...
    'finalSpatialInfo',currentInfo,'wallFaces',currentInfo.wallFaces, ...
    'wallTemperatureValueKind','prescribed','wallPressureValueKind','derived-from-final-fluid-state', ...
    'wallHeatFluxValueKind','derived-from-final-fluid-state','wallFieldsConverged',steady, ...
    'options',options,'cp',cp,'cv',cv,'k',k,'positivityClippingCount',0, ...
    'scope','真实全局SSP-RK2推进；没有网格/边界收敛或实验验证，未收敛终态是暂态，不计覆盖');
result.options.residualFunction=func2str(operator);

    function [p,cpv,cvv,kv]=state(q)
        f=reshape(q,n,5);f=f(fluid,:);
        assert(all(isfinite(f),'all') && isreal(f) && all(f(:,1)>0), ...
            'RocketAeroHeat:NSAdvance:State','流体rho/保守量无效；不裁密度。');
        rho=f(:,1);u=f(:,2:4)./rho;internal=f(:,5)-.5*rho.*sum(u.^2,2);
        assert(all(isfinite(internal) & internal>0),'RocketAeroHeat:NSAdvance:InternalEnergy', ...
            'RK候选内能非正；需要整步减dt而非裁能量。');
        cvv=gas.R/(gas.gamma-1);cpv=gas.gamma*cvv;T=internal./(rho*cvv);
        assert(all(isfinite(T) & T>=gas.temperatureRange(1) & T<=gas.temperatureRange(2)), ...
            'RocketAeroHeat:NSAdvance:Temperature','RK候选T超出显式物性域；不裁温。');
        p=[rho,u,T];kv=gas.mu*cpv/gas.Pr;
        if isfield(gas,'k'),assert(abs(gas.k/kv-1)<1e-12,'RocketAeroHeat:NSAdvance:Pr','k/Pr合同不一致。');kv=gas.k;end
    end
    function value=frequency(p)
        a=sqrt(gas.gamma*gas.R*p(:,5));value=max(sum((abs(p(:,2:4))+a)./h,2));
    end
    function value=diffusive(p)
        diffusion=max((4/3)*gas.mu./p(:,1),k./(p(:,1)*cv));
        value=2*max(diffusion)*sum(1./h.^2);
    end
    function [rate,info]=evaluate(q)
        [rate,info]=operator(q,grid,reshape(solid,shape),gas,boundary);
        assert(isequal(size(rate),size(q)),'RocketAeroHeat:NSAdvance:Operator','残差算子输出尺寸不符。');
        values=reshape(rate,n,5);assert(all(isfinite(values(fluid,:)),'all'), ...
            'RocketAeroHeat:NSAdvance:Operator','流体空间残差非有限。');
        assert(all(isnan(values(~fluid,:)),'all'),'RocketAeroHeat:NSAdvance:Solid','固体残差须NaN。');
        assert(isfield(info,'netBoundaryFlux') && numel(info.netBoundaryFlux)==5 && ...
            all(isfinite(info.netBoundaryFlux)),'RocketAeroHeat:NSAdvance:Operator','算子须给有限五量净边界通量。');
        info.netBoundaryFlux=reshape(info.netBoundaryFlux,1,5);
        assert(isfield(info,'wallFaces'),'RocketAeroHeat:NSAdvance:Operator','算子须保存最终墙面派生量。');
    end
    function q=combine(a,wa,b,wb,rate,dt)
        aa=reshape(a,n,5);bb=reshape(b,n,5);rr=reshape(rate,n,5);qq=nan(n,5);
        qq(fluid,:)=wa*aa(fluid,:)+wb*bb(fluid,:)+dt*rr(fluid,:);q=reshape(qq,[shape,5]);
    end
    function v=total(q)
        f=reshape(q,n,5);v=sum(f(fluid,:),1)*volume;
    end
    function v=normResidual(r)
        f=reshape(r,n,5);v=max(abs(f(fluid,:))./scale,[],'all');
    end
    function value=reachedTarget(t)
        value=isfinite(options.targetTime) && t>=options.targetTime;
    end
end
function o=defaults(o)
% 模块作用：拒绝不合理CFL和缺失停止目标；步数上限不代表收敛。
d=struct('residualFunction',@compressibleNavierStokesResidual3D,'initialTime',0, ...
    'targetTime',Inf,'steadyResidualTolerance',NaN,'steadyRequiredSteps',5,'minimumSteadySteps',5, ...
    'maxSteps',5000,'convectiveCFL',.45,'diffusiveCFL',.45,'maxRetries',12, ...
    'minimumTimeStep',1e-14,'maximumTimeStep',Inf,'conservationTolerance',1e-9,'residualScale',[]);
for key=fieldnames(d)'
    if ~isfield(o,key{1}),o.(key{1})=d.(key{1});end
end
assert(isa(o.residualFunction,'function_handle'),'RocketAeroHeat:NSAdvance:Operator','residualFunction须为函数句柄。');
for key={'convectiveCFL','diffusiveCFL'}
    v=o.(key{1});assert(isscalar(v) && isfinite(v) && v>0 && v<=1, ...
        'RocketAeroHeat:NSAdvance:CFL','CFL须在(0,1]；不自动裁CFL。');
end
for key={'maxSteps','maxRetries','steadyRequiredSteps','minimumSteadySteps'}
    v=o.(key{1});assert(isscalar(v) && isfinite(v) && v==floor(v) && v>=0, ...
        'RocketAeroHeat:NSAdvance:Options','计数控制须非负有限整数。');
end
assert(o.maxSteps>0 && o.steadyRequiredSteps>0,'RocketAeroHeat:NSAdvance:Options','步数和稳态持续步数须正。');
validateattributes(o.minimumTimeStep,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(o.conservationTolerance,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(o.initialTime,{'numeric'},{'real','finite','scalar','nonnegative'});
assert(isscalar(o.maximumTimeStep) && o.maximumTimeStep>0,'RocketAeroHeat:NSAdvance:Options','最大dt须正。');
o.steadyEnabled=isfinite(o.steadyResidualTolerance);
if o.steadyEnabled,assert(o.steadyResidualTolerance>0,'RocketAeroHeat:NSAdvance:Options','稳态残差门限须正。');end
assert(isscalar(o.targetTime) && o.targetTime>=o.initialTime && ...
    (isfinite(o.targetTime) || o.steadyEnabled),'RocketAeroHeat:NSAdvance:Target', ...
    '须给有限targetTime或真正稳态残差门限；不能只给步数称收敛。');
end
