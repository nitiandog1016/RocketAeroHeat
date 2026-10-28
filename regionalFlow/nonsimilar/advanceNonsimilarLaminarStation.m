function result=advanceNonsimilarLaminarStation(initial,history,gas,options,forcing)
%ADVANCENONSIMILARLAMINARSTATION 真实上游剖面驱动的一次隐式流向推进。
% 模块作用：自写有限差分Newton和回溯线搜索，联合解当前站u/T，
% 保留完整u_s、H_s与r_s连续耦合。壁u=0/T=Tw，外界u=Ue/T=Te；
% 外界v不另施加。失败返回NaN热流，不拿初猜或旧局部BVP冒充本步解。
% 本模块只证明离散方程一步求解；不代替驻点启动、域/网格/薄层与
% 实验验证。制造源必须显式传入；当前草稿始终不计生产面元覆盖。
% 正式迁移身份：受限条件研究模块；核心数值保持v61冻结草稿原式。
if nargin<4,options=struct;end
if nargin<5,forcing=struct;end
o=struct('maximumIterations',20,'relativeResidualTolerance',1e-9, ...
    'finiteDifferenceStep',2e-6,'minimumLineSearch',2^-16,'minimumRcond',1e-13);
keys=fieldnames(options);
for j=1:numel(keys)
    if ~isfield(o,keys{j}),error('RocketAeroHeat:NonSimilar:Options','未知选项%s。',keys{j});end
    o.(keys{j})=options.(keys{j});
end
validateattributes(o.maximumIterations,{'numeric'},{'scalar','integer','nonnegative','finite'});
for key={'relativeResidualTolerance','finiteDifferenceStep','minimumLineSearch','minimumRcond'}
    validateattributes(o.(key{1}),{'numeric'},{'scalar','finite','positive','<',1});
end
N=numel(initial.n);ids=2:N-1;M=numel(ids);
% 模块作用：给定边界在所有试探上精确施加，不能列作求得壁温。
c=initial;c.n=c.n(:);c.u=c.u(:);c.T=c.T(:);
c.u([1,end])=[0;c.Ue];c.T([1,end])=[c.Tw;c.Te];
us=max(c.Ue,1);ts=max([c.Te,c.Tw,1]);
z=[c.u(ids)/us;c.T(ids)/ts];historyNorm=[];lineHistory=[];reason='超过Newton工作上限';
result=struct('available',false,'converged',false,'physicalSolution',false, ...
    'eligibleForSurfaceCoverage',false,'profile',[], ...
    'wallHeatFluxIntoWall',NaN,'wallShearStress',NaN, ...
    'wallTemperatureInput',c.Tw,'wallTemperatureSolved',NaN, ...
    'statusCode',uint8(1),'status',reason,'iterations',0, ...
    'residualHistory',[],'lineSearchHistory',[],'residual',[], ...
    'method','物理n非相似薄层方程的流向隐式Newton一步', ...
    'scope','实验草稿；非驻点启动证明、非全区域/稳定性/实验精度验证');
for iteration=0:o.maximumIterations
    [f,detail,trial]=evaluate(z);normF=norm(f,inf);historyNorm(end+1)=normF; %#ok<AGROW>
    if normF<=o.relativeResidualTolerance
        result.available=true;result.converged=true;result.statusCode=uint8(0);
        result.status='离散非相似方程本步已收敛；物理区域有效性尚未验证';
        result.profile=trial;result.residual=detail;
        result.wallHeatFluxIntoWall=detail.wallHeatFluxIntoWall;
        result.wallShearStress=detail.wallShearStress;break;
    end
    if iteration==o.maximumIterations,break;end
    J=zeros(numel(f),numel(z));
    for col=1:numel(z)
        dz=o.finiteDifferenceStep*max(abs(z(col)),1);
        if col>M
            above=2500/ts-z(col);below=z(col)-80/ts;
            if dz<=above
                % 正方向试探仍在原物性域。
            elseif dz<=below
                dz=-dz;
            elseif above>=below
                dz=.5*above;
            else
                dz=-.5*below;
            end
            if abs(dz)<16*eps(max(abs(z(col)),1))
                error('RocketAeroHeat:NonSimilar:Jacobian','物性域内没有可解析的温度差分距离。');
            end
        end
        zp=z;zp(col)=zp(col)+dz;
        fp=evaluate(zp);J(:,col)=(fp-f)/dz;
    end
    rc=rcond(J);
    if ~isfinite(rc)||rc<o.minimumRcond,reason='Newton雅可比未解析或病态';break;end
    direction=-(J\f);alpha=1;accepted=false;
    while alpha>=o.minimumLineSearch
        candidate=z+alpha*direction;
        temps=candidate(M+1:end)*ts;vels=candidate(1:M)*us;
        if all(temps>=80&temps<=2500)&&all(vels>=0)
            fc=evaluate(candidate);
            if norm(fc,inf)<=(1-1e-4*alpha)*normF
                z=candidate;accepted=true;break;
            end
        end
        alpha=alpha/2;
    end
    lineHistory(end+1)=alpha; %#ok<AGROW>
    if ~accepted,reason='物性/无逆流约束内未取得下降的Newton步';break;end
end
result.iterations=numel(historyNorm)-1;result.residualHistory=historyNorm;
result.lineSearchHistory=lineHistory;
if ~result.available,result.status=reason;end

    function [f,d,t]=evaluate(zz)
        t=c;t.u(ids)=zz(1:M)*us;t.T(ids)=zz(M+1:end)*ts;
        d=nonsimilarLaminarResidual(t,history,gas,forcing);f=d.scaledResidual;
    end
end
