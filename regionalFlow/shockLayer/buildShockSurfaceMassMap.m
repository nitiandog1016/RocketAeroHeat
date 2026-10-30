function map=buildShockSurfaceMassMap(sampler,q,apex,input)
%BUILDSHOCKSURFACEMASSMAP 一般凸共同激波图的投影质量连续式及同流线入口映射。
% 模块一：dq/dtau=grad(f)，沿同一梯度线解div(a*grad(f))=1；a的量纲为m。
% 模块二：先向驻点回溯几何/面积，再从驻点方向稳定推进a与log(psi)。
% 模块三：只保存数值网格与Hermite导数，不保存闭包或把旧矩阵当作新入口。
% 范围：均匀上游、显式驻点、沿回溯线Hessian正定的单值图；不求驻点整体流场。
if nargin<4,input=struct;end
assert(isa(sampler,'function_handle'),'ShockLayer:MassSampler','需要同一共同曲面的采样器。');
validateattributes(q,{'double'},{'real','finite','size',[2 1]});
validateattributes(apex,{'double'},{'real','finite','size',[2 1]});
backTime=25;if isfield(input,'backTime'),backTime=input.backTime;end
relativeTolerance=2e-11;if isfield(input,'relativeTolerance'),relativeTolerance=input.relativeTolerance;end
gridRateStep=.015;if isfield(input,'gridRateStep'),gridRateStep=input.gridRateStep;end
validateattributes(backTime,{'double'},{'real','finite','scalar','>=',12,'<=',45});
validateattributes(relativeTolerance,{'double'},{'real','finite','scalar','positive','<=',1e-8});
validateattributes(gridRateStep,{'double'},{'real','finite','scalar','positive','<=',.05});
start=sampler(apex);H0=start.graphHessian;g0=start.graphGradient;
assert(all(size(H0)==[2 2])&&norm(H0-H0','fro')<1e-11*norm(H0,'fro')&&min(eig(H0))>0, ...
    'ShockLayer:MassConvexity','显式驻点Hessian必须对称正定。');
referenceRate=min(eig(H0));lengthScale=norm(q-apex);
assert(norm(g0)<1e-12,'ShockLayer:MassApex','给定横向驻点的梯度不为零。');
assert(lengthScale>0,'ShockLayer:StagnationLimit','有限法向列不能替代驻点整体极限。');
evaluations=0;largestRate=trace(H0)/referenceRate;
options=odeset('RelTol',relativeTolerance,'AbsTol',[1e-13 1e-13 1e-12 1e-13], ...
    'MaxStep',.12);
geometrySolution=ode45(@geometryRHS,[0 backTime],[(q-apex)/lengthScale;0;0],options);
tail=geometrySolution.y(:,end);tailQ=apex+lengthScale*tail(1:2);last=sampler(tailQ);
tailHessianIndicator=norm(last.graphHessian-H0,'fro')/norm(H0,'fro');
assert(tailHessianIndicator<1e-8,'ShockLayer:MassTail','回溯尾部尚未分辨驻点线性极限。');
% I是从当前点向驻点累积的投影面积因子；余项使用明确的驻点线性极限。
integralA=(tail(4)+exp(tail(3))*referenceRate/trace(H0))/referenceRate;
massOptions=odeset('RelTol',relativeTolerance,'AbsTol',[1e-13 1e-12],'MaxStep',.12);
massSolution=ode45(@massRHS,[backTime 0],[referenceRate/trace(H0);0],massOptions);
current=massSolution.y(:,end);coefficient=current(1)/referenceRate;
twoPass=abs(coefficient-integralA)/max(coefficient,realmin);
assert(coefficient>0&&isfinite(coefficient)&&twoPass<1e-8, ...
    'ShockLayer:MassContinuity','稳定质量推进与独立面积积分没有一致。');

%% 一、对数质量比严格单调；Hermite节点按实际最大曲率速率加密
segments=ceil(backTime*largestRate/gridRateStep);
assert(segments<=50000,'ShockLayer:MassBudget','当前曲率跨度超出显式质量映射存储预算。');
grid=unique([linspace(0,backTime,segments+1) geometrySolution.x massSolution.x]);
state=deval(geometrySolution,grid);a=deval(massSolution,grid);coordinate=apex+lengthScale*state(1:2,:);
coordinateDerivative=zeros(size(coordinate));
for j=1:numel(grid)
    sample=sampler(coordinate(:,j));coordinateDerivative(:,j)=-sample.graphGradient/referenceRate;
end
logRatio=a(2,:)-current(2);logDerivative=-1./a(1,:);
assert(all(a(1,:)>0)&&all(diff(logRatio)<0)&&abs(logRatio(1))<1e-10, ...
    'ShockLayer:MassContinuity','质量坐标必须从当前点向驻点严格递减。');
local=sampler(q);massPerDensitySpeed=coefficient*norm(local.graphGradient);
map=struct('method','共同凸激波图div(a*grad(f))=1；同一梯度线投影流量与log(psi)输运', ...
    'coordinate',q,'apexCoordinate',apex,'apexHessian',H0,'coefficientM',coefficient, ...
    'projectedMassVectorPerDensitySpeedM',coefficient*local.graphGradient, ...
    'massPerDensitySpeedM',massPerDensitySpeed,'referenceRatePerM',referenceRate, ...
    'traceParameter',grid,'sourceCoordinateGridM',coordinate,'coordinateDerivative',coordinateDerivative, ...
    'logMassRatio',logRatio,'logMassDerivative',logDerivative,'coefficientGridM',a(1,:)/referenceRate, ...
    'logProjectedAreaRatio',state(3,:),'projectedAreaIntegralM',integralA, ...
    'twoPassMassRelativeDifference',twoPass,'tailHessianRelativeIndicator',tailHessianIndicator, ...
    'tailProjectedAreaRatio',exp(tail(3)),'tailCoordinateM',tailQ,'backTime',backTime, ...
    'relativeTolerance',relativeTolerance,'gridRateStep',gridRateStep,'RHSGeometryEvaluations',evaluations, ...
    'interpolationMethod','单调log(psi)上的分段三次Hermite；区间内二分反演', ...
    'tailApproximationIsRigorousBound',false,'stagnationGlobalMatchingSolved',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0);

    function value=geometryRHS(~,y)
        localQ=apex+lengthScale*y(1:2);sample=checked(localQ);
        value=[-sample.graphGradient/(referenceRate*lengthScale); ...
            -trace(sample.graphHessian)/referenceRate;exp(y(3))];
    end
    function value=massRHS(t,y)
        location=deval(geometrySolution,t);sample=checked(apex+lengthScale*location(1:2));
        assert(y(1)>0,'ShockLayer:MassContinuity','稳定推进中的质量系数非正。');
        value=[trace(sample.graphHessian)/referenceRate*y(1)-1;-1/y(1)];
    end
    function sample=checked(location)
        evaluations=evaluations+1;sample=sampler(location);H=sample.graphHessian;
        assert(all(isfinite(H),'all')&&norm(H-H','fro')<1e-10*max(1,norm(H,'fro'))&&min(eig(H))>0, ...
            'ShockLayer:MassConvexity','共同曲面回溯进入非凸/无效区域，不能套旧质量映射。');
        largestRate=max(largestRate,trace(H)/referenceRate);
    end
end
