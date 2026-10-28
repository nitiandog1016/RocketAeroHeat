function table = solveAxisymmetricEntropyShockTable(freestream,geometry,cfg)
%SOLVEAXISYMMETRICENTROPYSHOCKTABLE 内部生成零攻角球锥的近似激波质量表。
% 模块作用：由解析球锥母线、修正牛顿压力及 CR-112180 式65–67，
% 计算与物面弧长 S 关联的 Gamma(S) 与 Psi_sh(S)，供熵层质量匹配。
% 这采用报告式64/66的物面度量近似，不输出真实空间的三维激波面。
% 不读取外部 CFD，不覆盖既有网格、生产外流或热流。
% 输入 geometry.noseRadius、baseRadius、coneAngleDeg、tipX 均为 SI；
% 主轴X、鼻尖+X、零攻角来流-X。球鼻与直锥相切；不含柱段或尾底。
% 物理压力 pe=p_inf+(p_stag-p_inf)*mu_n^2；辅助 pb 则用正激波静压
% p2代替p_stag，只用于式65，二者分别存储，禁止混为同一压力。
% 鼻端按附录B16的B=1正则展开启动，不给Psi添加任意常数。
% 定比热使用原公式；NASA7 是同守恒方程的固定组分热物性推广，
% 不是原报告的平衡空气复现，几何/热流实验精度尚未验证。

%% 1. 输入与范围；只支持已明确给定的零攻角球锥几何
if nargin<3,cfg=struct;end
assert(isstruct(geometry)&&isscalar(geometry)&&isstruct(cfg)&&isscalar(cfg), ...
    'RocketAeroHeat:EntropyShock:Input','geometry 和 cfg 必须为标量结构体。');
for key={'noseRadius','baseRadius','coneAngleDeg','tipX'}
    name=key{1};assert(isfield(geometry,name),'RocketAeroHeat:EntropyShock:Geometry', ...
        '缺少几何量 geometry.%s。',name);
    validateattributes(geometry.(name),{'numeric'},{'real','finite','scalar'});
end
Rn=geometry.noseRadius;Rb=geometry.baseRadius;alpha=geometry.coneAngleDeg;
assert(Rn>0&&Rb>Rn*cosd(alpha)&&alpha>=15&&alpha<=45, ...
    'RocketAeroHeat:EntropyShock:Geometry','要求正球鼻半径、15~45度相切球锥及大于切点的尾部半径。');
for key={'alpha','beta'}
    name=key{1};assert(isfield(freestream,name), ...
        'RocketAeroHeat:EntropyShock:Flight','需明确 freestream.%s，不能推断零攻角。',name);
    validateattributes(freestream.(name),{'numeric'},{'real','finite','scalar'});
end
opt=struct;
if isfield(cfg,'externalFlow')&&isfield(cfg.externalFlow,'entropyShock')
    opt=cfg.externalFlow.entropyShock;
    assert(isstruct(opt)&&isscalar(opt),'RocketAeroHeat:EntropyShock:Options','entropyShock 必须为标量结构体。');
end
count=option('stationCount',161);nPressure=option('pressureTableNodes',257);
startRatio=option('noseStartRatio',1e-3);relTol=option('relativeTolerance',1e-8);
absTol=option('absoluteTolerance',1e-11);
validateattributes(count,{'numeric'},{'scalar','integer','>=',41,'<=',10001});
validateattributes(nPressure,{'numeric'},{'scalar','integer','>=',65,'<=',4097});
validateattributes(startRatio,{'numeric'},{'scalar','finite','>',0,'<=',.02});
validateattributes(relTol,{'numeric'},{'scalar','finite','>=',1e-11,'<=',1e-3});
validateattributes(absTol,{'numeric'},{'scalar','finite','>',0,'<=',1e-5});
table=struct('available',false,'status','尚未获得有效近似激波质量表', ...
    'method','CR-112180 式65–67 物面关联近似 + 附录B16鼻端展开', ...
    'scope','连续流工程假设下零攻角相切球锥15~45度，物面关联近似；未核对Kn及实验精度，不含柱段/底面/分离再附', ...
    'sourceURL','https://ntrs.nasa.gov/api/citations/19730005231/downloads/19730005231.pdf', ...
    'sourcePages','印刷24–28、39–40、42–46；式64–67、A9、B16，B=1', ...
    'geometry',geometry,'upstreamState',freestream, ...
    'upstreamStateSource','当前内部来流输入；同一热物性求正激波与完整总焓', ...
    'shockSurfaceAvailable',false,'defaultSolverEnabled',false, ...
    'geometryVerified',false,'analyticGeometrySpecified',true,'experimentallyValidated',false, ...
    'massMatchingSolved',false,'massCoordinate',[],'shockCrossingAngleDeg',[]);
if hypot(freestream.alpha,freestream.beta)>1e-8
    table.status='该物面关联激波表仅适用于明确的零攻角/零侧滑工况';return;
end

%% 2. 同热物性的法向激波锚点与角度—静压表；不切换物性填补域外
anchor=recoverShockEntropyEdgeState(freestream,90,freestream.pressure,cfg);
table.thermalPropertyMode=anchor.thermalPropertyMode;
table.thermalScope=standaloneAirThermalScope(anchor.thermalPropertyMode);
if ~anchor.available||anchor.upstreamMachPhysical<3
    table.status='需同热物性可用正激波总状态及物理来流Mach>=3';return;
end
p1=freestream.pressure;p2=anchor.shock.pressure;ps=anchor.totalPressure;
Q=freestream.density*freestream.speed^2;
wMax=1-1/anchor.upstreamMachPhysical^2;
wNodes=unique([linspace(0,wMax,nPressure),1e-6,1e-5,1e-4]);
gammaNodes=acosd(sqrt(wNodes));gammaNodes(end)=anchor.machAngleDeg;
jumps=recoverShockEntropyEdgeState(freestream,gammaNodes,p1,cfg);
if ~all(jumps.shock.valid)
    table.status='部分激波角的同物性前向跃迁失败；没有用定比热回退';return;
end
pNodes=jumps.shock.pressure;
if strcmp(anchor.thermalPropertyMode,'caloricallyPerfect')
    D=2*Q/(freestream.gamma+1);derivativeError=0;
    derivativeMethod='定比热正激波解析 D=2*rho_inf*U_inf^2/(gamma+1)';
else
    % Psh=p2-D*cos(Gamma)^2+O(cos^4Gamma)，三档四分之一步长外推到鼻端。
    probe=[1e-3,2.5e-4,6.25e-5];
    near=recoverShockEntropyEdgeState(freestream,acosd(sqrt(probe)),p1,cfg);
    if ~all(near.shock.valid),table.status='鼻端激波压力导数的NASA7校核失败';return;end
    d=(p2-near.shock.pressure)./probe;
    first=(4*d(2)-d(1))/3;D=(4*d(3)-d(2))/3;
    derivativeError=abs(D-first)/abs(D);
    derivativeMethod='同NASA7前向跃迁的cos^2(Gamma)导数，两次Richardson外推';
end
if ~(isfinite(D)&&D>0&&derivativeError<1e-4)
    table.status='正激波压力导数未通过收敛检查';return;
end
a0=(p2-p1)/Q;d0=D/Q;A=a0/(1+2*d0);
table.noseExpansion=struct('A',A,'cosAngleGradient',sqrt(2*A)/Rn, ...
    'shockCurvatureMinusStandOff',Rn/sqrt(2*A), ...
    'pressureDerivative',D,'pressureDerivativeRelativeError',derivativeError, ...
    'pressureDerivativeMethod',derivativeMethod, ...
    'method','附录B16 B=1；sinGamma=1-A*(S/Rn)^2；Psi/(rhoU Rn^2)=t^2/2-(A+1/6)t^4/4', ...
    'startRatio',startRatio);

%% 3. 解析母线及真正的鼻端极限；S/h用于报告的近似度量而非三维激波坐标
join=pi/2-deg2rad(alpha);rt=cosd(alpha);
tEnd=join+(Rb/Rn-rt)/sind(alpha);
% 两段共享切点只存一次，显式写端点，避免 logspace 的末位舍入产生
% 极近的重复站点并破坏严格递增质量坐标。stationCount 是实际总站数。
nNose=ceil((count-1)/2);nCone=count-1-nNose;
noseNodes=logspace(log10(startRatio),log10(join),nNose);
noseNodes(1)=startRatio;noseNodes(end)=join;
coneNodes=linspace(join,tEnd,nCone+1);
t=[0,noseNodes,coneNodes(2:end)]';
q=t(2:end).^2;q0=q(1);
u0=1-A*q0;psi0=q0/2-(A+1/6)*q0^2/4;
if u0<=sqrt(1-wMax)||u0>=1||psi0<=0
    table.status='鼻端展开启动长度超出有效范围';return;
end
settings=odeset('RelTol',relTol,'AbsTol',[absTol,absTol], ...
    'MaxStep',max(.15,tEnd^2/200));
try
    [actual,state]=ode45(@rhs,q,[u0;psi0],settings);
catch problem
    table.status=['近似激波质量积分失败：' problem.message];return;
end
if numel(actual)~=numel(q)||any(~isfinite(state(:)))|| ...
        any(state(:,1)<=sqrt(1-wMax)|state(:,1)>1)||any(diff(state(:,2))<=0)|| ...
        any(diff(state(:,1))>32*eps(max(abs(state(:,1)))))
    table.status='积分状态/质量单调性不适用；没有外推或裁剪激波强度';return;
end
u=[1;state(:,1)];psi=[0;state(:,2)];
[H,mu,x]=meridian(t);
Gamma=asind(u);w=1-u.^2;
shockP=interp1(wNodes,pNodes,w,'pchip');
physical=p1+(ps-p1)*mu.^2;auxiliary=p1+(p2-p1)*mu.^2;
table.available=true;
table.status='近似Gamma/Psi表已内部生成；质量匹配与边界层需另行计算';
table.surfaceArcLength=Rn*t;table.streamlineMetric=Rn*H;
table.surfaceX=x;table.surfaceRadius=Rn*H;
table.shockCrossingAngleDeg=Gamma;
table.massCoordinate=freestream.density*freestream.speed*Rn^2*psi;
table.massCoordinateUnit='kg/s，轴对称每弧度流管；不是整周总质量流率';
table.physicalPressure=physical;table.auxiliaryShockODEPressure=auxiliary;
table.shockStaticPressure=shockP;
table.massCoordinateDerivative=freestream.density*freestream.speed*Rn*H.*u;
table.sinAngleDerivative=zeros(size(t));
table.sinAngleDerivative(2:end)=(auxiliary(2:end)-shockP(2:end))./Q.* ...
    H(2:end)./psi(2:end)/Rn;
table.stagnationPressure=ps;table.normalShockStaticPressure=p2;
table.massCoordinateSource='CR-112180 式66：物面h=r和弧长S的近似流管积分；内部ODE数值解';
table.shockShapeSource=['CR-112180 式65/67 + B16鼻端闭合；物面关联Gamma，' ...
    '不是真实空间shock mesh；' derivativeMethod];
table.physicalPressureSource='同热物性正激波驻点总压 + 解析相切球锥修正牛顿压力';
table.streamlineMetricSource='零攻角解析轴对称母线h=r(S)，每弧度流管';
table.auxiliaryPressureScope='只用于激波角ODE；不可标为物面压力或边界层外缘压力';
table.numerics=struct('stationCount',numel(t),'pressureTableNodes',numel(wNodes), ...
    'relativeTolerance',relTol,'absoluteTolerance',absTol,'noseStartRatio',startRatio, ...
    'pressureInterpolation','同热物性前向跃迁Psh(cos²Gamma)的PCHIP数值插值');
for key={'physicalPressure','streamlineMetric','massCoordinate','shockCrossingAngleDeg'}
    name=key{1};sourceName=[name 'Source'];
    if strcmp(name,'shockCrossingAngleDeg'),sourceName='shockShapeSource';end
    table.fieldMetadata.(name)=struct('valueKind','computed','source',table.(sourceName), ...
        'scope',table.scope,'valid',true(size(t)));
end

    function value=option(name,fallback)
        value=fallback;if isfield(opt,name),value=opt.(name);end
    end
    function [H,mu,x]=meridian(tt)
        spherical=tt<=join;H=sin(min(tt,join));mu=cos(min(tt,join));
        x=geometry.tipX-Rn*(1-cos(min(tt,join)));
        downstream=tt(~spherical)-join;
        H(~spherical)=rt+downstream*sind(alpha);mu(~spherical)=sind(alpha);
        x(~spherical)=geometry.tipX-Rn*(1-sind(alpha))-Rn*downstream*cosd(alpha);
    end
    function change=rhs(qq,y)
        % q=(S/Rn)^2。有限展开在q0启动，未改变物理鼻端Psi(0)=0。
        tt=sqrt(qq);[H,mu]=meridian(tt);
        ww=1-y(1)^2;
        if y(2)<=0||ww<0||ww>wMax||y(1)<=0
            error('RocketAeroHeat:EntropyShock:Domain','ODE试探状态超出物理激波/正质量域。');
        end
        psh=interp1(wNodes,pNodes,ww,'pchip');pb=p1+(p2-p1)*mu^2;
        change=[(pb-psh)/Q*H/(2*tt*y(2));y(1)*H/(2*tt)];
    end
end
