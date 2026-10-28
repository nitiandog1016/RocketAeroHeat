function out=nonsimilarLaminarResidual(current,history,gas,forcing)
%NONSIMILARLAMINARRESIDUAL 物理弧长/壁法向坐标的非相似层流单站残差。
% 模块作用：保留真实上游剖面对当前站的影响；不将局部相似BVP改名。
% 来源：NASA TP-3271式1--4的层流极限。s沿物面，n朝流体，r为
% 物面轴对称流管半径。忽略O(delta/R)法向曲率/度量修正，故需薄层。
% 连续：d_s(r*rho*u)+d_n(r*rho*v)=0；壁面无吹吸。
% 动量：rho*u*u_s+rho*v*u_n+p_s-d_n(mu*u_n)=forcing.momentum。
% 能量：rho*u*H_s+rho*v*H_n-d_n(k*T_n+mu*u*u_n)=forcing.energy，
% H=h(T)+u^2/2；NASA7的焓零点不被抹除，k=mu*cp/Pr。
% n网格固定于物理空间；相同eta编号不能作为相同n。v由连续积分，
% 不另外给外缘v=0。当前站/上游站须含s,r,p,n,u,T；当前另给p_s、
% Ue、Te、He、Tw和逐字段sources。上游必须是真实剖面，非厚度标量。
% 一个上游站使用后退一阶；两个上游站使用不等步长后退二阶。
% forcing仅用于制造解；给定力/热源不能当成无源球鼻计算。
% 返回数值空间残差/壁导数，不声称已解区域、稳定性或实验精度。
% 正式迁移身份：受限条件研究模块；核心数值保持v61冻结草稿原式。
if nargin<4,forcing=struct;end
checkGas(gas);checkStation(current,true);
if ~isstruct(history)||isempty(history)||numel(history)>2
    error('RocketAeroHeat:NonSimilar:History','需要1或2个真实上游剖面，按最近到最早排列。');
end
n=double(current.n(:));u=double(current.u(:));T=double(current.T(:));N=numel(n);
if N<5||numel(u)~=N||numel(T)~=N||n(1)~=0||any(~isfinite(n))||any(diff(n)<=0)
    error('RocketAeroHeat:NonSimilar:Grid','需要从n=0开始、至少5点的严格递增物理n网格。');
end
profiles=[{current},num2cell(history(:)')];
for j=1:numel(profiles)
    entry=profiles{j};checkStation(entry,j==1);
    if ~isequal(double(entry.n(:)),n)
        error('RocketAeroHeat:NonSimilar:Grid','各站必须使用完全相同的物理n网格；禁止隐式外推。');
    end
    if numel(entry.u)~=N||numel(entry.T)~=N
        error('RocketAeroHeat:NonSimilar:Size','各站u/T必须与物理n网格一致。');
    end
end
s=cellfun(@(z)double(z.s),profiles);ds=s(1)-s(2);
if any(diff(s)>=0)
    error('RocketAeroHeat:NonSimilar:History','历史弧长必须严格小于当前站，且最近站在前。');
end
if numel(history)==1
    weights=[1,-1]/ds;order=1;
else
    oldStep=s(2)-s(3);
    weights=[(2*ds+oldStep)/(ds*(ds+oldStep)), ...
        -(ds+oldStep)/(ds*oldStep),ds/(oldStep*(ds+oldStep))];order=2;
end
U=zeros(N,numel(profiles));HH=U;density=U;mass=U;
for j=1:numel(profiles)
    entry=profiles{j};[~,h]=thermo(double(entry.T(:)),gas);
    U(:,j)=double(entry.u(:));
    HH(:,j)=h+.5*U(:,j).^2;
    density(:,j)=entry.p./(gas.R*double(entry.T(:)));
    mass(:,j)=entry.r*density(:,j).*U(:,j);
end
[cp,h]=thermo(T,gas);rho=density(:,1);H=HH(:,1);
if strcmp(gas.transport,'sutherland')
    mu=1.458e-6*T.^(3/2)./(T+110.4);
else
    mu=gas.constantViscosity*ones(N,1);
end
k=mu.*cp/gas.Prandtl;
massS=mass*weights';uS=U*weights';HS=HH*weights';
normalMass=-cumtrapz(n,massS);v=normalMass./(current.r*rho);
dy=diff(n);dual=.5*(n(3:end)-n(1:end-2));
left=n(2:end-1)-n(1:end-2);right=n(3:end)-n(2:end-1);
uN=derivative(u);HN=derivative(H);
muFace=.5*(mu(1:end-1)+mu(2:end));kFace=.5*(k(1:end-1)+k(2:end));
shear=muFace.*diff(u)./dy;
energyFlux=kFace.*diff(T)./dy+shear.*(.5*(u(1:end-1)+u(2:end)));
momDiff=diff(shear)./dual;energyDiff=diff(energyFlux)./dual;
ids=2:N-1;
momConv=rho(ids).*u(ids).*uS(ids)+rho(ids).*v(ids).*uN;
energyConv=rho(ids).*u(ids).*HS(ids)+rho(ids).*v(ids).*HN;
fm=source('momentum');fe=source('energy');
momentum=momConv+current.p_s-momDiff-fm(ids);
energy=energyConv-energyDiff-fe(ids);
% 模块作用：壁面导数用真实n前三点的二次插值，不用单元中心距离。
wn=derivativeWeights(n(1:3),0);uWallSlope=wn*u(1:3);TWallSlope=wn*T(1:3);
q=k(1)*TWallSlope;tau=mu(1)*uWallSlope;
[~,he]=thermo(current.Te,gas);edgeH=he+.5*current.Ue^2;
if abs(edgeH-current.He)>1e-7*max([abs(he),.5*current.Ue^2,1])
    error('RocketAeroHeat:NonSimilar:EdgeEnergy','外缘Ue/Te/He未使用同一物性或总焓。');
end
rhoE=current.p/(gas.R*current.Te);[cpW,~]=thermo(current.Tw,gas);
if strcmp(gas.transport,'sutherland'),muW=1.458e-6*current.Tw^1.5/(current.Tw+110.4);
else,muW=gas.constantViscosity;end
speedScale=max(current.Ue,1);tempScale=max([current.Te,current.Tw,1]);
momScale=max([rhoE*speedScale^2/max(current.s,ds),abs(current.p_s),muW*speedScale/n(end)^2,1]);
energyScale=max([rhoE*speedScale*(cpW*tempScale+.5*speedScale^2)/max(current.s,ds),muW*cpW/gas.Prandtl*tempScale/n(end)^2,1]);
finiteAudit=[weights(:);massS;uS;HS;rho;H;h;cp;mu;k;normalMass;v; ...
    momentum;energy;shear;energyFlux;q;tau;momScale;energyScale];
if any(~isfinite(finiteAudit))
    error('RocketAeroHeat:NonSimilar:NumericalOverflow','流向/法向差分或物性运算产生非有限量，不能作为残差/解。');
end
out=struct('momentum',momentum,'energy',energy, ...
    'scaledResidual',[momentum/momScale;energy/energyScale], ...
    'momentumScale',momScale,'energyScale',energyScale, ...
    'wallBoundaryResidual',[u(1);T(1)-current.Tw], ...
    'edgeBoundaryResidual',[u(end)-current.Ue;T(end)-current.Te], ...
    'density',rho,'enthalpy',h,'totalEnthalpy',H,'heatCapacity',cp, ...
    'viscosity',mu,'conductivity',k,'normalMassFluxPerRadian',normalMass, ...
    'streamwiseMetricMassDerivative',massS,'normalVelocity',v, ...
    'velocityStreamwiseDerivative',uS,'totalEnthalpyStreamwiseDerivative',HS, ...
    'velocityNormalDerivativeInterior',uN,'totalEnthalpyNormalDerivativeInterior',HN, ...
    'momentumConvection',momConv,'energyConvection',energyConv, ...
    'momentumDiffusion',momDiff,'energyDiffusion',energyDiff, ...
    'shearAtNormalFaces',shear,'totalEnthalpyDiffusiveFluxAtNormalFaces',energyFlux, ...
    'wallHeatFluxIntoWall',q,'wallShearStress',tau, ...
    'wallTemperatureInput',current.Tw,'wallTemperatureSolved',NaN, ...
    'streamwiseDifferenceOrder',order,'streamwiseWeights',weights, ...
    'normalCoordinate','物理n [m]，从壁朝流体；不是eta', ...
    'curvatureScope','使用物面r(s)的薄层轴对称度量；忽略O(delta/R)法向曲率项', ...
    'normalTopVelocityPrescribed',false,'forcingPresent',any(fm~=0)|any(fe~=0), ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'method','TP-3271物理薄层方程，保留流向后退差分与连续积分法向质量流', ...
    'thermalPropertyMode',gas.thermalModel,'sourceEquations','NASA TP-3271式1--4，层流极限');
out.continuityIntervalBalance=diff(normalMass)+.5*dy.*(massS(1:end-1)+massS(2:end));
out.integratedNormalMassBalance=normalMass(end)+trapz(n,massS);
% 模块作用：独立保守积分账；非保守点式差分不满足精确离散乘积法则。
% 必须报告此账并做网格/步长收敛，不能用Newton残差代替能量守恒。
wt=derivativeWeights(n(end-2:end),n(end));
topShear=mu(end)*(wt*u(end-2:end));
topEnergyFlux=k(end)*(wt*T(end-2:end))+u(end)*topShear;
energyStorage=trapz(n,(mass.*HH)*weights');
energyBoundary=normalMass(end)*H(end)-current.r*topEnergyFlux+current.r*q;
energyForcing=current.r*trapz(n,fe);
out.conservativeEnergyIntegral=struct('streamwiseFluxDerivativeWPerM',energyStorage, ...
    'normalBoundaryBalanceWPerM',energyBoundary,'forcingWPerM',energyForcing, ...
    'defectWPerM',energyStorage+energyBoundary-energyForcing, ...
    'relativeDefect',abs(energyStorage+energyBoundary-energyForcing)/ ...
        max([abs(energyStorage),abs(energyBoundary),abs(energyForcing),1]), ...
    'scope','独立保守积分诊断；不等于点式Newton残差，需随ds/dn收敛');
if any(~isfinite([energyStorage,energyBoundary,energyForcing, ...
        out.conservativeEnergyIntegral.defectWPerM,out.conservativeEnergyIntegral.relativeDefect]))
    error('RocketAeroHeat:NonSimilar:NumericalOverflow','保守能量积分发生非有限运算。');
end
out.topGradients=struct('velocityPerS',wt*u(end-2:end), ...
    'temperatureKPerM',wt*T(end-2:end),'energyFluxWPerM2',topEnergyFlux);

    function d=derivative(f)
        d=-right./(left.*(left+right)).*f(1:end-2)+ ...
            (right-left)./(left.*right).*f(2:end-1)+ ...
            left./(right.*(left+right)).*f(3:end);
    end
    function f=source(key)
        f=zeros(N,1);
        if isfield(forcing,key)
            if ~isnumeric(forcing.(key))||~isreal(forcing.(key))
                error('RocketAeroHeat:NonSimilar:Forcing','制造源必须是实数数值向量。');
            end
            f=double(forcing.(key)(:));
            if ~isreal(f)||numel(f)~=N||any(~isfinite(f))
                error('RocketAeroHeat:NonSimilar:Forcing','制造源项必须为有限N点向量。');
            end
        end
    end
end

function checkGas(g)
required={'R','gamma','thermalModel','Prandtl','transport'};
for i=1:numel(required)
    if ~isfield(g,required{i}),error('RocketAeroHeat:NonSimilar:Gas','缺少gas.%s。',required{i});end
end
if ~isnumeric(g.R)||~isreal(g.R)||~isscalar(g.R)||~isfinite(g.R)||g.R<=0|| ...
        ~isnumeric(g.gamma)||~isreal(g.gamma)||~isscalar(g.gamma)||~isfinite(g.gamma)||g.gamma<=1|| ...
        ~isnumeric(g.Prandtl)||~isreal(g.Prandtl)||~isscalar(g.Prandtl)||~isfinite(g.Prandtl)||g.Prandtl<=0|| ...
        ~ismember(g.thermalModel,{'thermallyPerfect','caloricallyPerfect'})|| ...
        ~ismember(g.transport,{'sutherland','constant'})
    error('RocketAeroHeat:NonSimilar:Gas','热物性/输运合同不合法。');
end
if strcmp(g.thermalModel,'thermallyPerfect')&&(abs(g.R-287.05)>5||abs(g.gamma-1.4)>.05)
    error('RocketAeroHeat:NonSimilar:Gas','NASA7仅接受项目固定组分干空气。');
end
if strcmp(g.transport,'constant')&&(~isfield(g,'constantViscosity')|| ...
        ~isnumeric(g.constantViscosity)||~isreal(g.constantViscosity)|| ...
        ~isscalar(g.constantViscosity)||~isfinite(g.constantViscosity)||g.constantViscosity<=0)
    error('RocketAeroHeat:NonSimilar:Gas','常物性测试必须明确给出正黏度。');
end
end

function checkStation(c,isCurrent)
required={'s','r','p','n','u','T','sources'};
if isCurrent,required=[required,{'p_s','Ue','Te','He','Tw'}];end
for i=1:numel(required)
    if ~isfield(c,required{i}),error('RocketAeroHeat:NonSimilar:Input','剖面缺少%s。',required{i});end
end
for key={'s','r','p'}
    z=c.(key{1});
    if ~isnumeric(z)||~isreal(z)||~isscalar(z)||~isfinite(z)||z<=0
        error('RocketAeroHeat:NonSimilar:Input','s/r/p须为有限正标量；s=0须另用解析驻点启动。');
    end
end
if ~isnumeric(c.n)||~isreal(c.n)||~isvector(c.n)|| ...
        ~isnumeric(c.u)||~isreal(c.u)||~isvector(c.u)|| ...
        ~isnumeric(c.T)||~isreal(c.T)||~isvector(c.T)|| ...
        any(c.u(:)<0)||any(~isfinite(c.u(:)))||any(~isfinite(c.T(:)))||any(c.T(:)<80|c.T(:)>2500)
    error('RocketAeroHeat:NonSimilar:Property','剖面非有限或温度超出80--2500K，禁止裁温。');
end
if isCurrent
    for key={'p_s','Ue','Te','He','Tw'}
        z=c.(key{1});
        if ~isnumeric(z)||~isreal(z)||~isscalar(z)||~isfinite(z)
            error('RocketAeroHeat:NonSimilar:Input','当前站边界量必须为有限标量。');
        end
    end
    if c.Ue<=0||c.Te<80||c.Te>2500||c.Tw<80||c.Tw>2500
        error('RocketAeroHeat:NonSimilar:Property','边界速度/温度越域，不作静默修正。');
    end
end
if ~isstruct(c.sources)||~isscalar(c.sources)
    error('RocketAeroHeat:NonSimilar:Source','各站须保留逐字段来源。');
end
for key={'geometry','pressure','profile'}
    if ~isfield(c.sources,key{1})||strlength(string(c.sources.(key{1})))==0
        error('RocketAeroHeat:NonSimilar:Source','缺少%s来源。',key{1});
    end
end
if isCurrent
    for key={'edge','wallTemperature'}
        if ~isfield(c.sources,key{1})||strlength(string(c.sources.(key{1})))==0
            error('RocketAeroHeat:NonSimilar:Source','缺少%s来源。',key{1});
        end
    end
end
end

function [cp,h]=thermo(T,g)
if strcmp(g.thermalModel,'thermallyPerfect')
    [cp,h]=airThermalProperties(T);
else
    cp=g.gamma*g.R/(g.gamma-1)+zeros(size(T));h=cp.*T;
end
if any(~isfinite(cp(:)))||any(~isfinite(h(:)))||any(cp(:)<=0)
    error('RocketAeroHeat:NonSimilar:Property','原热物性没有可用值。');
end
end

function w=derivativeWeights(n,target)
L=n(end)-n(1);z=(n-target)/L;
w=([ones(1,3);z';(z.^2)']\[0;1;0])'/L;
end
