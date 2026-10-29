function state=shockSurfaceJumpState(upstream,gas,sample,normalOrientation)
%SHOCKSURFACEJUMPSTATE 从同一光滑激波面求三维守恒跃迁及面内热力导数。
% 模块一：法向和导数由曲面一/二阶几何计算，不按面元另选二维转角。
% 模块二：切向速度连续，法向质量/动量/总焓以同一cp积分求非平凡压缩根。
% 模块三：隐式微分压缩根，输出激波面上的真实两向导数和三维梯度。
% 范围：均匀上游、静止光滑激波、理想固定组分气体；不求激波形状或壁面状态。
required={'pressure','temperature','density','velocityVector','method'};
assert(isstruct(upstream)&&isscalar(upstream)&&all(isfield(upstream,required)), ...
    'ShockLayer:Upstream','需要均匀上游p/T/rho、三维速度及来源。');
assert(isstruct(sample)&&isscalar(sample)&&all(isfield(sample, ...
    {'position','firstDerivatives','secondDerivatives'})), ...
    'ShockLayer:Geometry','需要共同激波面的坐标及一/二阶导数。');
assert(isscalar(normalOrientation)&&any(normalOrientation==[-1 1]), ...
    'ShockLayer:Orientation','法向取向须显式为+1或-1。');
assert(strlength(string(upstream.method))>0,'ShockLayer:Source','上游来源不得留空。');
for name={'pressure','temperature','density'}
    validateattributes(upstream.(name{1}),{'double'},{'real','finite','scalar','positive'});
end
validateattributes(upstream.velocityVector,{'double'},{'real','finite','size',[3 1]});
validateattributes(sample.position,{'double'},{'real','finite','size',[3 1]});
validateattributes(sample.firstDerivatives,{'double'},{'real','finite','size',[3 2]});
validateattributes(sample.secondDerivatives,{'double'},{'real','finite','size',[3 3]});
p1=upstream.pressure;T1=upstream.temperature;rho1=upstream.density;u1=upstream.velocityVector;
assert(abs(p1-rho1*gas.R*T1)<=1e-10*p1,'ShockLayer:EOS','上游状态方程不一致。');
t1=surfaceFlowThermodynamicState(T1,p1,gas);H0=t1.enthalpy+dot(u1,u1)/2;
upper=surfaceFlowThermodynamicState(gas.temperatureRangeK(2),p1,gas);
assert(H0<=upper.enthalpy,'ShockLayer:Temperature','总焓超出声明温域，不退回常gamma。');

%% 一、共同曲面的法向、法向导数与法向Mach
A=sample.firstDerivatives;B=sample.secondDerivatives;
lengths=vecnorm(A,2,1);assert(all(lengths>0),'ShockLayer:Geometry','曲面坐标方向退化。');
unitA=A./lengths;metric=unitA'*unitA;
assert(rcond(metric)>1e-10,'ShockLayer:Geometry','曲面坐标线性相关，不能投影修复。');
c=cross(A(:,1),A(:,2));jacobian=norm(c);n=normalOrientation*c/jacobian;
dc=[cross(B(:,1),A(:,2))+cross(A(:,1),B(:,2)), ...
    cross(B(:,2),A(:,2))+cross(A(:,1),B(:,3))];
dn=normalOrientation*(eye(3)-n*n')*dc/jacobian;
w1=dot(u1,n);a1=sqrt(t1.gammaLocal*gas.R*T1);Mn1=w1/a1;
assert(Mn1>1+1e-5,'ShockLayer:NormalMach', ...
    '要求正向法向Mach>1+1e-5；不翻转法向，不将无激波点补作自由流。');

%% 二、密度比为唯一未知；温域外只用单调焓严格确定残差符号
normalH=t1.enthalpy+w1^2/2;
lo=1+min(1e-4,.25*(Mn1-1));hi=1e4;
assert(energy(lo)>0&&energy(hi)<0,'ShockLayer:CompressionRoot','非平凡压缩根未括住。');
for k=1:72
    mid=(lo+hi)/2;if energy(mid)>0,lo=mid;else,hi=mid;end
end
r=(lo+hi)/2;w2=w1/r;rho2=rho1*r;p2=p1+rho1*w1^2*(1-1/r);
T2=p2/(rho2*gas.R);t2=surfaceFlowThermodynamicState(T2,p2,gas);
u2=u1+(w2-w1)*n;speed=norm(u2);S=t2.entropy;
assert(r>1&&p2>p1&&T2>T1&&w2/sqrt(t2.gammaLocal*gas.R*T2)<1, ...
    'ShockLayer:CompressionRoot','未得到压缩且法向亚声速的激波后状态。');
assert(S>t1.entropy,'ShockLayer:Entropy','没有可分辨的物理熵增。');

%% 三、从几何法向导数隐式微分，不对不同点另配激波角或重置总状态
dw1=u1'*dn;
dTdr=-T1/r^2+(w1^2/gas.R)*(2-r)/r^3;
dTdw=2*w1/gas.R*(r-1)/r^2;
Fr=t2.specificHeat*dTdr-w1^2/r^3;
Fw=t2.specificHeat*dTdw+w1*(1/r^2-1);
assert(isfinite(Fr)&&abs(Fr)>1e-10*max(1,w1^2), ...
    'ShockLayer:DerivativeResolution','压缩根导数接近退化，不能输出可靠入口梯度。');
dr=-(Fw/Fr)*dw1;
dT=dTdr*dr+dTdw*dw1;
dp=rho1*(2*w1*(1-1/r)*dw1+(w1^2/r^2)*dr);
drho=rho1*dr;dw2=dw1/r-w1/r^2*dr;
du=n*(dw2-dw1)+(w2-w1)*dn;
dH=t2.specificHeat*dT+u2'*du;
dS=t2.specificHeat/T2*dT-gas.R/p2*dp;
dspeed=(u2'/speed)*du;
grad=@(v)unitA*(metric\(v'./lengths'));
physicalGradients=struct('pressure',grad(dp),'temperature',grad(dT),'density',grad(drho), ...
    'speed',grad(dspeed),'totalEnthalpy',grad(dH),'entropy',grad(dS));
partial=struct('normal',dn,'normalSpeed',dw1,'densityRatio',dr,'pressure',dp, ...
    'temperature',dT,'density',drho,'velocityVector',du,'speed',dspeed,'totalEnthalpy',dH,'entropy',dS);

%% 四、回代完整三维Rankine-Hugoniot；激波后温度与边界层外缘温度分列
mass1=rho1*w1;mass2=rho2*dot(u2,n);
momentum1=mass1*u1+p1*n;momentum2=mass2*u2+p2*n;
residual=struct('mass',(mass2-mass1)/mass1, ...
    'momentumVector',(momentum2-momentum1)/max(p1,rho1*dot(u1,u1)), ...
    'totalEnthalpy',(t2.enthalpy+dot(u2,u2)/2-H0)/max(1,abs(H0)), ...
    'tangentialVelocityVector',(eye(3)-n*n')*(u2-u1)/norm(u1), ...
    'differentiatedTotalEnthalpy',dH/max(1,abs(H0)));
assert(max(abs([residual.mass;residual.momentumVector;residual.totalEnthalpy; ...
    residual.tangentialVelocityVector]))<1e-9,'ShockLayer:Conservation','三维跃迁守恒未通过。');
anchor=struct('available',true,'pressure',p2,'temperature',T2,'density',rho2,'speed',speed, ...
    'method','同一三维激波面法向守恒跃迁，切向速度连续；同cp积分约定');
total=surfaceFlowTotalStateFromAnchor(anchor,gas);
state=struct('available',true,'position',sample.position,'normal',n, ...
    'pressure',p2,'temperature',T2,'density',rho2,'velocityVector',u2,'speed',speed, ...
    'upstreamNormalMach',Mn1,'downstreamNormalMach',w2/sqrt(t2.gammaLocal*gas.R*T2), ...
    'densityRatio',r,'totalEnthalpy',t2.enthalpy+speed^2/2,'entropy',S, ...
    'entropyIncrease',S-t1.entropy,'partialDerivatives',partial, ...
    'physicalGradients',physicalGradients,'inletTotalState',total,'residual',residual, ...
    'method',anchor.method,'stateLocation','给定共同激波面紧后；不是箭体壁面或边界层外缘', ...
    'gas',gas,'shockShapeSolved',false,'bodyMatched',false,'boundaryLayerEdgeAvailable',false, ...
    'pressureOnBodySolved',false,'heatFluxSolved',false,'wallTemperatureSolved',false, ...
    'productionQualified',false,'newProductionSolvedAreaM2',0);

    function value=energy(ratio)
        temp=T1/ratio+(w1^2/gas.R)*(ratio-1)/ratio^2;
        if temp>gas.temperatureRangeK(2)
            % h(temp)>h(Tmax)>=H0>=normalH，故残差严格为正；没有外推物性。
            value=1;return;
        elseif temp<gas.temperatureRangeK(1)
            % temp<Tmin<=T1且ratio>1，焓和法向动能均下降，故严格为负。
            value=-1;return;
        end
        thermal=surfaceFlowThermodynamicState(temp,p1,gas);
        value=thermal.enthalpy+.5*(w1/ratio)^2-normalH;
    end
end
