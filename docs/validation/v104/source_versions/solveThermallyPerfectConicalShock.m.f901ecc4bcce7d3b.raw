function state = solveThermallyPerfectConicalShock(freestream,coneHalfAngleDeg)
%SOLVETHERMALLYPERFECTCONICALSHOCK 固定组分干空气的附体切锥弱解。
% 模块作用：在锥形激波处满足质量、法向动量和总焓守恒，随后积分
% 轴对称无粘锥流方程，搜索使锥面极向速度为零的最小激波角。
% 局部声速、焓和熵均来自同一干空气模型（80~2500 K，低温振动冻结）。
% 仅适用于零攻角直圆锥附体弱支；不描述弯曲前体、脱体激波、
% 三维激波相互作用、化学反应和激波层粘性流动。

%% 1. 校验来流，并建立显式无解结果
validateattributes(coneHalfAngleDeg,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<',90},mfilename,'coneHalfAngleDeg');
required={'Mach','gamma','R','pressure','temperature','density','speed'};
for k=1:numel(required)
    name=required{k};
    if ~isfield(freestream,name)
        error('RocketAeroHeat:ThermalCone:Input','缺少 freestream.%s。',name);
    end
    validateattributes(freestream.(name),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['freestream.' name]);
end
selectAirThermalModel(freestream,struct());
R=freestream.R; T1=freestream.temperature;
U1=freestream.speed; p1=freestream.pressure;
rho1=freestream.density;
if abs(rho1*R*T1/p1-1)>1e-5 || ...
        abs(U1/(freestream.Mach*sqrt(freestream.gamma*R*T1))-1)>1e-5
    error('RocketAeroHeat:ThermalCone:State','来流状态方程或 Mach 与速度不一致。');
end
state=struct('available',false,'attached',false,'branch','weak', ...
    'method','固定组分干空气 NASA7 激波守恒 + 变比热 Taylor-Maccoll 轴对称锥流', ...
    'scope',['零攻角直圆锥、附体弱支、80~2500 K；80~150 K 为振动冻结近似；不含脱体波、' ...
        '三维波系、粘性激波层或高温化学反应'], ...
    'thermalPropertyMode','thermallyPerfect', ...
    'numericalMethod','原始物性 + 去除平凡根的激波能量方程 + ode45 锥面事件 + 近 Mach 波对数括区', ...
    'integrationRelativeTolerance',1e-9,'shootingBracketDeg',[NaN NaN], ...
    'coneHalfAngleDeg',coneHalfAngleDeg,'status','尚未求解');
blank={'shockAngleDeg','upstreamMachThermal','Mach','pressure', ...
    'temperature','density','speed','speedOfSound','totalPressure', ...
    'totalTemperature','upstreamTotalPressureThermal','pressureRatio', ...
    'totalPressureRatio','radialVelocityAtWall','polarVelocityAtWall', ...
    'entropyIncrease','shockMassResidual','shockMomentumResidual', ...
    'energyResidual','isentropicResidual'};
for k=1:numel(blank), state.(blank{k})=NaN; end

%% 2. 缓存来流热力锚点；激波和锥内积分统一使用原始物性
% 弱激波不能混用插值焓与多项式焓：极小的插值误差会淹没真实焓差。
% 比热同样直接求值，避免近 Mach 波处声速误差放大到锥流奇异分母。
persistent cachedKey cachedTable
key=[T1,U1,R];
if isempty(cachedKey) || ~isequal(cachedKey,key)
    table=struct;
    [cp1,h1,s1]=airThermalProperties(T1);
    table.cp1=cp1; table.h1=h1; table.s1=s1;
    table.h0=h1+U1^2/2;
    table.T0=airRecoveryTemperature(T1,U1,1);
    cachedKey=key; cachedTable=table;
end
table=cachedTable;
if ~isfinite(table.T0) || ~isfinite(table.cp1) || table.cp1<=R
    state.status='来流或总焓超出干空气 80~2500 K 热物性域'; return;
end
[~,~,s0]=airThermalProperties(table.T0);
a1=sqrt(table.cp1/(table.cp1-R)*R*T1);
M1=U1/a1;
if M1<=1
    state.status='温度相关声速下的来流不满足超声速条件'; return;
end
p01=p1*exp((s0-table.s1)/R);
state.upstreamMachThermal=M1;
state.upstreamTotalPressureThermal=p01;
state.totalTemperature=table.T0;
if coneHalfAngleDeg==0
    state.available=true; state.attached=true;
    state.status='零半锥角：NASA7 自由来流，无激波';
    state.shockAngleDeg=asind(1/M1);
    state.Mach=M1; state.pressure=p1; state.temperature=T1;
    state.density=rho1; state.speed=U1; state.speedOfSound=a1;
    state.totalPressure=p01; state.pressureRatio=1;
    state.totalPressureRatio=1;
    state.radialVelocityAtWall=U1; state.polarVelocityAtWall=0;
    state.entropyIncrease=0; state.shockMassResidual=0;
    state.shockMomentumResidual=0; state.energyResidual=0;
    state.isentropicResidual=0;
    return;
end

%% 3. 扫描附体弱支并反求激波角；无根时保持 NaN
thetaC=deg2rad(coneHalfAngleDeg);
betaMin=asin(a1/U1);
lower=max(betaMin,thetaC)+1e-10;
upper=pi/2-1e-5;
if lower>=upper
    state.status='半锥角超出附体激波数值范围'; return;
end
fraction=linspace(0,1,36).^2;
angles=lower+(upper-lower)*fraction;
% 小半锥角的弱支紧贴 Mach 波；增加对数节点，不跳过极弱激波根。
nearMach=lower+logspace(-10,-2,30);
angles=unique([angles,nearMach(nearMach<upper)]);
residual=nan(size(angles));
pair=[];
for k=1:numel(angles)
    residual(k)=wallResidual(angles(k),thetaC,freestream,table);
    if k>1 && all(isfinite(residual(k-1:k))) && prod(residual(k-1:k))<=0
        pair=k-1;break; % 首个连续有效括区对应弱支，不再扫描强支。
    end
end
if isempty(pair)
    state.status='未找到附体锥流弱支数值根；需检查模型与数值适用域，不能单独据此判定脱体';
    return;
end
try
    beta=fzero(@(b)wallResidual(b,thetaC,freestream,table), ...
        [angles(pair),angles(pair+1)], ...
        optimset('TolX',1e-11,'Display','off'));
catch
    state.status='锥面无穿透条件的激波角求根失败'; return;
end
[wall,shock,reachedAngle]=integrateCone(beta,thetaC,freestream,table);
if ~shock.valid || any(~isfinite(wall)) || abs(wall(2))>1e-3 || ...
        abs(reachedAngle-thetaC)>1e-8
    state.status='变比热锥流积分未满足锥面无穿透条件'; return;
end

%% 4. 用同一熵常数恢复锥面静压、总压并校验守恒
Tw=wall(3);
[cpw,hw,sw]=airThermalProperties(Tw);
pw=shock.pressure*exp((sw-shock.entropyStandard)/R);
rw=pw/(R*Tw);
Uw=hypot(wall(1),wall(2));
aw=sqrt(cpw/(cpw-R)*R*Tw);
p02=shock.pressure*exp((s0-shock.entropyStandard)/R);
[~,h2,s2]=airThermalProperties(shock.temperature);
entropyIncrease=(s2-table.s1)-R*log(shock.pressure/p1);
massResidual=(shock.density*shock.normalSpeed-rho1*U1*sin(beta))/ ...
    (rho1*U1*sin(beta));
momentumResidual=(shock.pressure+shock.density*shock.normalSpeed^2 ...
    -p1-rho1*(U1*sin(beta))^2)/(p1+rho1*(U1*sin(beta))^2);
energyResidual=(hw+Uw^2/2-table.h0)/max(1,abs(table.h0));
isentropicResidual=(sw-s2)-R*log(pw/shock.pressure);
shockEnergyResidual=(h2+ ...
    (shock.normalSpeed^2+(U1*cos(beta))^2)/2-table.h0)/ ...
    max(1,abs(table.h0));
if ~isfinite(pw) || ~isfinite(p02) || pw<=0 || p02<=pw || ...
        p02>=p01 || entropyIncrease<=0 || ...
        max(abs([massResidual,momentumResidual,energyResidual, ...
        isentropicResidual/max(1,abs(sw)),shockEnergyResidual]))>2e-6
    state.status='激波跃迁或锥内总焓/熵守恒校核失败'; return;
end
state.available=true; state.attached=true;
state.status='NASA7 附体弱锥形激波；锥面无穿透与守恒校核通过';
state.shockAngleDeg=rad2deg(beta);
state.Mach=Uw/aw; state.pressure=pw; state.temperature=Tw;
state.density=rw; state.speed=Uw; state.speedOfSound=aw;
state.totalPressure=p02; state.pressureRatio=pw/p1;
state.totalPressureRatio=p02/p01;
state.radialVelocityAtWall=wall(1);
state.polarVelocityAtWall=wall(2);
state.entropyIncrease=entropyIncrease;
state.shockMassResidual=massResidual;
state.shockMomentumResidual=momentumResidual;
state.energyResidual=energyResidual;
state.isentropicResidual=isentropicResidual;
state.shootingBracketDeg=rad2deg(angles(pair:pair+1));
end

function value=wallResidual(beta,thetaC,freestream,table)
[wall,shock,reachedAngle]=integrateCone(beta,thetaC,freestream,table);
value=NaN;
if shock.valid && all(isfinite(wall))
    if reachedAngle>thetaC
        % 已先遇到无穿透锥面，角度偏差给出正残差；不再向实体内积分。
        value=freestream.speed*(reachedAngle-thetaC);
    else
        value=wall(2);
    end
end
end

function [wall,shock,reachedAngle]=integrateCone(beta,thetaC,free,table)
% 激波后切向速度连续；自适应反向积分处理近 Mach 波处的快速变化。
% 极向速度为零即遇锥面，终止射击，避免越过物理边界后的奇异积分。
shock=normalJump(beta,free,table);
wall=[NaN;NaN;NaN];reachedAngle=NaN;
if ~shock.valid, return; end
wall=[free.speed*cos(beta);-shock.normalSpeed;shock.temperature];
betaMin=asin(sqrt(table.cp1/(table.cp1-free.R)*free.R*free.temperature)/free.speed);
options=odeset('RelTol',1e-9,'AbsTol',[1e-7 1e-7 1e-8], ...
    'InitialStep',min((beta-thetaC)/10,min(1e-5,max(beta-betaMin,1e-10))), ...
    'MaxStep',(beta-thetaC)/40,'Events',@coneWallEvent);
try
    [theta,values,eventAngle]=ode45(@(t,y)coneDerivative(t,y,free.R), ...
        [beta thetaC],wall,options);
    wall=values(end,:)';reachedAngle=theta(end);
    if reachedAngle>thetaC+1e-12 && isempty(eventAngle)
        wall=[NaN;NaN;NaN]; % 中途数值停止不算到达锥面。
    end
catch
    wall=[NaN;NaN;NaN];reachedAngle=NaN;
end
end

function [value,terminal,direction]=coneWallEvent(~,wall)
value=wall(2);terminal=1;direction=0;
end

function derivative=coneDerivative(theta,v,R)
% 轴对称 Euler 方程：dVr/dtheta=Vtheta；连续与极向动量消元给
% dVtheta/dtheta；等熵 dh=cp dT 给温度导数。
Vr=v(1); Vt=v(2); T=v(3);
cp=airThermalProperties(T);
if ~isfinite(cp) || cp<=R || theta<=0
    derivative=[NaN;NaN;NaN]; return;
end
a2=cp/(cp-R)*R*T;
denominator=a2-Vt^2;
if denominator<=1e-9*a2
    derivative=[NaN;NaN;NaN]; return;
end
dVt=(Vt^2*Vr-a2*(2*Vr+Vt/tan(theta)))/denominator;
derivative=[Vt;dVt;-Vt*(dVt+Vr)/cp];
end

function shock=normalJump(beta,free,table)
% 激波法向一维守恒；以密度比求非平凡压缩根。
shock=struct('valid',false,'normalSpeed',NaN,'pressure',NaN, ...
    'temperature',NaN,'density',NaN,'entropyStandard',NaN);
Un1=free.speed*sin(beta);
if Un1<=sqrt(table.cp1/(table.cp1-free.R)*free.R*free.temperature)
    return;
end
lower=1+1e-10; upper=50;
if shockEnergy(lower,Un1,free)<=0 || ...
        shockEnergy(upper,Un1,free)>=0
    return;
end
for k=1:52
    mid=(lower+upper)/2;
    if shockEnergy(mid,Un1,free)>0
        lower=mid;
    else
        upper=mid;
    end
end
r=(lower+upper)/2;
T2=shockTemperature(r,Un1,free);
if T2<80 || T2>2500, return; end
Un2=Un1/r;
p2=free.pressure+free.density*Un1*(Un1-Un2);
[cp2,~,s2]=airThermalProperties(T2);
if ~isfinite(cp2) || cp2<=free.R || ...
        Un2>=sqrt(cp2/(cp2-free.R)*free.R*T2)
    return;
end
shock.valid=true; shock.normalSpeed=Un2;
shock.pressure=p2; shock.temperature=T2;
shock.density=free.density*r; shock.entropyStandard=s2;
end

function residual=shockEnergy(r,Un1,free)
T2=shockTemperature(r,Un1,free);
if T2>2500, residual=1; return; end
if T2<80, residual=-1; return; end
% 消去 r=1 的平凡无激波根及大数焓差相减。下式等于原能量
% 残差乘 r²/(r-1)，不会改变 r>1 的物理解；近等温用比热极限。
cpMean=airEnthalpySecantCp(T2,free.temperature);
residual=cpMean*(Un1^2/free.R-free.temperature*r)-0.5*Un1^2*(r+1);
end

function T2=shockTemperature(r,Un1,free)
T2=free.temperature/r+Un1^2/free.R*(r-1)/r^2;
end
