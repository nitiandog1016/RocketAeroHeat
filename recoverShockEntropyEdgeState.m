function edge = recoverShockEntropyEdgeState(freestream,shockCrossingAngleDeg,physicalPressure,cfg)
%RECOVERSHOCKENTROPYEDGESTATE 已知激波穿越角和物理静压的独立外缘恢复内核。
% 模块作用：将 Gamma（激波与上游速度的夹角）转成法向速度，只对
% 法向分量执行 Rankine-Hugoniot 跃迁，切向速度保持连续；由完整速度
% 的总焓与激波后熵恢复局部 p、T、rho、U、Mach、p0、T0。
% 输入物理静压和穿越角并非本内核求解；Gamma 不是壁面转角，不能
% 当成切楔转角。几何识别、激波形状、流管度量、质量匹配均在调用方。
% 当前独立模块不接入 solveExternalFlow，不产生三维速度方向或热流。
% 来源：NASA CR-112180 (1972), 式68–69 的穿越激波熵恢复；完全气体
% 跃迁见 NACA Report 1135。默认固定组分 NASA7 干空气 80~2500 K，
% 域外保持 NaN，不静默切换定比热。显式定比热模式沿用项目配置。
%
% 调用：edge=recoverShockEntropyEdgeState(fs,GammaDeg,p,cfg)。角度和
% 压力可为同尺寸数组，其中一个标量时广播。可在 cfg.externalFlow
% 中设置 shockCrossingAngleSource、physicalPressureSource 来源文字。

%% 1. 校验上游状态与数组契约；结构错误报错，逐点域外返回 NaN
if nargin<4, cfg=struct; end
if ~isstruct(freestream) || ~isscalar(freestream) || ...
        ~isstruct(cfg) || ~isscalar(cfg)
    error('RocketAeroHeat:ShockEntropyEdge:Input','上游状态和 cfg 必须为标量结构体。');
end
required={'Mach','gamma','R','pressure','temperature','density','speed'};
for j=1:numel(required)
    key=required{j};
    if ~isfield(freestream,key)
        error('RocketAeroHeat:ShockEntropyEdge:Input','缺少 freestream.%s。',key);
    end
    validateattributes(freestream.(key),{'numeric'}, ...
        {'real','finite','scalar','positive'},mfilename,['freestream.' key]);
end
if freestream.gamma<=1
    error('RocketAeroHeat:ShockEntropyEdge:State','上游 gamma 必须大于 1。');
end
R=freestream.R; T1=freestream.temperature; p1=freestream.pressure;
rho1=freestream.density; U1=freestream.speed;
if abs(rho1*R*T1/p1-1)>1e-5 || ...
        abs(U1/(freestream.Mach*sqrt(freestream.gamma*R*T1))-1)>1e-5
    error('RocketAeroHeat:ShockEntropyEdge:State','上游密度、名义 Mach 与速度不一致。');
end
validateattributes(shockCrossingAngleDeg,{'numeric'},{'real','nonempty'});
validateattributes(physicalPressure,{'numeric'},{'real','nonempty'});
angle=double(shockCrossingAngleDeg); pInput=double(physicalPressure);
if isscalar(angle), angle=repmat(angle,size(pInput)); end
if isscalar(pInput), pInput=repmat(pInput,size(angle)); end
if ~isequal(size(angle),size(pInput))
    error('RocketAeroHeat:ShockEntropyEdge:Size','穿越角和物理静压须同尺寸或其中一个为标量。');
end
[mode,thermalScope]=selectAirThermalModel(freestream,cfg);
thermalScope=standaloneAirThermalScope(mode);
pressureSource=inputSource(cfg,'physicalPressureSource','输入物理静压；本内核未计算压力分布');
angleSource=inputSource(cfg,'shockCrossingAngleSource','输入激波穿越角；本内核未计算激波几何或质量匹配');
blank=nan(size(angle));
names={'pressure','temperature','density','speed','Mach','speedOfSound', ...
    'totalPressure','totalTemperature','entropyIncrease','energyResidual','entropyResidual'};
edge=struct('available',false,'algebraicallyValid',false(size(angle)), ...
    'shockCrossingAngleDegInput',angle,'physicalPressureInput',pInput, ...
    'statusCode',ones(size(angle),'uint8'), ...
    'statusNames',{{'已计算（独立热力恢复）','激波穿越角或上游 Mach 不适用', ...
        '输入静压非正、非有限或超过同模型总压','热物性范围外', ...
        '前向激波求根失败或强度低于数值精度','守恒或熵校核失败'}}, ...
    'thermalPropertyMode',mode,'thermalScope',thermalScope, ...
    'method',['穿越角的前向法向激波 + 切向速度连续 + 同总焓/熵外缘恢复；' mode], ...
    'scope',['独立单点/批量热力内核；Gamma 为激波角；不验证几何、波形、' ...
        '流管质量匹配或压力来源；不求三维流向、分离再附或热流'], ...
    'pressureSource',pressureSource,'shockCrossingAngleSource',angleSource, ...
    'sourceURL',{{'https://ntrs.nasa.gov/api/citations/19730005231/downloads/19730005231.pdf', ...
        'https://www.grc.nasa.gov/www/k-12/airplane/Images/naca1135.pdf'}}, ...
    'geometryVerified',false,'massMatchingSolved',false, ...
    'status','尚无通过校核的独立恢复状态');
for j=1:numel(names), edge.(names{j})=blank; end
shockNames={'pressure','temperature','density','speed','totalPressure', ...
    'normalSpeedBefore','normalSpeedAfter','tangentialSpeedBefore','tangentialSpeedAfter', ...
    'normalMachBefore','normalMachAfter','deflectionDeg','entropyIncrease', ...
    'massResidual','normalMomentumResidual','energyResidual','entropyIdentityResidual'};
shock=struct('valid',false(size(angle)),'machWaveLimit',false(size(angle)), ...
    'normalKernelMethod',{repmat({''},size(angle))});
for j=1:numel(shockNames), shock.(shockNames{j})=blank; end
edge.shock=shock;
edge.fieldMetadata=fieldMetadata(names,edge.method,edge.scope,pressureSource);
edge=syncFieldMetadata(edge,names);

%% 2. 同一热物性的上游声速、完整流动总焓/总温/总压
if strcmp(mode,'thermallyPerfect')
    [cp1,h1,s1]=airThermalProperties(T1);
    T0=airRecoveryTemperature(T1,U1,1);
    if ~isfinite(cp1) || cp1<=R || ~isfinite(T0)
        edge.statusCode(:)=3; edge.status='上游静温或完整速度总焓超出干空气物性域';
        edge=syncFieldMetadata(edge,names); return;
    end
    gamma1=cp1/(cp1-R); [~,h0,s0]=airThermalProperties(T0);
    p01=p1*exp((s0-s1)/R);
else
    gamma1=freestream.gamma; cp1=gamma1*R/(gamma1-1);
    h1=cp1*T1; h0=h1+U1^2/2; T0=h0/cp1;
    s1=cp1*log(T1); s0=cp1*log(T0);
    p01=p1*(T0/T1)^(gamma1/(gamma1-1));
end
a1=sqrt(gamma1*R*T1); M1=U1/a1;
edge.upstreamMachPhysical=M1; edge.upstreamTotalPressure=p01;
edge.upstreamTotalTemperature=T0; edge.upstreamTotalEnthalpy=h0;
if M1<=1 || ~isfinite(p01) || p01<=0
    edge.status='同一热物性下的上游流动非超声速或总状态不可用';
    edge=syncFieldMetadata(edge,names); return;
end
mu=asind(1/M1); edge.machAngleDeg=mu;
geometric=isfinite(angle) & angle>=mu-1e-11 & angle<=90 & angle>0;
edge.statusCode(geometric)=4;

%% 3. 已知穿越角的前向跃迁：Gamma 只决定法向分量，保留切向能量
% 相同角度只计算一次。总温与总压必须由完整速度恢复，不能拿法向核
% 的法向总温/总压当作斜激波后全流总状态。
[uniqueAngle,~,mapping]=unique(angle(geometric));
ids=find(geometric);
for k=1:numel(uniqueAngle)
    beta=uniqueAngle(k); group=ids(mapping==k);
    Un1=U1*sind(beta); Ut1=U1*cosd(beta); Mn1=Un1/a1;
    machWave=abs(Mn1-1)<=64*eps;
    if Mn1<1 && ~machWave, edge.statusCode(group)=1; continue; end
    if machWave
        p2=p1; T2=T1; rho2=rho1; Un2=Un1;
        normalMethod='解析 Mach 波极限：零强度跃迁';
    else
        [jump,normalMethod]=forwardNormalJump(freestream,Un1,gamma1,mode);
        if ~jump.available, continue; end
        p2=jump.postShockPressure; T2=jump.postShockTemperature;
        rho2=jump.postShockDensity; Un2=jump.postShockSpeed;
    end
    U2=hypot(Un2,Ut1);
    if strcmp(mode,'thermallyPerfect')
        [cp2,h2,s2]=airThermalProperties(T2); gamma2=cp2/(cp2-R);
        p02=p2*exp((s0-s2)/R);
    else
        cp2=cp1; h2=cp2*T2; s2=cp2*log(T2); gamma2=gamma1;
        p02=p2*(1+U2^2/(2*cp2*T2))^(gamma2/(gamma2-1));
    end
    a2=sqrt(gamma2*R*T2);
    ds=(s2-s1)-R*log(p2/p1);
    residual=[(rho2*Un2-rho1*Un1)/(rho1*Un1), ...
        (p2+rho2*Un2^2-p1-rho1*Un1^2)/(p1+rho1*Un1^2), ...
        (h2+U2^2/2-h0)/max(1,abs(h0)), ...
        ds-R*log(p01/p02)];
    if machWave
        p02=p01; ds=0; residual(:)=0;
    end
    if ~all(isfinite([p02,a2,ds,residual])) || p02<=0 || ...
            (machWave && p02~=p01) || ...
            (~machWave && (ds<=0 || p02>=p01 || Un2/a2>=1)) || ...
            max(abs(residual(1:3)))>1e-8 || abs(residual(4))>1e-6
        edge.statusCode(group)=5; continue;
    end
    edge.shock.valid(group)=true; edge.shock.machWaveLimit(group)=machWave;
    edge.shock.normalKernelMethod(group)={normalMethod};
    values=[p2,T2,rho2,U2,p02,Un1,Un2,Ut1,Ut1,Mn1,Un2/a2, ...
        beta-atan2d(Un2,Ut1),ds,residual];
    for j=1:numel(shockNames), edge.shock.(shockNames{j})(group)=values(j); end
end

%% 4. 输入物理静压与穿越熵恢复外缘；给定量保留来源，域外不补零
validShock=edge.shock.valid;
p0=edge.shock.totalPressure;
validP=validShock & isfinite(pInput) & pInput>0 & pInput<=p0;
edge.statusCode(validShock & ~validP)=2;
if any(validP(:))
    if strcmp(mode,'thermallyPerfect')
        % 旧恢复核内部按列向量求根；展平输入后再恢复原数组尺寸，
        % 避免行向量静压与列向量温度隐式广播生成二维密度数组。
        recovered=recoverThermallyPerfectEdgeState(pInput(:),p0(:),T0,R);
        valid=validP & reshape(recovered.valid,size(angle));
        Te=reshape(recovered.temperature,size(angle));
        rhoe=reshape(recovered.density,size(angle));
        Ue=reshape(recovered.speed,size(angle)); Me=reshape(recovered.Mach,size(angle));
        ae=reshape(recovered.speedOfSound,size(angle));
        [~,he,se]=airThermalProperties(Te);
        entropyResidual=s0-se-R*log(p0./pInput);
    else
        valid=validP;
        Te=T0*(pInput./p0).^((gamma1-1)/gamma1);
        Ue=sqrt(max(0,2*cp1*(T0-Te)));
        rhoe=pInput./(R*Te); ae=sqrt(gamma1*R*Te); Me=Ue./ae;
        he=cp1*Te; se=cp1*log(Te);
        entropyResidual=s0-se-R*log(p0./pInput);
    end
    edge.statusCode(validP & ~valid)=3;
    atStagnation=valid & pInput==p0;
    Te(atStagnation)=T0; Ue(atStagnation)=0; Me(atStagnation)=0;
    rhoe(atStagnation)=p0(atStagnation)/(R*T0);
    he(atStagnation)=h0; entropyResidual(atStagnation)=0;
    energyResidual=(he+Ue.^2/2-h0)/max(1,abs(h0));
    good=valid & isfinite(energyResidual) & isfinite(entropyResidual) & ...
        abs(energyResidual)<=1e-8 & abs(entropyResidual)<=1e-6;
    edge.statusCode(valid & ~good)=5;
    edge.algebraicallyValid(good)=true; edge.statusCode(good)=0;
    fields={pInput,Te,rhoe,Ue,Me,ae,p0,repmat(T0,size(angle)), ...
        edge.shock.entropyIncrease,energyResidual,entropyResidual};
    for j=1:numel(names), edge.(names{j})(good)=fields{j}(good); end
end
edge=syncFieldMetadata(edge,names);
edge.available=any(edge.algebraicallyValid(:));
edge.status=sprintf(['独立恢复 %d/%d 个状态；物理静压与穿越角是输入；' ...
    '未求激波形状/质量匹配，未接入默认三维外流'],nnz(edge.algebraicallyValid),numel(angle));
end

function edge=syncFieldMetadata(edge,names)
% 每个出口同步逐字段有效掩码与状态；物性越界早退也不得留下初始状态。
for j=1:numel(names)
    edge.fieldMetadata.(names{j}).valid=edge.algebraicallyValid;
    edge.fieldMetadata.(names{j}).statusCode=edge.statusCode;
end
end

function [state,method]=forwardNormalJump(fs,normalSpeed,gammaPhysical,mode)
% 现有核中的 gamma 同时作干空气识别及名义 Mach/速度换算。
% 将 speed 换成法向速度并重算名义 Mach，不能只换 Mach 或只换 speed。
input=fs; input.speed=normalSpeed;
if strcmp(mode,'caloricallyPerfect')
    input.Mach=normalSpeed/sqrt(fs.gamma*fs.R*fs.temperature);
    cp=fs.gamma*fs.R/(fs.gamma-1);
    input.totalTemperature=fs.temperature+normalSpeed^2/(2*cp);
    input.totalPressure=fs.pressure*(input.totalTemperature/fs.temperature)^(fs.gamma/(fs.gamma-1));
    state=solveNormalShock(input);
    method='复用完全气体法向 Rankine-Hugoniot 核（仅法向速度）';
elseif gammaPhysical>1.35+1e-10 && gammaPhysical<1.45-1e-10
    input.gamma=gammaPhysical;
    input.Mach=normalSpeed/sqrt(gammaPhysical*fs.R*fs.temperature);
    state=solveThermallyPerfectNormalShock(input);
    method='复用 NASA7 法向守恒核（真实声速、仅法向速度）';
else
    % 高静温干空气的真实 gamma 可小于旧核气体识别阈值1.35。
    % 保留同一 NASA7 物性和密度比守恒方程，避免伪造名义 Mach 或
    % 使用错误声速扩大/缩小真实 Mach 波适用域；不退回定比热。
    state=thermalDensityJump(input,gammaPhysical);
    method='NASA7 前向密度比守恒求根（旧核 gamma 识别接口不适用）';
end
end

function state=thermalDensityJump(fs,gammaPhysical)
% 密度比 r>1 排除平凡根。式由质量/法向动量/状态方程严格消元，
% 只对非平凡压缩根求解；物性越界只用于括区符号，不输出外推状态。
state=struct('available',false);
T1=fs.temperature; Un=fs.speed; R=fs.R; rho1=fs.density;
Mn=Un/sqrt(gammaPhysical*R*T1);
if Mn-1<1e-8, return; end
[~,h1]=airThermalProperties(T1); hn=h1+Un^2/2;
lo=1+min(1e-4,0.25*(Mn-1)); hi=1e4;
if densityEnergy(lo,T1,Un,R,hn)<=0 || densityEnergy(hi,T1,Un,R,hn)>=0, return; end
for k=1:72
    mid=0.5*(lo+hi);
    if densityEnergy(mid,T1,Un,R,hn)>0, lo=mid; else, hi=mid; end
end
r=0.5*(lo+hi); Un2=Un/r;
p2=fs.pressure+rho1*Un*(Un-Un2); rho2=rho1*r; T2=p2/(rho2*R);
if T2<80 || T2>2500 || r<=1 || p2<=fs.pressure || T2<=T1, return; end
state.available=true; state.postShockPressure=p2;
state.postShockTemperature=T2; state.postShockDensity=rho2; state.postShockSpeed=Un2;
end

function residual=densityEnergy(r,T1,Un,R,hn)
T=T1/r+(Un^2/R)*(r-1)/r^2;
if T>2500, residual=1; elseif T<80, residual=-1;
else, [~,h]=airThermalProperties(T); residual=h+(Un/r)^2/2-hn; end
end

function source=inputSource(cfg,key,fallback)
source=fallback;
if isfield(cfg,'externalFlow') && isfield(cfg.externalFlow,key)
    value=cfg.externalFlow.(key);
    if ~(ischar(value) && isrow(value)) && ~(isstring(value) && isscalar(value))
        error('RocketAeroHeat:ShockEntropyEdge:Source','externalFlow.%s 须为来源文字。',key);
    end
    source=char(value);
end
end

function meta=fieldMetadata(names,method,scope,pressureSource)
% 每个场量保留方法与输入/求解区别，不把给定静压标为求得压力。
meta=struct;
for j=1:numel(names)
    meta.(names{j})=struct('method',method,'source','同一激波穿越熵与完整总焓', ...
        'valueKind','computed','scope',scope,'valid',[],'statusCode',[]);
end
meta.pressure.method='输入物理静压原值'; meta.pressure.source=pressureSource;
meta.pressure.valueKind='input';
meta.entropyIncrease.method='激波后标准熵差减 R*ln(p2/p1)';
meta.totalPressure.method='激波后熵与完整速度总焓的同模型等熵总压';
meta.totalTemperature.method='上游完整速度总焓反算总温，绝热激波守恒';
end
