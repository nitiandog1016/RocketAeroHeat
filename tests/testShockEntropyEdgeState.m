function testShockEntropyEdgeState()
%TESTSHOCKENTROPYEDGESTATE 独立前向穿越激波/外缘恢复的解析与守恒校核。
% 模块作用：验证已知激波角的法向速度、切向速度及完整总焓，检查
% 与壁面转角接口区别、NASA7 实际声速、Mach 波极限与真实域外 NaN。
% 制造的角度/物理静压仅用于内核校核，不是激波几何或质量匹配验证。
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
cfg=struct('heating',struct('airThermalModel','caloricallyPerfect'));
cfg.externalFlow=struct('physicalPressureSource','解析校核输入的静压', ...
    'shockCrossingAngleSource','解析校核输入的 Gamma；未求质量匹配');
fs=makeState(220,5000,4); before=fs;

%% 1. 独立完全气体解析：指定激波角，不调用现有激波函数作参考
angles=[asind(1/fs.Mach),25,40,60,90];
reference=caloricReference(fs,angles);
pressure=0.85*reference.p2;
edge=recoverShockEntropyEdgeState(fs,angles,pressure,cfg);
assert(isequaln(fs,before),'独立内核不得改写输入状态');
assert(edge.available && all(edge.algebraicallyValid));
cp=fs.gamma*fs.R/(fs.gamma-1);
T0=fs.temperature+fs.speed^2/(2*cp);
expectedT=T0*(pressure./reference.p0).^((fs.gamma-1)/fs.gamma);
expectedU=sqrt(2*cp*(T0-expectedT));
assert(max(abs(edge.shock.pressure./reference.p2-1))<1e-12);
assert(max(abs(edge.shock.temperature./reference.T2-1))<1e-12);
assert(max(abs(edge.shock.density./reference.rho2-1))<1e-12);
assert(max(abs(edge.shock.speed./reference.U2-1))<1e-12);
assert(max(abs(edge.totalPressure./reference.p0-1))<1e-12);
assert(max(abs(edge.temperature./expectedT-1))<1e-12);
assert(max(abs(edge.speed./expectedU-1))<1e-12);
assert(max(abs(edge.shock.tangentialSpeedBefore-edge.shock.tangentialSpeedAfter))==0);
assert(edge.shock.machWaveLimit(1) && edge.entropyIncrease(1)==0);
assert(all(diff(edge.totalPressure)<0));
assert(edge.shock.deflectionDeg(end)==0 && edge.shock.deflectionDeg(2)>0);
assert(max(abs(edge.shock.normalSpeedBefore-fs.speed*sind(angles)))<1e-12);
% 如果误把全速当法向速，25度和90度总压相同；本测试明确拒绝。
assert(edge.totalPressure(2)>1.5*edge.totalPressure(end));
assert(~edge.geometryVerified && ~edge.massMatchingSolved);
assert(strcmp(edge.fieldMetadata.pressure.valueKind,'input'));
assert(strcmp(edge.fieldMetadata.speed.valueKind,'computed'));
assert(strcmp(edge.fieldMetadata.pressure.source,cfg.externalFlow.physicalPressureSource));
assert(all(edge.fieldMetadata.temperature.valid));

%% 2. 从激波后静压恢复必须得到完整下游速度；驻点零速度是求得极限
atShock=recoverShockEntropyEdgeState(fs,angles,reference.p2,cfg);
assert(max(abs(atShock.temperature./reference.T2-1))<1e-12);
assert(max(abs(atShock.speed./reference.U2-1))<1e-12);
stag=recoverShockEntropyEdgeState(fs,angles,atShock.totalPressure,cfg);
assert(all(stag.algebraicallyValid) && all(stag.speed==0) && all(stag.Mach==0));
assert(all(stag.temperature==T0));
column=recoverShockEntropyEdgeState(fs,angles',reference.p2',cfg);
assert(isequal(size(column.pressure),size(angles')));
assert(max(abs(column.speed'-atShock.speed))<1e-10);
broadcast=recoverShockEntropyEdgeState(fs,40,[5000 6000;7000 8000],cfg);
assert(isequal(size(broadcast.pressure),[2 2]) && all(broadcast.algebraicallyValid(:)));

%% 3. 默认 NASA7：独立检查法向质量/动量、全速度能量和熵关系
thermal=struct; [cp1,h1,s1]=airThermalProperties(fs.temperature);
gPhysical=cp1/(cp1-fs.R); a1=sqrt(gPhysical*fs.R*fs.temperature);
mu=asind(a1/fs.speed); thermalAngles=[mu,mu+.05,30,60,90];
edgeT=recoverShockEntropyEdgeState(fs,thermalAngles,fs.pressure,thermal);
assert(edgeT.available && all(edgeT.algebraicallyValid));
assert(strcmp(edgeT.thermalPropertyMode,'thermallyPerfect'));
checkConservation(fs,edgeT,h1,s1);
assert(edgeT.shock.machWaveLimit(1) && edgeT.totalPressure(1)==edgeT.upstreamTotalPressure);
assert(abs(edgeT.temperature(1)-fs.temperature)<1e-9);
assert(abs(edgeT.speed(1)-fs.speed)<1e-8);
assert(edgeT.shock.entropyIncrease(2)>0 && ...
    edgeT.shock.entropyIncrease(2)<edgeT.shock.entropyIncrease(end));
assert(contains(edgeT.shock.normalKernelMethod{3},'复用 NASA7'));
recovered=recoverShockEntropyEdgeState(fs,thermalAngles,edgeT.shock.pressure,thermal);
assert(all(recovered.algebraicallyValid));
assert(max(abs(recovered.speed-edgeT.shock.speed))<1e-7);
assert(max(abs(recovered.temperature-edgeT.shock.temperature))<1e-8);
% NASA7 恢复核的行/列/矩阵输出契约；单点角度广播时密度不得发生二维扩张。
columnT=recoverShockEntropyEdgeState(fs,thermalAngles',fs.pressure,thermal);
assert(isequal(size(columnT.density),size(thermalAngles')));
assert(all(columnT.algebraicallyValid) && max(abs(columnT.speed'-edgeT.speed))<1e-8);
matrixPressure=[5000 6000;7000 8000];
matrixT=recoverShockEntropyEdgeState(fs,40,matrixPressure,thermal);
assert(all(matrixT.algebraicallyValid(:)));
assert(isequal(size(matrixT.temperature),[2 2]) && isequal(size(matrixT.speed),[2 2]) && ...
    isequal(size(matrixT.density),[2 2]));
assert(max(abs(matrixT.density(:).*fs.R.*matrixT.temperature(:)./matrixPressure(:)-1))<1e-12);

%% 4. 高静温弱激波：真实 gamma 低于旧核识别阈值，不能混用名义声速
hot=makeState(1200,9000,1.7);
[cpHot,hHot,sHot]=airThermalProperties(hot.temperature);
gHot=cpHot/(cpHot-hot.R); aHot=sqrt(gHot*hot.R*hot.temperature);
assert(gHot<1.35);
hotMu=asind(aHot/hot.speed);
hotAngles=[hotMu,asind(1.02*aHot/hot.speed),asind(1.1*aHot/hot.speed),90];
hotEdge=recoverShockEntropyEdgeState(hot,hotAngles,hot.pressure,thermal);
assert(all(hotEdge.algebraicallyValid));
assert(abs(hotEdge.machAngleDeg-hotMu)<1e-12);
assert(hotEdge.shock.machWaveLimit(1) && hotEdge.entropyIncrease(1)==0);
assert(contains(hotEdge.shock.normalKernelMethod{2},'密度比'));
assert(hotEdge.shock.normalMachBefore(2)>1 && hotEdge.shock.normalMachBefore(2)<1.03);
assert(hotEdge.shock.pressure(2)>hot.pressure && hotEdge.shock.temperature(2)>hot.temperature);
assert(hotEdge.totalPressure(2)<hotEdge.upstreamTotalPressure);
assert(hotEdge.shock.density(2)/hot.density<1.05);
checkConservation(hot,hotEdge,hHot,sHot);
% 独立温度根法复核弱激波，不复用被测函数的密度比二分。
Un1=hot.speed*sind(hotAngles(2)); Ut1=hot.speed*cosd(hotAngles(2));
normalEnthalpy=hHot+Un1^2/2;
rootT=fzero(@(T)temperatureEnergy(T,hot,Un1,normalEnthalpy), ...
    [hot.temperature*1.0001,hot.temperature*1.1]);
[~,independentH]=airThermalProperties(rootT);
independentUn2=sqrt(2*(normalEnthalpy-independentH));
assert(abs(hotEdge.shock.temperature(2)-rootT)<1e-6);
assert(abs(hotEdge.shock.speed(2)-hypot(independentUn2,Ut1))<1e-6);
fprintf('High-temperature weak shock: gamma %.6f, Mn %.6f, T2 %.6f K, p0/p01 %.9f.\n', ...
    gHot,hotEdge.shock.normalMachBefore(2),hotEdge.shock.temperature(2), ...
    hotEdge.totalPressure(2)/hotEdge.upstreamTotalPressure);

%% 5. 逐点域外/物性越界/数值失败 NaN；结构与状态错误须报错
badAngles=[-1,mu-.1,91,NaN,Inf];
bad=recoverShockEntropyEdgeState(fs,badAngles,fs.pressure,thermal);
assert(~bad.available && all(bad.statusCode==1)); assertNaNEdge(bad);
badPressure=recoverShockEntropyEdgeState(fs,40,[-1,0,NaN,Inf,1e99],thermal);
assert(~badPressure.available && all(badPressure.statusCode==2)); assertNaNEdge(badPressure);
tooCold=recoverShockEntropyEdgeState(fs,40,1e-8,thermal);
assert(~tooCold.available && tooCold.statusCode==3); assertNaNEdge(tooCold);
tooHot=recoverShockEntropyEdgeState(makeState(220,5000,10),40,5000,thermal);
assert(~tooHot.available && tooHot.statusCode==3); assertNaNEdge(tooHot);
% 上游物性早退与外缘恢复越界都须同步逐字段状态，不能遗留“角度不适用”。
upstreamCold=recoverShockEntropyEdgeState(makeState(70,5000,4),[30 40],5000,thermal);
assert(~upstreamCold.available && all(upstreamCold.statusCode==3)); assertNaNEdge(upstreamCold);
subsonic=recoverShockEntropyEdgeState(makeState(220,5000,.5),40,5000,cfg);
assert(~subsonic.available); assertNaNEdge(subsonic);
% 数值分辨率以下的非零激波不可假装为零强度 Mach 波。
weakAngle=asind((1+1e-9)*aHot/hot.speed);
weakFail=recoverShockEntropyEdgeState(hot,weakAngle,hot.pressure,thermal);
assert(~weakFail.available && weakFail.statusCode==4); assertNaNEdge(weakFail);
expectError(@()recoverShockEntropyEdgeState(fs,[30 40],[5000;6000],cfg), ...
    'RocketAeroHeat:ShockEntropyEdge:Size');
broken=fs; broken.speed=fs.speed*1.1;
expectError(@()recoverShockEntropyEdgeState(broken,40,5000,cfg), ...
    'RocketAeroHeat:ShockEntropyEdge:State');
broken=rmfield(fs,'density');
expectError(@()recoverShockEntropyEdgeState(broken,40,5000,cfg), ...
    'RocketAeroHeat:ShockEntropyEdge:Input');
expectError(@()recoverShockEntropyEdgeState(fs,40+1i,5000,cfg),'');
badCfg=cfg;badCfg.externalFlow.physicalPressureSource=42;
expectError(@()recoverShockEntropyEdgeState(fs,40,5000,badCfg), ...
    'RocketAeroHeat:ShockEntropyEdge:Source');
fprintf('Shock entropy edge kernel passed: independent caloric formulas; full velocity/enthalpy; NASA7 and hot weak-shock conservation; Mach wave; NaN/error rejection.\n');
fprintf('Standalone thermodynamic checks only; shock geometry, mass matching and experimental accuracy remain unvalidated.\n');
end

function fs=makeState(T,p,M)
% 只构造与名义 gamma/状态方程一致的上游；NASA7 实际 Mach 由内核另算。
gamma=1.4; R=287.05;
fs=struct('gamma',gamma,'R',R,'Mach',M,'temperature',T,'pressure',p, ...
    'density',p/(R*T),'speed',M*sqrt(gamma*R*T));
end

function ref=caloricReference(fs,beta)
% NACA 1135 法向分量闭式公式；保持切向速度，形成完整下游总状态。
g=fs.gamma; Mn=fs.Mach*sind(beta); r=(g+1)*Mn.^2./((g-1)*Mn.^2+2);
ratio=1+2*g/(g+1)*(Mn.^2-1);
ref.p2=fs.pressure*ratio; ref.rho2=fs.density*r;
ref.T2=fs.temperature*ratio./r;
Un2=fs.speed*sind(beta)./r; Ut=fs.speed*cosd(beta);
ref.U2=hypot(Un2,Ut);
cp=g*fs.R/(g-1);T0=fs.temperature+fs.speed^2/(2*cp);
ref.p0=ref.p2.*(T0./ref.T2).^(g/(g-1));
end

function checkConservation(fs,edge,h1,s1)
% 直接从输出原始量算守恒，不依赖内核自报 residual。
s=edge.shock; [~,h2,s2]=airThermalProperties(s.temperature);
mass0=fs.density*s.normalSpeedBefore;
momentum0=fs.pressure+fs.density*s.normalSpeedBefore.^2;
assert(max(abs((s.density.*s.normalSpeedAfter-mass0)./mass0))<1e-10);
assert(max(abs((s.pressure+s.density.*s.normalSpeedAfter.^2-momentum0)./momentum0))<1e-10);
assert(max(abs(h2+s.speed.^2/2-h1-fs.speed^2/2))<2e-6);
entropy=(s2-s1)-fs.R*log(s.pressure/fs.pressure);
assert(max(abs(entropy-s.entropyIncrease))<1e-6);
assert(max(abs(s.entropyIncrease-fs.R*log(edge.upstreamTotalPressure./s.totalPressure)))<1e-6);
assert(max(abs(s.tangentialSpeedBefore-s.tangentialSpeedAfter))==0);
assert(max(abs(edge.entropyResidual))<1e-6 && max(abs(edge.energyResidual))<1e-8);
end

function residual=temperatureEnergy(T,fs,Un1,normalEnthalpy)
% 以温度作为独立未知，能量式先给 Un2，再用法向动量和质量消元。
[~,h]=airThermalProperties(T); Un2=sqrt(2*(normalEnthalpy-h));
mass=fs.density*Un1;rho2=mass/Un2;p2=rho2*fs.R*T;
residual=p2+rho2*Un2^2-fs.pressure-fs.density*Un1^2;
end

function assertNaNEdge(edge)
fields={'pressure','temperature','density','speed','Mach','speedOfSound', ...
    'totalPressure','totalTemperature','entropyIncrease','energyResidual','entropyResidual'};
for j=1:numel(fields)
    key=fields{j}; assert(all(isnan(edge.(key)(:))));
    assert(isequal(edge.fieldMetadata.(key).valid,edge.algebraicallyValid));
    assert(isequal(edge.fieldMetadata.(key).statusCode,edge.statusCode));
end
end

function expectError(action,id)
failed=false;
try, action(); catch errorValue
    failed=isempty(id) || strcmp(errorValue.identifier,id);
end
assert(failed,'预期输入错误没有报错或标识不符');
end
