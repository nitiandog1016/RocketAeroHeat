function testThermallyPerfectNormalShock()
% 默认 Mach 5 鼻端：变比热正激波应守恒、增熵且与定比热有可解释差异。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();');
fs=calcFreestream(cfg.flight);

%% 1. 独立单点跃迁的质量、动量、总焓和总压损失
shock=solveThermallyPerfectNormalShock(fs);
old=solveNormalShock(fs);
assert(shock.available);
assert(strcmp(shock.thermalPropertyMode,'thermallyPerfect'));
assert(shock.postShockMach>0 && shock.postShockMach<1);
assert(shock.postShockPressure>fs.pressure);
assert(shock.postShockTemperature>fs.temperature);
assert(shock.totalPressureRatio>0 && shock.totalPressureRatio<1);
assert(shock.entropyIncrease>0);
assert(max(abs([shock.massResidual,shock.momentumResidual, ...
    shock.energyResidual]))<1e-10);
assert(shock.stagnationTemperature<old.stagnationTemperature);
assert(shock.postShockTemperature<old.postShockTemperature);
assert(abs(shock.stagnationTemperature- ...
    airRecoveryTemperature(fs.temperature,fs.speed,1))<1e-8);
[~,hBefore,sBefore]=airThermalProperties(fs.temperature);
[~,hAfter,sAfter]=airThermalProperties(shock.postShockTemperature);
assert(abs(hBefore+fs.speed^2/2-hAfter- ...
    shock.postShockSpeed^2/2)<1e-6);
assert(abs(shock.entropyIncrease-(sAfter-sBefore)+ ...
    fs.R*log(shock.pressureRatio))<1e-7);
% 弱激波的非平凡压缩根可能位于速度二次方程的快根上，不能漏解。
[cpWeak,~]=airThermalProperties(fs.temperature);
gammaWeak=cpWeak/(cpWeak-fs.R);
weak=fs;
weak.speed=1.02*sqrt(gammaWeak*fs.R*fs.temperature);
weak.Mach=weak.speed/sqrt(fs.gamma*fs.R*fs.temperature);
weakJump=solveThermallyPerfectNormalShock(weak);
assert(weakJump.available && weakJump.postShockMach<1);
assert(weakJump.postShockTemperature>weak.temperature);
assert(weakJump.densityRatio>1 && weakJump.densityRatio<1.1);
assert(weakJump.totalPressureRatio<1);
assert(max(abs([weakJump.massResidual,weakJump.momentumResidual, ...
    weakJump.energyResidual]))<1e-9);

%% 2. 适用域与旧版选项；不能在高温时静默外推
hot=cfg.flight; hot.Mach=10;
tooHot=solveThermallyPerfectNormalShock(calcFreestream(hot));
assert(~tooHot.available && contains(tooHot.status,'上限'));
cfg.heating.airThermalModel='caloricallyPerfect';
assert(strcmp(selectAirThermalModel(fs,cfg),'caloricallyPerfect'));
fprintf('Thermal normal shock passed: T2 %.1f K vs caloric %.1f K, p02/p01 %.4f.\n', ...
    shock.postShockTemperature,old.postShockTemperature, ...
    shock.totalPressureRatio);
end
