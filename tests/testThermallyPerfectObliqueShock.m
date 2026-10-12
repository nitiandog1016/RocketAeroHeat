function testThermallyPerfectObliqueShock()
% 二维干空气附体弱斜激波：守恒、熵增、附体极限及等熵反算回归。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();');
fs=calcFreestream(cfg.flight);

%% 1. 非零转角的几何关系与独立守恒校核
shock=solveThermallyPerfectObliqueShock(fs,10);
assert(shock.available && shock.attached);
assert(strcmp(shock.thermalPropertyMode,'thermallyPerfect'));
assert(shock.shockAngleDeg>shock.machAngleDeg);
assert(shock.shockAngleDeg<shock.peakShockAngleDeg);
assert(shock.pressure>fs.pressure && shock.temperature>fs.temperature);
assert(shock.totalPressureRatio>0 && shock.totalPressureRatio<1);
assert(shock.entropyIncrease>0);
assert(max(abs([shock.massResidual,shock.momentumResidual, ...
    shock.energyResidual]))<1e-10);
beta=deg2rad(shock.shockAngleDeg);
Un1=fs.speed*sin(beta);
Un2=shock.speed*sin(beta-deg2rad(10));
Ut1=fs.speed*cos(beta);
Ut2=shock.speed*cos(beta-deg2rad(10));
assert(abs(Ut1-Ut2)<1e-7);
assert(abs(fs.density*Un1-shock.density*Un2)< ...
    1e-9*fs.density*Un1);
assert(abs(fs.pressure+fs.density*Un1^2- ...
    shock.pressure-shock.density*Un2^2)< ...
    1e-9*(fs.pressure+fs.density*Un1^2));
[~,h1,s1]=airThermalProperties(fs.temperature);
[~,h2,s2]=airThermalProperties(shock.temperature);
assert(abs(h1+fs.speed^2/2-h2-shock.speed^2/2)<1e-6);
assert(abs(shock.entropyIncrease-(s2-s1)+ ...
    fs.R*log(shock.pressure/fs.pressure))<1e-7);

%% 2. 从总状态独立恢复静态；无效域不得制造边界层外缘值
edge=recoverThermallyPerfectEdgeState( ...
    [shock.pressure;shock.pressure*2;NaN], ...
    [shock.totalPressure;shock.pressure;shock.totalPressure], ...
    shock.totalTemperature,fs.R);
assert(edge.valid(1));
assert(abs(edge.temperature(1)-shock.temperature)<1e-8);
assert(abs(edge.speed(1)-shock.speed)<1e-7);
assert(abs(edge.Mach(1)-shock.Mach)<1e-8);
assert(abs(edge.entropyResidual(1))<1e-7);
assert(~edge.valid(2) && ~edge.valid(3));

%% 3. 零转角、附体极限、脱体拒绝与热物性上限
zero=solveThermallyPerfectObliqueShock(fs,0);
assert(zero.available && zero.attached && zero.pressure==fs.pressure);
assert(zero.totalPressureRatio==1 && zero.entropyIncrease==0);
assert(shock.maximumAttachedDeflectionDeg>10);
% 从 Mach 角到附体极限，弱激波角度节点都应有连续热力解。
for turn=0.5:0.5:floor(shock.maximumAttachedDeflectionDeg*2)/2
    node=solveThermallyPerfectObliqueShock(fs,turn);
    assert(node.available && node.attached);
    assert(node.totalPressureRatio>0 && node.totalPressureRatio<1);
end
atLimit=solveThermallyPerfectObliqueShock( ...
    fs,shock.maximumAttachedDeflectionDeg);
assert(atLimit.available && atLimit.attached);
assert(abs(atLimit.shockAngleDeg-atLimit.peakShockAngleDeg)<1e-5);
detached=solveThermallyPerfectObliqueShock( ...
    fs,min(89,shock.maximumAttachedDeflectionDeg+1));
assert(~detached.available && ~detached.attached);
assert(isnan(detached.pressure) && contains(detached.status,'脱体'));
hot=cfg.flight; hot.Mach=10;
outside=solveThermallyPerfectObliqueShock(calcFreestream(hot),10);
assert(~outside.available && contains(outside.status,'热物性'));
fprintf('Thermal oblique shock passed: beta %.4f deg, p2/p1 %.4f, T2 %.1f K, max turn %.3f deg.\n', ...
    shock.shockAngleDeg,shock.pressureRatio,shock.temperature, ...
    shock.maximumAttachedDeflectionDeg);
end
