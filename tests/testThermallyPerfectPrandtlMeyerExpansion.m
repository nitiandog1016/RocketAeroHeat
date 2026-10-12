function testThermallyPerfectPrandtlMeyerExpansion()
% 干空气二维等熵膨胀：独立角度积分、总焓/熵守恒与低温边界。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();');
fs=calcFreestream(cfg.flight);

%% 1. 零转角保持静态，10° 背风转角单调膨胀
zero=solveThermallyPerfectPrandtlMeyerExpansion(fs,0);
turn=solveThermallyPerfectPrandtlMeyerExpansion(fs,10);
caloric=solvePrandtlMeyerExpansion(fs,10);
assert(zero.available && zero.pressure==fs.pressure);
assert(zero.temperature==fs.temperature && zero.speed==fs.speed);
assert(turn.available && turn.Mach>zero.Mach);
assert(turn.pressure<fs.pressure && turn.temperature<fs.temperature);
assert(turn.speed>fs.speed && turn.totalPressure==zero.totalPressure);
assert(turn.totalTemperature==zero.totalTemperature);
assert(abs(turn.pressureRatio/caloric.pressureRatio-1)<0.02);
assert(abs(turn.entropyResidual)<1e-8);
assert(abs(turn.energyResidual)<1e-10);

%% 2. 用独立自适应积分校核给定转角，不能只检验输出趋势
[~,h1,s1]=airThermalProperties(fs.temperature);
[~,h2,s2]=airThermalProperties(turn.temperature);
h0=h1+0.5*fs.speed^2;
angle=integral(@(T)turnIntegrand(T,h0,fs.R), ...
    turn.temperature,fs.temperature,'RelTol',1e-10,'AbsTol',1e-11);
assert(abs(angle-deg2rad(10))<2e-6);
assert(abs(h2+0.5*turn.speed^2-h0)<1e-6);
assert(abs((s2-s1)-fs.R*log(turn.pressure/fs.pressure))<1e-7);
assert(abs(turn.pressure-turn.density*fs.R*turn.temperature)<1e-8);

%% 3. 跨过 150 K 接缝时，膨胀角与外缘状态反算仍连续
cold=solveThermallyPerfectPrandtlMeyerExpansion(fs,15);
assert(cold.available && cold.temperature>=80 && cold.temperature<150);
coldAngle=integral(@(T)turnIntegrand(T,h0,fs.R), ...
    cold.temperature,fs.temperature,'RelTol',1e-10,'AbsTol',1e-11);
assert(abs(coldAngle-deg2rad(15))<2e-6);
recovered=recoverThermallyPerfectEdgeState( ...
    cold.pressure,cold.totalPressure,cold.totalTemperature,fs.R);
assert(recovered.valid && abs(recovered.temperature-cold.temperature)<1e-8);
assert(abs(recovered.speed-cold.speed)<1e-6);
assert(abs(recovered.entropyResidual)<1e-7);

%% 4. 超出 80 K 下限或总温上限时明确无解
assert(turn.maximumAdditionalTurnDeg>10);
outside=solveThermallyPerfectPrandtlMeyerExpansion( ...
    fs,min(89,turn.maximumAdditionalTurnDeg+1));
assert(~outside.available && isnan(outside.temperature));
assert(contains(outside.status,'80 K'));
hot=cfg.flight; hot.Mach=10;
hotState=solveThermallyPerfectPrandtlMeyerExpansion( ...
    calcFreestream(hot),5);
assert(~hotState.available && contains(hotState.status,'热物性'));
fprintf('Thermal Prandtl-Meyer passed: turn 10 deg, M2 %.5f, T2 %.2f K, p2/p1 %.5f.\n', ...
    turn.Mach,turn.temperature,turn.pressureRatio);
end

function value=turnIntegrand(T,h0,R)
% 与主程序的离散累积独立：MATLAB 自适应积分验证特征关系。
[cp,h]=airThermalProperties(T);
U=sqrt(2*(h0-h));
sound=sqrt(cp./(cp-R).*R.*T);
value=sqrt((U./sound).^2-1).*cp./U.^2;
end
