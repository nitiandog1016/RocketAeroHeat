function testPrandtlMeyerExpansion()
% 零转折、膨胀降压、总状态守恒和适用极限。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); fs=calcFreestream(cfg.flight);
zero=solvePrandtlMeyerExpansion(fs,0);
assert(zero.available && zero.Mach==fs.Mach);
assert(zero.pressure==fs.pressure && zero.totalPressure==fs.totalPressure);
turn=solvePrandtlMeyerExpansion(fs,10);
assert(turn.available && turn.Mach>fs.Mach);
assert(turn.pressure<fs.pressure && turn.temperature<fs.temperature);
assert(abs(turn.totalPressure-fs.totalPressure)<1e-12*fs.totalPressure);
assert(abs(turn.totalTemperature-fs.totalTemperature)<1e-12*fs.totalTemperature);
assert(abs(turn.pressure-turn.density*fs.R*turn.temperature)<1e-10*turn.pressure);
fprintf('Prandtl-Meyer passed: M=%.1f, turn=10 deg -> M=%.5f, p2/p1=%.5f.\n', ...
    fs.Mach,turn.Mach,turn.pressureRatio);
end
