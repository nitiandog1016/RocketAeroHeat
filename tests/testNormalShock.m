function testNormalShock()
% 用 Rankine-Hugoniot 解析比值与三项守恒检查局部驻点法向激波。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();');
fs=calcFreestream(cfg.flight);
s=solveNormalShock(fs);
g=fs.gamma; M=fs.Mach;
expectedP=1+2*g/(g+1)*(M^2-1);
expectedRho=(g+1)*M^2/((g-1)*M^2+2);
assert(s.available && abs(s.pressureRatio-expectedP)<1e-12);
assert(abs(s.densityRatio-expectedRho)<1e-12);
assert(s.postShockMach<1 && s.postShockMach>0);
assert(s.totalPressureRatio<1 && s.totalPressureRatio>0);
assert(s.entropyIncrease>0);
assert(abs(s.postShockTotalTemperature/fs.totalTemperature-1)<1e-12);
assert(max(abs([s.massResidual s.momentumResidual s.energyResidual]))<1e-12);
assert(isnan(s.shockStandOffDistance));
bad=fs; bad.Mach=1;
failed=false; try, solveNormalShock(bad); catch, failed=true; end
assert(failed);
bad=fs; bad.density=2*fs.density;
failed=false; try, solveNormalShock(bad); catch, failed=true; end
assert(failed);
fprintf('Normal shock passed: p2/p1=%.3f, M2=%.4f, p02/p01=%.4f.\n', ...
    s.pressureRatio,s.postShockMach,s.totalPressureRatio);
end
