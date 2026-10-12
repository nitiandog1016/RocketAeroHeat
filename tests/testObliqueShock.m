function testObliqueShock()
% 解析关系和边界工况：零转折、附体弱解、附体极限及脱体拒绝。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); fs=calcFreestream(cfg.flight);
zero=solveObliqueShock(fs,0);
assert(zero.attached && zero.Mach==fs.Mach);
assert(zero.pressure==fs.pressure && zero.totalPressure==fs.totalPressure);
assert(zero.entropyIncrease==0);
assert(abs(zero.shockAngleDeg-asind(1/fs.Mach))<1e-11);

small=solveObliqueShock(fs,10);
assert(small.attached && strcmp(small.branch,'weak'));
assert(small.shockAngleDeg>zero.shockAngleDeg);
assert(small.shockAngleDeg<small.peakShockAngleDeg);
assert(small.pressure>fs.pressure && small.temperature>fs.temperature);
assert(small.Mach<fs.Mach && small.totalPressure<fs.totalPressure);
assert(small.entropyIncrease>0);
assert(abs(small.totalTemperature-fs.totalTemperature)<1e-10*fs.totalTemperature);
% 质量守恒：法向密度与速度比互为倒数。
assert(abs(small.normalMachBefore*sqrt(fs.gamma*fs.R*fs.temperature)/ ...
    (small.normalMachAfter*sqrt(fs.gamma*fs.R*small.temperature)) ...
    -small.densityRatio)<1e-10);
% 直接代回 theta-beta-M 方程，不能只验证代码自行输出的压力比。
beta=deg2rad(small.shockAngleDeg); M=fs.Mach; g=fs.gamma;
theta=atan(2*cot(beta)*(M^2*sin(beta)^2-1)/ ...
    (M^2*(g+cos(2*beta))+2));
assert(abs(rad2deg(theta)-10)<1e-9);

limit=small.maximumAttachedDeflectionDeg;
atLimit=solveObliqueShock(fs,limit);
assert(atLimit.attached && abs(atLimit.shockAngleDeg- ...
    atLimit.peakShockAngleDeg)<1e-6);
detached=solveObliqueShock(fs,min(89,limit+1));
assert(~detached.attached && isnan(detached.pressure));
assert(contains(detached.status,'脱体'));
rejected=false;
try, solveObliqueShock(fs,-1); catch, rejected=true; end
assert(rejected);
fprintf('Oblique shock passed: M=%.1f, theta=10 deg -> beta=%.5f deg, p2/p1=%.5f; thetaMax=%.5f deg.\n', ...
    M,small.shockAngleDeg,small.pressureRatio,limit);
end
