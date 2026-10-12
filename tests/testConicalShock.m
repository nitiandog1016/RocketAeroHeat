function testConicalShock()
% 零角极限、附体锥形激波、锥面边界条件与二维楔解对照。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); fs=calcFreestream(cfg.flight);
zero=solveConicalShock(fs,0);
assert(zero.attached && zero.Mach==fs.Mach && zero.pressure==fs.pressure);
cone=solveConicalShock(fs,10);
wedge=solveObliqueShock(fs,10);
assert(cone.attached && contains(cone.method,'Taylor-Maccoll'));
assert(cone.shockAngleDeg>asind(1/fs.Mach));
assert(cone.shockAngleDeg<wedge.shockAngleDeg);
assert(cone.pressure>fs.pressure && cone.pressure<wedge.pressure);
assert(cone.totalPressure<fs.totalPressure && cone.totalPressure>0);
assert(abs(cone.polarVelocityAtWall)<1e-3);
assert(abs(cone.totalTemperature-fs.totalTemperature)<1e-10*fs.totalTemperature);
assert(abs(cone.pressure-cone.density*fs.R*cone.temperature)<1e-8*cone.pressure);
large=solveConicalShock(fs,80);
assert(~large.attached && isnan(large.pressure));
fprintf('Conical shock passed: M=%.1f, cone=10 deg -> beta=%.5f deg, pwall/pinf=%.5f.\n', ...
    fs.Mach,cone.shockAngleDeg,cone.pressureRatio);
end
