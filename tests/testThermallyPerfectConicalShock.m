function testThermallyPerfectConicalShock()
% 变比热切锥：核对弱支、锥面无穿透、激波跃迁及锥内总焓/熵。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();');
fs=calcFreestream(cfg.flight);

%% 1. 零角极限与自由来流的 NASA7 总状态一致
free=solveThermallyPerfectConicalShock(fs,0);
assert(free.available && free.attached);
assert(abs(free.pressure-fs.pressure)<1e-9);
assert(abs(free.speed-fs.speed)<1e-9);
assert(free.totalPressure>0 && free.totalTemperature>fs.temperature);
assert(strcmp(free.thermalPropertyMode,'thermallyPerfect'));

%% 2. 锥面解应落在 Mach 角和同角二维楔弱激波之间
cone=solveThermallyPerfectConicalShock(fs,10);
wedge=solveThermallyPerfectObliqueShock(fs,10);
assert(cone.available && cone.attached && wedge.available);
assert(cone.shockAngleDeg>asind(1/cone.upstreamMachThermal));
assert(cone.shockAngleDeg<wedge.shockAngleDeg);
assert(cone.pressure>fs.pressure && cone.pressure<wedge.pressure);
assert(cone.totalPressure<free.totalPressure && ...
    cone.totalPressure>cone.pressure);
assert(abs(cone.polarVelocityAtWall)<1e-3);
assert(abs(cone.pressure-cone.density*fs.R*cone.temperature)<1e-8);
assert(abs(cone.energyResidual)<1e-7);
assert(abs(cone.shockMassResidual)<1e-9);
assert(abs(cone.shockMomentumResidual)<1e-9);
assert(abs(cone.isentropicResidual)<1e-8);
assert(cone.entropyIncrease>0);
[~,h1,s1]=airThermalProperties(fs.temperature);
[~,hw,sw]=airThermalProperties(cone.temperature);
assert(abs(hw+cone.speed^2/2-h1-fs.speed^2/2)<1e-2);
assert(abs(sw-s1-fs.R*log(cone.pressure/fs.pressure) ...
    -cone.entropyIncrease)<1e-6);

%% 3. 变比热锥流可以接同模型的二维膨胀与压缩转角
upstream=struct('Mach',cone.speed/sqrt( ...
    fs.gamma*fs.R*cone.temperature),'gamma',fs.gamma,'R',fs.R, ...
    'pressure',cone.pressure,'temperature',cone.temperature, ...
    'density',cone.density,'speed',cone.speed, ...
    'totalPressure',cone.totalPressure, ...
    'totalTemperature',cone.totalTemperature);
expand=solveThermallyPerfectPrandtlMeyerExpansion(upstream,5);
assert(expand.available && expand.pressure<cone.pressure);
assert(abs(expand.totalPressure-cone.totalPressure)< ...
    1e-7*cone.totalPressure);
upstream.Mach=expand.speed/sqrt(fs.gamma*fs.R*expand.temperature);
upstream.pressure=expand.pressure;
upstream.temperature=expand.temperature;
upstream.density=expand.density;
upstream.speed=expand.speed;
upstream.totalPressure=expand.totalPressure;
upstream.totalTemperature=expand.totalTemperature;
compress=solveThermallyPerfectObliqueShock(upstream,3);
assert(compress.available && compress.attached);
assert(compress.totalPressure<expand.totalPressure);
assert(abs(compress.totalTemperature-cone.totalTemperature)<1e-5);

%% 4. 脱体角与总焓越界均须显式无解
detached=solveThermallyPerfectConicalShock(fs,80);
assert(~detached.attached && isnan(detached.pressure));
hot=cfg.flight; hot.Mach=10;
outside=solveThermallyPerfectConicalShock(calcFreestream(hot),10);
assert(~outside.available && isnan(outside.pressure));
fprintf('Thermal conical shock passed: beta %.5f deg, pwall/pinf %.5f.\n', ...
    cone.shockAngleDeg,cone.pressureRatio);
end
