function testStagnationHeating()
% 验证半径依赖、公式值、缺失输入和超出适用攻角。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); fs=calcFreestream(cfg.flight);
m.stagnation=struct('representativeFaceID',1, ...
    'representativePosition',[0 0 0], ...
    'hasPhysicalRadius',false,'physicalRadius',NaN, ...
    'radiusSource','未提供','numberComponents',1);
m=calcStagnationHeating(m,fs,cfg);
assert(~m.stagnationHeating.available && isnan(m.stagnationHeating.heatFlux));
m.stagnation.hasPhysicalRadius=true;
m.stagnation.physicalRadius=0.01;
m.stagnation.radiusSource='测试给定物理半径';
m=calcStagnationHeating(m,fs,cfg);
expected=1.83e-4*sqrt(fs.density/0.01)*fs.speed^3;
assert(m.stagnationHeating.available);
assert(abs(m.stagnationHeating.heatFlux-expected)<1e-10*expected);
assert(abs(m.stagnationHeating.referenceHeatFlux-expected)<1e-10*expected);
hotWall=cfg; hotWall.flight.Twall=500;
corrected=calcStagnationHeating(m,fs,hotWall);
thermalTotal=airRecoveryTemperature(fs.temperature,fs.speed,1);
[~,hTotal]=airThermalProperties(thermalTotal);
[~,hHot]=airThermalProperties(500);
[~,hReference]=airThermalProperties(300);
expectedHot=expected*(hTotal-hHot)/(hTotal-hReference);
assert(abs(corrected.stagnationHeating.heatFlux-expectedHot)<1e-10*expected);
assert(corrected.stagnationHeating.heatFlux<expected);
oldMode=hotWall;
oldMode.heating.airThermalModel='caloricallyPerfect';
old=calcStagnationHeating(m,fs,oldMode);
oldExpected=expected*(fs.totalTemperature-500)/ ...
    (fs.totalTemperature-300);
assert(abs(old.stagnationHeating.heatFlux-oldExpected)<1e-10*expected);
m.stagnation.physicalRadius=0.04;
m=calcStagnationHeating(m,fs,cfg);
assert(abs(m.stagnationHeating.heatFlux-expected/2)<1e-10*expected);
high=fs; high.alpha=10;
m=calcStagnationHeating(m,high,cfg);
assert(~m.stagnationHeating.available && isnan(m.stagnationHeating.heatFlux));
fprintf('Stagnation heating passed: Rn=0.01 m -> q=%.3f W/cm2.\n',expected/1e4);
end
