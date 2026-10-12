function testWallFeedback()
% 壁温反馈应改变 h，且集中热容与一维导热均满足离散能量守恒。
addpath(fileparts(fileparts(mfilename('fullpath'))));
m.F=[1 2 3];
m.heating.heatTransferCoefficient=100;
m.heating.recoveryTemperature=1000;
m.heating.eligibleFace=true;
m.heating.laminarFace=true;
m.heating.turbulentFace=false;
m.heating.Prandtl=0.72;
m.heating.wallTemperature=300;
m.boundaryLayerEdge.temperature=300;
m.boundaryLayerEdge.pressure=101325;
m.boundaryLayerEdge.freestream.R=287;
cfg.flight.Twall=300;
material=struct('density',1000,'specificHeat',1000, ...
    'thickness',0.01,'thermalConductivity',10, ...
    'emissivity',0,'radiationTemperature',300);
t=[0 10 20];
warm=calcWallThermalResponse(m,t,material,cfg);
w=warm.wallThermal;
assert(w.wallTemperatureFeedback && w.lumpedFace);
assert(abs(w.heatTransferCoefficientHistory(1,1)-100)<1e-10);
assert(w.heatTransferCoefficientHistory(1,2)>100);
C=material.density*material.specificHeat*material.thickness;
for j=2:numel(t)
    residual=C/(t(j)-t(j-1))*(w.temperature(1,j)-w.temperature(1,j-1)) ...
        -w.convectiveHeatFlux(1,j);
    assert(abs(residual)<1e-8*max(1,w.convectiveHeatFlux(1,j)));
end
borderline=material;
borderline.thermalConductivity=100*borderline.thickness/(2*0.099);
conservative=calcWallThermalResponse(m,t,borderline,cfg);
assert(conservative.wallThermal.Biot<0.1);
assert(conservative.wallThermal.maximumEstimatedBiot>0.1);
assert(conservative.wallThermal.oneDimensionalFace);
material.thermalConductivity=0.1;
thick=calcWallThermalResponse(m,t,material,cfg);
w=thick.wallThermal;
assert(w.oneDimensionalFace && w.wallTemperatureFeedback);
profile=w.oneDimensional.finalDepthTemperature;
dx=material.thickness/(size(profile,2)-1);
stored=material.density*material.specificHeat*dx*( ...
    sum(profile(2:end-1))+0.5*(profile(1)+profile(end)) ...
    -300*(size(profile,2)-1));
inputEnergy=sum(diff(t).*w.convectiveHeatFlux(1,2:end));
assert(abs(stored-inputEnergy)<1e-8*max(1,inputEnergy));
stale=cfg; stale.flight.Twall=301;
failed=false;
try, calcWallThermalResponse(m,t,material,stale); catch, failed=true; end
assert(failed);
fprintf('Wall feedback passed: h(0)=%.3f, h(10 s)=%.3f W/(m2 K).\n', ...
    warm.wallThermal.heatTransferCoefficientHistory(1,1), ...
    warm.wallThermal.heatTransferCoefficientHistory(1,2));
end
