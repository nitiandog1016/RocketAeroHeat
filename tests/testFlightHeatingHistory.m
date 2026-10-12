function testFlightHeatingHistory()
% 默认箭体真实网格：同一网格逐站重算内置外流，并驱动变工况壁温。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);
[raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg);
mesh=analyzeLocalGeometry(mesh,cfg);
flight=struct('time',[0 0.2 0.4], ...
    'altitude',[30000 20000 30000], ...
    'Mach',[5 5 5], ...
    'alpha',[0 0 0],'beta',[0 0 0]);
history=calcFlightHeatingHistory(mesh,flight,cfg);
assert(isequal(size(history.heatTransferCoefficient),[size(mesh.F,1) 3]));
assert(isequal(size(history.faceReasonCode),[size(mesh.F,1) 3]));
assert(history.summary.inputStations==3 && ...
    history.summary.evaluationStations==3);
assert(isequal(history.faceReasonCode==0,history.eligibleFace));
assert(size(history.summary.reasonAreaByStation,1)==numel(history.reasonNames));
assert(all(history.freestreamKnudsenByStation<0.01));
assert(all(abs(sum(history.summary.reasonAreaByStation,1)- ...
    sum(mesh.localGeometry.area))<1e-9));
assert(all(history.summary.calculatedAreaFractionByStation>0.9));
assert(all(history.axisymmetricMarchAppliedByStation));
assert(all(history.summary.axialMarchFacesByStation>10000));
assert(history.summary.validAcrossHistory>1000);
assert(all(history.stagnation.available));
assert(history.stagnation.heatTransferCoefficient(2)> ...
    history.stagnation.heatTransferCoefficient(1));
assert(history.summary.noseFacesByStation(1)>10);
% 长间隔自动插入中间计算站；完全相同的稳态工况只求解一次。
constantFlight=struct('time',[0 2], ...
    'altitude',[30000 30000],'Mach',[5 5], ...
    'alpha',[0 0],'beta',[0 0]);
constantHistory=calcFlightHeatingHistory(mesh,constantFlight,cfg);
assert(isequal(constantHistory.time,[0 1 2]));
assert(isequal(constantHistory.reusedSteadyState,[false true true]));
assert(constantHistory.summary.reusedSteadyStates==2);
assert(isequaln(constantHistory.heatTransferCoefficient(:,1), ...
    constantHistory.heatTransferCoefficient(:,3)));
varyingFlight=constantFlight;
varyingFlight.alpha=[0 10];
varyingHistory=calcFlightHeatingHistory(mesh,varyingFlight,cfg);
assert(isequal(varyingHistory.time,[0 1 2]));
assert(abs(varyingHistory.flight.alpha(2)-5)<1e-12);
assert(abs(varyingHistory.freestream{2}.alpha-5)<1e-12);
assert(~any(varyingHistory.reusedSteadyState));
% 整段高空稀薄流不应伪造连续介质壁面边界条件。
rareFlight=struct('time',[0 0.2], ...
    'altitude',[80000 80000],'Mach',[5 5], ...
    'alpha',[0 0],'beta',[0 0]);
rareHistory=calcFlightHeatingHistory(mesh,rareFlight,cfg);
assert(~any(rareHistory.eligibleFace,'all'));
assert(all(rareHistory.freestreamKnudsenByStation>0.01));
assert(any(rareHistory.faceReasonCode==10,'all'));
assert(~any(rareHistory.stagnation.available));
material=struct('density',2700,'specificHeat',900,'thickness',0.003, ...
    'thermalConductivity',180,'emissivity',0.8, ...
    'radiationTemperature',cfg.flight.Twall);
updated=calcWallThermalHistory(mesh,history,material,cfg);
updated=calcAerothermalExposure(updated);
wall=updated.wallThermal;
exposure=updated.aerothermalExposure;
assert(wall.summary.calculatedFaces>1000);
assert(exposure.summary.calculatedFaces==wall.summary.calculatedFaces);
assert(exposure.summary.resolvedAreaFraction>0.9);
assert(exposure.summary.totalNetEnergyResolved>0);
assert(exposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(abs(exposure.summary.totalConvectiveEnergyResolved- ...
    exposure.summary.totalRadiativeEnergyResolved- ...
    exposure.summary.totalNetEnergyResolved)<1e-6);
assert(all(isnan(exposure.cumulativeNetHeat(~wall.validFace,:)),'all'));
assert(all(isnan(exposure.peakConvectiveHeatFlux(~wall.validFace))));
assert(max(abs(exposure.peakConvectiveHeatFlux(wall.validFace)- ...
    max(wall.convectiveHeatFlux(wall.validFace,:),[],2)))<1e-8);
assert(abs(exposure.summary.totalNetEnergyResolved- ...
    sum(mesh.localGeometry.area(wall.validFace).* ...
    exposure.cumulativeNetHeat(wall.validFace,end)))<1e-7);
assert(wall.stagnationPoint.available);
assert(wall.stagnationPoint.temperature(end)>cfg.flight.Twall);
assert(isfield(wall,'coupling') && ~any(wall.lumpedFace));
assert(all(wall.oneDimensionalFace==wall.validFace));
% 耦合求解保存厚度温度，储热必须沿厚度积分，不能用表面温度代替。
depth=wall.oneDimensional.depth;
stored=material.density*material.specificHeat* ...
    trapz(depth,wall.oneDimensional.finalDepthTemperature-cfg.flight.Twall,2);
assert(max(abs(stored-exposure.cumulativeNetHeat(wall.validFace,end)))<1e-3);
tilted=flight; tilted.alpha=[0 10 0];
tiltedHistory=calcFlightHeatingHistory(mesh,tilted,cfg);
assert(tiltedHistory.stagnation.available(2));
assert(~tiltedHistory.axisymmetricMarchAppliedByStation(2));
assert(tiltedHistory.summary.noseFacesByStation(2)>1000);
% 低温冻结段补齐多数背风区域；低于 80 K 的残余面元仍明确留空。
assert(tiltedHistory.summary.validAcrossHistory>0.8*size(mesh.F,1));
assert(tiltedHistory.summary.validAcrossHistoryAreaFraction>0.8);
assert(tiltedHistory.summary.calculatedAreaFractionByStation(2)>0.9);
assert(tiltedHistory.summary.reasonAreaByStation(9,2)>0);
assert(tiltedHistory.summary.reasonAreaByStation(9,2)< ...
    0.05*sum(mesh.localGeometry.area));
tiltedWall=calcWallThermalHistory(mesh,tiltedHistory,material,cfg);
assert(~tiltedWall.wallThermal.stagnationPoint.available);
assert(contains(tiltedWall.wallThermal.stagnationPoint.status,'移动'));
assert(tiltedWall.wallThermal.summary.calculatedFaces>0.8*size(mesh.F,1));
fprintf('Flight heating history passed: valid across history %d/%d, nose %d faces.\n', ...
    history.summary.validAcrossHistory,size(mesh.F,1), ...
    history.summary.noseFacesByStation(1));
end
