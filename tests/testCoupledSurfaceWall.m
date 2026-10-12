function testCoupledSurfaceWall()
%TESTCOUPLEDSURFACEWALL 检查内置外流、沿程边界层和导热的实际接通。
% 模块作用：小时间段真实箭体端到端检查，核对末时刻热流确实来自
% 末时刻温度场的整条边界层重算，而不是旧局部幂指数缩放。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg)); [raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
fs=calcFreestream(cfg.flight); mesh=modifiedNewtonianPressure(mesh,fs,cfg);
mesh=findStagnationRegion(mesh,fs,cfg); mesh=solveExternalFlow(mesh,fs,cfg);
mesh=traceSurfaceStreamlines(mesh,cfg); mesh=calcSurfaceHeating(mesh,fs,cfg);
mat=struct('density',2700,'specificHeat',900,'thickness',0.003, ...
    'thermalConductivity',180,'emissivity',0.8,'radiationTemperature',300);
tic; out=calcWallThermalResponse(mesh,[0 0.05 0.1],mat,cfg); elapsed=toc;
w=out.wallThermal;
assert(isfield(w,'coupling') && contains(w.method,'Picard'));
assert(nnz(w.validFace)>1000 && ~any(w.lumpedFace));
assert(all(w.oneDimensionalFace==w.validFace));
again=recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,w.temperature(:,end));
actual=again.heating.heatFlux;
assert(all(isfinite(actual(w.validFace))));
assert(all(abs(actual(w.validFace)-w.convectiveHeatFlux(w.validFace,end)) ...
    <=1e-9*max(1,abs(actual(w.validFace)))));
assert(all(w.coupling.faceMethodIDHistory(w.validFace,end)== ...
    again.heating.faceMethodID(w.validFace)));
out=calcAerothermalExposure(out);
assert(out.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
oldCfg=cfg; oldCfg.heating.wallCoupling='local';
old=calcWallThermalResponse(mesh,[0 0.05 0.1],mat,oldCfg);
assert(~isfield(old.wallThermal,'coupling'));
assert(w.stagnationPoint.available);
% 驻点曲线与面元共用已接受时间站，自动子步不能只作用于面元矩阵。
assert(isequal(w.time,w.stagnationPoint.time));
assert(numel(w.stagnationPoint.temperature)==size(w.temperature,2));
fprintf('Coupled surface wall passed: %d faces, %d lost, %.2f s, max energy residual %.3g J/m2.\n', ...
    nnz(w.validFace),nnz(w.coupling.lostFace),elapsed, ...
    out.aerothermalExposure.summary.maximumEnergyBalanceResidual);
end
