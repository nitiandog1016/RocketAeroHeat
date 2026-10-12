function testCoupledFlightWall()
%TESTCOUPLEDFLIGHTWALL 真实 STL 的恒定/变工况沿程壁温耦合回归。
% 检查恒定轨迹与同外流固定工况一致，变工况末步热流确实来自
% 当前飞行条件及非均匀壁温，厚度储热守恒且驻点迁移不混成材料点。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg)); [raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
mat=struct('density',2700,'specificHeat',900,'thickness',0.003, ...
    'thermalConductivity',180,'emissivity',0.8,'radiationTemperature',300);
t=[0 .02 .04];
f=struct('time',t,'altitude',cfg.flight.altitude*ones(1,3), ...
    'Mach',cfg.flight.Mach*ones(1,3),'alpha',zeros(1,3),'beta',zeros(1,3));
h=calcFlightHeatingHistory(mesh,f,cfg);
out=calcWallThermalHistory(mesh,h,mat,cfg); w=out.wallThermal;
assert(isfield(w,'coupling') && w.coupling.externalStateSolves==1);
assert(w.stagnationPoint.available && isequal(w.stagnationPoint.time,w.time));
condition=cfg.flight; condition.alpha=0; condition.beta=0;
station=calcFlightStationState(mesh,condition,cfg);
fixed=calcWallThermalResponse(station,t,mat,cfg);
assert(isequal(w.validFace,fixed.wallThermal.validFace));
assert(max(abs(w.temperature(w.validFace,:)- ...
    fixed.wallThermal.temperature(w.validFace,:)),[],'all')<1e-8);
assert(max(abs(w.oneDimensional.finalDepthTemperature- ...
    fixed.wallThermal.oneDimensional.finalDepthTemperature),[],'all')<1e-8);
f.altitude=f.altitude+[0 -100 -200]; f.alpha=[0 5 10];
h=calcFlightHeatingHistory(mesh,f,cfg);
out=calcWallThermalHistory(mesh,h,mat,cfg); w=out.wallThermal;
assert(nnz(w.validFace)>1000 && ~w.stagnationPoint.available);
assert(contains(w.stagnationPoint.status,'移动'));
assert(w.coupling.externalStateSolves>=3);
assert(all(ismember(t,w.time)));
condition.altitude=f.altitude(end); condition.alpha=f.alpha(end);
station=calcFlightStationState(mesh,condition,cfg);
again=recomputeSurfaceHeatingAtWallTemperature(station, ...
    station.externalFlow.freestream,cfg,w.temperature(:,end));
assert(max(abs(again.heating.heatFlux(w.validFace)- ...
    w.convectiveHeatFlux(w.validFace,end)))<1e-6);
for j=1:numel(w.time)
    m=w.coupling.boundaryMetadata{j};
    assert(m.time==w.time(j));
    assert(abs(m.flight.alpha-interp1(t,f.alpha,w.time(j)))<1e-12);
end
out=calcAerothermalExposure(out);
assert(out.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
localCfg=cfg; localCfg.heating.wallCoupling='local';
local=calcWallThermalHistory(mesh,h,mat,localCfg);
assert(~isfield(local.wallThermal,'coupling'));
assert(~local.wallThermal.stagnationPoint.available);
bad=mesh; bad.F([1 2],:)=bad.F([2 1],:);
try
    calcWallThermalHistory(bad,h,mat,cfg);
    error('Test:MissingError','应拒绝面元重排');
catch problem
    assert(strcmp(problem.identifier,'RocketAeroHeat:FlightCoupling:Mesh'));
end
fprintf('Coupled flight passed: fixed equivalence, changing state, final heat, energy and moving stagnation.\n');
end
