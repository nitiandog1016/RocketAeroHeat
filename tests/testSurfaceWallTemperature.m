function testSurfaceWallTemperature(outputFolder)
%TESTSURFACEWALLTEMPERATURE 验证非均匀壁温的沿程传递和鼻端焓差。
% 模块一：统一壁温场必须复现原标量调用，包括球鼻热流。
% 模块二：仅改变上游壁温，下游自身壁温不变，边界层仍须响应。
% 模块三：未知壁温留空；非均匀热边界不能误用作旧统一初温历程。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg)); [raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
fs=calcFreestream(cfg.flight);
mesh=modifiedNewtonianPressure(mesh,fs,cfg);
mesh=findStagnationRegion(mesh,fs,cfg);
mesh=solveExternalFlow(mesh,fs,cfg); mesh=traceSurfaceStreamlines(mesh,cfg);
mesh=calcSurfaceHeating(mesh,fs,cfg);
n=size(mesh.F,1); initial=cfg.flight.Twall;
uniform=recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,initial*ones(n,1));
assert(isequaln(uniform.heating.heatFlux,mesh.heating.heatFlux));
assert(isequaln(uniform.heating.turbulentBoundaryLayer.shapeFactor, ...
    mesh.heating.turbulentBoundaryLayer.shapeFactor));

%% 同一温度场用标量和向量输入，不能因驻点参考温度不同而改变鼻端结果
hotCfg=cfg; hotCfg.flight.Twall=800;
scalar=calcSurfaceHeating(mesh,fs,hotCfg);
vector=recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,800*ones(1,n));
active=scalar.heating.eligibleFace & vector.heating.eligibleFace;
assert(isequal(scalar.heating.eligibleFace,vector.heating.eligibleFace));
assert(all(abs(scalar.heating.heatFlux(active)-vector.heating.heatFlux(active)) ...
    <1e-10*max(1,abs(scalar.heating.heatFlux(active)))));
assert(nnz(vector.heating.noseFace)>0);
assert(isequal(size(vector.heating.heatFlux),[n 1]));
assert(vector.stagnationHeating.wallTemperature==initial);

%% 上游热斑会影响下游积分状态，不能退化成仅由当地温度决定的局部缩放
distance=mesh.surfaceStreamlines.distance;
heated=distance>0.4 & distance<0.9;
temperature=initial*ones(n,1); temperature(heated)=900;
changed=recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,temperature);
b0=mesh.heating.turbulentBoundaryLayer; b1=changed.heating.turbulentBoundaryLayer;
downstream=distance>1 & mesh.heating.turbulentFace & ...
    changed.heating.turbulentFace & b0.sourceID==2 & b1.sourceID==2;
responding=downstream & abs(b1.momentumThickness-b0.momentumThickness)>1e-9;
assert(nnz(responding)>100);
assert(all(changed.heating.wallTemperatureByFace(responding)==initial));
assert(any(abs(changed.heating.heatFlux(responding)-mesh.heating.heatFlux(responding))>1e-3));
assert(isequaln(mesh.externalFlow,changed.externalFlow));
assert(isequaln(mesh.surfaceStreamlines,changed.surfaceStreamlines));
assert(all(changed.boundaryLayer.wallTemperatureByFace==temperature));

%% 未知壁温以及超范围壁温均不能被球鼻延拓或平板回退填补
missing=temperature;
noseID=find(mesh.heating.noseFace,1); sideID=find(responding,1);
missing([noseID sideID])=NaN;
gap=recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,missing);
assert(all(isnan(gap.heating.heatFlux([noseID sideID]))));
assert(~any(gap.heating.eligibleFace([noseID sideID])));
assert(gap.heatingCoverage.reasonCode(sideID)==13);
allMissing=recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,nan(n,1));
assert(~any(allMissing.heating.eligibleFace));
assert(all(isnan(allMissing.heating.heatFlux)));
% 给定 Tw=Taw 的鼻端面元，q=0 但换热系数应有有限极限。
equilibrium=temperature;
equilibrium(noseID)=mesh.stagnationHeating.adiabaticWallTemperature;
neutral=recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,equilibrium);
assert(neutral.heating.heatFlux(noseID)==0);
assert(isfinite(neutral.heating.heatTransferCoefficient(noseID)));

%% 尺寸、物理取值以及旧壁温求解接口的误用应给出清楚错误
for bad={ones(n-1,1),zeros(n,1),inf(n,1)}
    failed=false;
    try, recomputeSurfaceHeatingAtWallTemperature(mesh,fs,cfg,bad{1}); catch, failed=true; end
    assert(failed);
end
material=struct('density',1000,'specificHeat',1000,'thickness',0.01, ...
    'thermalConductivity',10,'emissivity',0,'radiationTemperature',300);
failed=false;
try, calcWallThermalResponse(changed,[0 1],material,cfg);
catch problem, failed=strcmp(problem.identifier,'RocketAeroHeat:WallThermal:StaleHeating'); end
assert(failed);
if nargin>0
    displayCfg=cfg; displayCfg.pressureGUI.visible='off';
    displayCfg.heating.wallTemperatureByFace=temperature;
    fig=pressureGUI(changed,displayCfg);
    cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
    ui=fig.UserData; ui.tabs.SelectedTab=ui.heatTab;
    ui.heatChoice.Value=20;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]); drawnow;
    assert(isequal(ui.heatSurface.FaceVertexCData,temperature));
    exportapp(fig,fullfile(outputFolder,'prescribed_wall_temperature_gui.png'));
    ui.heatChoice.Value=1;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]); drawnow;
    exportapp(fig,fullfile(outputFolder,'nonuniform_wall_heating_gui.png'));
end
fprintf('Surface wall field passed: %d downstream faces respond to upstream heating; uniform, nose and missing-state checks passed.\n',nnz(responding));
end
