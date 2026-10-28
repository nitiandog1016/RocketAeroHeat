function testSurfaceHeating(outputFolder)
% 真实箭体回归：检查热力极限、算法标注、未计算区域及统一窗口显示。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);
[raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
fs=calcFreestream(cfg.flight);
mesh=modifiedNewtonianPressure(mesh,fs,cfg);
mesh=findStagnationRegion(mesh,fs,cfg);
mesh=calcSurfaceHeating(mesh,fs,cfg);
h=mesh.heating; e=mesh.boundaryLayerEdge;
c=mesh.heatingCoverage;
% 此直接调用未追踪流线，应明确采用局部平板回退。
assert(isempty(h.turbulentBoundaryLayer));
assert(contains(h.turbulentMethod,'局部平板'));
assert(isequal(c.reasonCode==0,h.eligibleFace));
assert(sum(c.countByReason)==size(mesh.F,1));
assert(abs(sum(c.areaByReason)-c.summary.totalArea)<1e-10);
assert(c.summary.computedAreaFraction>0.9);
assert(contains(h.method,'Eckert') && contains(e.method,'等熵'));
assert(all(isnan(h.heatFlux(h.unresolvedFace))));
assert(nnz(h.eligibleFace)>100 && nnz(h.unresolvedFace)>0);
assert(all(isfinite(h.heatFlux(h.eligibleFace))));
assert(nnz(h.noseFace)>10);
assert(all(h.faceMethodID(h.eligibleFace)>0));
assert(all(h.noseFace<=h.noseDistribution.faceMask));
assert(mesh.continuumRegime.freestreamKnudsen<0.01);
assert(contains(h.continuumMethod,'平均自由程'));
assert(all(abs(h.heatFlux(h.noseFace)- ...
    h.noseDistribution.stagnationHeatFlux.* ...
    cosd(h.noseDistribution.polarAngleDeg(h.noseFace)))<1e-8));
assert(all(abs(h.heatFlux(h.laminarFace)-h.laminarHeatFlux(h.laminarFace))<1e-9));
assert(all(abs(h.heatFlux(h.turbulentFace)-h.turbulentHeatFlux(h.turbulentFace))<1e-9));
% Cp=0 处，边界状态必须退回自由来流。
[~,id]=min(abs(mesh.pressure.Cp));
assert(abs(e.Mach(id)-fs.Mach)<1e-10);
assert(abs(e.temperature(id)-fs.temperature)<1e-9);
assert(abs(e.speed(id)-fs.speed)<1e-9);
% 导入逐面元三维速度后，热流改用有向图流动距离，不能再默默退回子午弧长。
flow=repmat(fs.flowVelocityVector,size(mesh.F,1),1);
normal=mesh.localGeometry.normal;
flow=flow-sum(flow.*normal,2).*normal;
data=struct('faceCenter',mesh.localGeometry.center, ...
    'faceNormal',normal,'pressure',fs.pressure*ones(size(mesh.F,1),1), ...
    'temperature',fs.temperature*ones(size(mesh.F,1),1), ...
    'density',fs.pressure/(fs.R*fs.temperature)*ones(size(mesh.F,1),1), ...
    'velocity',flow,'method','测试用解析投影流','caseID','test-projected-flow', ...
    'freestream',fs);
imported=importExternalFlowState(mesh,fs,data);
imported=traceSurfaceStreamlines(imported);
imported=calcSurfaceHeating(imported,fs,cfg);
assert(imported.surfaceStreamlines.summary.connectedFaces>100);
assert(imported.heating.summary.calculatedFaces>0);
assert(contains(imported.heating.streamlineDistanceAssumption,'Dijkstra'));
assert(isstruct(imported.heating.turbulentBoundaryLayer));
assert(contains(imported.heating.turbulentMethod,'Ludwieg-Tillmann'));
assert(all(isfinite(imported.heating.turbulentBoundaryLayer.skinFrictionCoefficient( ...
    imported.heating.turbulentFace))));
% 受控端到端回归：上游外缘热力占位不是有效流场；即使下游
% Reynolds较小且当地状态有效，也不能越过该未知状态求热流。
parent=imported.surfaceStreamlines.predecessorFace;
candidate=find(imported.heating.turbulentFace & parent>0 & ...
    imported.localGeometry.center(:,1)<-1);
candidate=candidate(imported.heating.turbulentFace(parent(candidate)));
assert(~isempty(candidate));
child=candidate(1); bridge=parent(child);
gap=imported; edge=gap.externalFlow.edge;
edge.density(child)=2e5*edge.dynamicViscosity(child)/ ...
    (edge.speed(child)*gap.surfaceStreamlines.distance(child));
edge.pressure(child)=edge.density(child)*fs.R*edge.temperature(child);
edge.thermalExpansionOutOfRangeFace=false(size(mesh.F,1),1);
edge.thermalExpansionOutOfRangeFace(bridge)=true;
gap.externalFlow.edge=edge;
gap=calcSurfaceHeating(gap,fs,cfg);
assert(isnan(gap.heating.heatFlux(bridge)));
assert(gap.heating.edgeReynolds(child)<cfg.heating.transitionRe);
assert(gap.continuumRegime.directStateUnresolvedFace(bridge));
assert(gap.continuumRegime.downstreamStateUnresolvedFace(child));
assert(isnan(gap.continuumRegime.knudsen(bridge)));
assert(isnan(gap.heating.heatFlux(child)) && ~gap.heating.eligibleFace(child));
% 人工受控减速外缘仅用于筛选回归，不声称它是实际箭体外流解。
% 强减速后恢复速度，已分离路径不能再由平板回退或鼻端延拓补热流。
decelerated=imported;
x=decelerated.localGeometry.center(:,1);
factor=1-0.8*min(1,max(0,(-x-0.9)/0.15));
factor(x<-1.4)=1;
decelerated.externalFlow.edge.speed=imported.externalFlow.edge.speed.*factor;
decelerated.externalFlow.edge.velocityVector=imported.externalFlow.edge.velocityVector.*factor;
decelerated.externalFlow.edge.Mach=imported.externalFlow.edge.Mach.*factor;
separationCfg=cfg; separationCfg.heating.transitionRe=1;
decelerated=calcSurfaceHeating(decelerated,fs,separationCfg);
state=decelerated.heating.turbulentBoundaryLayer;
assert(nnz(state.possibleSeparationFace)>0);
blocked=state.separationUnresolvedFace;
assert(any(blocked & x<-1.4));
assert(all(isnan(decelerated.heating.heatFlux(blocked))));
assert(all(isnan(decelerated.heating.turbulentHeatFlux(blocked))));
assert(~any(decelerated.heating.eligibleFace(blocked)));
assert(~any(decelerated.heating.noseFace(blocked)));
assert(any(decelerated.heatingCoverage.reasonCode==11));
assert(isequal(decelerated.heatingCoverage.reasonCode==0, ...
    decelerated.heating.eligibleFace));
% 原固定 H 模型在同一人工工况下不会检测湍流分离，作为回归对照。
legacyCfg=separationCfg; legacyCfg.heating.turbulentModel='fixedShape';
legacy=calcSurfaceHeating(decelerated,fs,legacyCfg);
assert(any(isfinite(legacy.heating.heatFlux(blocked))));
stale=fs; stale.alpha=10;
failed=false;
try, calcBoundaryLayerEdgeState(mesh,stale,cfg); catch, failed=true; end
assert(failed);
% 壁温高于所有恢复温度时，对流热流方向应反转。
hot=cfg; hot.flight.Twall=1e5;
hot.heating.airThermalModel='caloricallyPerfect';
warm=calcSurfaceHeating(mesh,fs,hot);
assert(all(warm.heating.heatFlux(warm.heating.eligibleFace)<0));
hot.heating.airThermalModel='thermallyPerfect';
outside=calcSurfaceHeating(mesh,fs,hot);
assert(~any(outside.heating.eligibleFace));
assert(any(outside.heating.thermalPropertyUnavailableFace));
% 80 km 已超出默认连续介质阈值，热流和驻点经验式必须留空。
highCfg=cfg; highCfg.flight.altitude=80000;
highFs=calcFreestream(highCfg.flight);
highMesh=modifiedNewtonianPressure(mesh,highFs,highCfg);
highMesh=findStagnationRegion(highMesh,highFs,highCfg);
highMesh=calcSurfaceHeating(highMesh,highFs,highCfg);
assert(highMesh.continuumRegime.freestreamKnudsen>0.01);
assert(~highMesh.stagnationHeating.available);
assert(~any(highMesh.heating.eligibleFace));
assert(any(highMesh.heatingCoverage.reasonCode==10));
% 同一压力窗口中核对选择、着色和算法说明。
cfg.pressureGUI.visible='off';
fig=pressureGUI(mesh,cfg); cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
ui=fig.UserData;
assert(ui.maxFlightStep.Value==1);
assert(strcmp(ui.loadsPanel.Visible,'off'));
assert(isequal(ui.contentGrid.ColumnWidth,{'1x',0}));
ui.loadToggle.ButtonPushedFcn(ui.loadToggle,[]);
assert(strcmp(ui.loadsPanel.Visible,'on'));
ui.loadToggle.ButtonPushedFcn(ui.loadToggle,[]);
assert(strcmp(ui.loadsPanel.Visible,'off'));
assert(ui.axes.CameraViewAngle<10);
assert(ui.heatAxis.CameraViewAngle<8);
if nargin>0
    ui.tabs.SelectedTab=ui.mapTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'pressure_gui.png'));
end
assert(isgraphics(ui.heatTab) && any(contains(string(ui.heatDescription.Value),'Reynolds')));
assert(any(contains(string(ui.heatDescription.Value),'转捩状态连续性')));
assert(ui.heatScale.Value==2);
assert(strcmp(ui.heatDescription.Visible,'off'));
assert(ui.heatGrid.RowHeight{4}==0);
assert(contains(ui.heatSummary.Text,'Eckert'));
if nargin>0
    ui.tabs.SelectedTab=ui.heatTab; drawnow;
    assert(ui.heatAxis.Position(3)>0.7);
    exportapp(fig,fullfile(outputFolder,'heating_initial_gui.png'));
    ui.tabs.SelectedTab=ui.mapTab; drawnow;
end
robustLimits=ui.heatAxis.CLim;
ui.heatScale.Value=1;
ui.heatScale.ValueChangedFcn(ui.heatScale,[]);
assert(ui.heatAxis.CLim(1)<=robustLimits(1));
assert(ui.heatAxis.CLim(2)>=robustLimits(2));
assert(isequaln(ui.heatSurface.FaceVertexCData,ui.heating.heatFlux/1e4));
ui.heatScale.Value=2;
ui.heatScale.ValueChangedFcn(ui.heatScale,[]);
ui.heatDetailsButton.ButtonPushedFcn(ui.heatDetailsButton,[]);
assert(strcmp(ui.heatDescription.Visible,'on'));
assert(ui.heatGrid.RowHeight{4}==190);
ui.heatDetailsButton.ButtonPushedFcn(ui.heatDetailsButton,[]);
assert(strcmp(ui.heatDescription.Visible,'off'));
assert(ui.heatGrid.RowHeight{4}==0);
assert(any(contains(string(ui.heatDescription.Value),'驻点法向激波')));
assert(any(contains(string(ui.heatDescription.Value),'Billig')));
assert(strcmp(ui.shockPreviewCheckbox.Enable,'on'));
ui.shockPreviewCheckbox.Value=true;
ui.shockPreviewCheckbox.ValueChangedFcn(ui.shockPreviewCheckbox,[]);
assert(strcmp(ui.shockPreview.Visible,'on'));
shockDiameter=max(ui.shockPreview.Vertices(:,2))- ...
    min(ui.shockPreview.Vertices(:,2));
assert(diff(xlim(ui.heatAxis))<4*shockDiameter);
ui.tabs.SelectedTab=ui.heatTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'bow_shock_preview_gui.png'));
end
ui.shockPreviewCheckbox.Value=false;
ui.shockPreviewCheckbox.ValueChangedFcn(ui.shockPreviewCheckbox,[]);
assert(strcmp(ui.shockPreview.Visible,'off'));
assert(ui.heatAxis.CameraViewAngle<8);
ui.noseZoomButton.ButtonPushedFcn(ui.noseZoomButton,[]);
assert(diff(xlim(ui.heatAxis))<4*shockDiameter);
ui.resetButton.ButtonPushedFcn(ui.resetButton,[]);
assert(diff(xlim(ui.heatAxis))>2);
assert(ui.heatAxis.CameraViewAngle<8);
% 在鼻端局部视图切换显示倍率时，应恢复整箭体，不保留旧放大相机。
ui.noseZoomButton.ButtonPushedFcn(ui.noseZoomButton,[]);
ui.displayScale.Value=5;
ui.displayScale.ValueChangedFcn(ui.displayScale,[]);
assert(diff(xlim(ui.heatAxis))>2);
assert(isequal(ui.heatAxis.DataAspectRatio,[5 1 1]));
assert(isequaln(ui.heatSurface.FaceVertexCData,ui.heating.heatFlux/1e4));
ui.displayScale.Value=3;
ui.displayScale.ValueChangedFcn(ui.displayScale,[]);
assert(isequaln(ui.heatSurface.FaceVertexCData,ui.heating.heatFlux/1e4));
ui.heatChoice.Value=10;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequal(ui.heatSurface.FaceVertexCData,double(ui.heating.faceMethodID)));
ui.heatChoice.Value=11;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequal(ui.heatSurface.FaceVertexCData, ...
    double(ui.heatingCoverage.reasonCode)));
ui.heatChoice.Value=14;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequal(ui.heatSurface.FaceVertexCData, ...
    double(ui.continuumRegime.applicabilityID)));
ui.heatChoice.Value=15;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequal(ui.heatSurface.FaceVertexCData, ...
    double(ui.heating.turbulentBoundaryLayer.sourceID)));
if nargin>0
    ui.tabs.SelectedTab=ui.heatTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'turbulent_sources_gui.png'));
end
ui.heatChoice.Value=16;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequaln(ui.heatSurface.FaceVertexCData, ...
    ui.heating.turbulentBoundaryLayer.skinFrictionCoefficient));
if nargin>0
    ui.tabs.SelectedTab=ui.heatTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'turbulent_cf_gui.png'));
end
ui.heatChoice.Value=17;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequaln(ui.heatSurface.FaceVertexCData, ...
    ui.heating.turbulentBoundaryLayer.adversePressureGradientParameter));
ui.heatChoice.Value=18;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequaln(ui.heatSurface.FaceVertexCData, ...
    ui.heating.turbulentBoundaryLayer.shapeFactor));
assert(contains(ui.heatSummary.Text,'Head'));
if nargin>0
    ui.tabs.SelectedTab=ui.heatTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'turbulent_shape_factor_gui.png'));
end
ui.heatChoice.Value=19;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
state=ui.heating.turbulentBoundaryLayer;
expected=zeros(size(ui.heating.heatFlux));
expected(state.separationUnresolvedFace)=2;
expected(state.possibleSeparationFace)=1;
expected(state.numericalUnresolvedFace)=3;
assert(isequal(ui.heatSurface.FaceVertexCData,expected));
if nargin>0
    ui.tabs.SelectedTab=ui.heatTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'turbulent_separation_gui.png'));
end
ui.heatChoice.Value=20;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequaln(ui.heatSurface.FaceVertexCData,ui.heating.wallTemperatureByFace));
assert(contains(ui.heatSummary.Text,'逐面元壁温'));
assert(any(contains(string(ui.heatDescription.Value),'Ludwieg-Tillmann')));
assert(any(contains(string(ui.heatDescription.Value),'连续介质筛选')));
if nargin>0
    ui.tabs.SelectedTab=ui.heatTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'continuum_regime_gui.png'));
end
assert(any(contains(string(ui.heatDescription.Value),'热流面积覆盖')));
ui.tabs.SelectedTab=ui.heatTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'heating_coverage_gui.png'));
end
ui.heatChoice.Value=10;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
ui.tabs.SelectedTab=ui.heatTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'nose_heat_method_gui.png'));
end
assert(contains(ui.boundaryLayerEdge.method,'内置工程外流'));
ui.heatChoice.Value=4;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequaln(ui.heatSurface.FaceVertexCData,ui.boundaryLayerEdge.temperature));
ui.heatChoice.Value=6;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequal(ui.heatSurface.FaceVertexCData, ...
    double(ui.boundaryLayerEdge.faceMethodID)));
ui.heatChoice.Value=12;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
assert(isequal(ui.heatSurface.FaceVertexCData, ...
    double(ui.boundaryLayerEdge.thermalPropertyModeByFace)));
ui.heatChoice.Value=7;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
expectedTransition=ui.heating.transition.intermittency;
expectedTransition(ui.heating.transition.onsetFace)=2;
assert(isequaln(ui.heatSurface.FaceVertexCData,expectedTransition));
ui.tabs.SelectedTab=ui.heatTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'external_methods_gui.png'));
end
ui.heatChoice.Value=1;
ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
ui.tabs.SelectedTab=ui.heatTab; drawnow;
assert(ui.heatAxis.Position(3)>0.7);
if nargin>0
    assert(ui.heatAxis.Position(3)>0.7);
    assert(ui.heatColorbar.Position(1)>0.85);
    exportapp(fig,fullfile(outputFolder,'heating_gui.png'));
end
assert(strcmp(ui.streamlineCheckbox.Enable,'on'));
ui.streamlineCheckbox.Value=true;
ui.streamlineCheckbox.ValueChangedFcn(ui.streamlineCheckbox,[]);
viewDirection=ui.heatAxis.CameraPosition-ui.heatAxis.CameraTarget;
viewDirection=viewDirection/norm(viewDirection);
visibleNormal=mesh.localGeometry.normal*viewDirection'>0.5;
pathIDs=find(ui.surfaceStreamlines.connectedFace & ...
    mesh.localGeometry.center(:,1)<-1 & visibleNormal);
assert(~isempty(pathIDs));
picked=pathIDs(1);
hit=struct('IntersectionPoint',mesh.localGeometry.center(picked,:));
ui.heatSurface.ButtonDownFcn(ui.heatSurface,hit);
assert(strcmp(ui.streamlineLine.Visible,'on'));
assert(numel(ui.streamlineLine.XData)>2);
assert(any(contains(string(ui.heatDescription.Value),'准流线：面元')));
ui.tabs.SelectedTab=ui.heatTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'streamline_gui.png'));
end
ui.streamlineCheckbox.Value=false;
ui.streamlineCheckbox.ValueChangedFcn(ui.streamlineCheckbox,[]);
assert(strcmp(ui.streamlineLine.Visible,'off'));
wallInput=[1000 1000 0.01 10 0.8 300 10 11];
for k=1:8, ui.wallFields(k).Value=wallInput(k); end
ui.runWallButton.ButtonPushedFcn(ui.runWallButton,[]);
assert(isfield(fig.UserData,'wallThermal'));
assert(fig.UserData.wallThermal.summary.calculatedFaces>0);
assert(isfield(fig.UserData.wallThermal,'coupling'));
assert(isequal(fig.UserData.wallThermal.time, ...
    fig.UserData.wallThermal.stagnationPoint.time));
assert(fig.UserData.wallThermal.stagnationPoint.available);
assert(any(contains(string(ui.wallDescription.Value),'驻点单点末温')));
assert(any(fig.UserData.wallThermal.temperature(:,end)>300,'all'));
% 同工况重建压力页时，已算壁温和热量应继续显示。
savedMesh=evalin('base','mesh');
restored=pressureGUI(savedMesh,cfg);
assert(isfield(restored.UserData,'wallThermal') && ...
    isfield(restored.UserData,'aerothermalExposure'));
assert(strcmp(restored.UserData.wallTimeSlider.Enable,'on'));
delete(restored);
wall=fig.UserData.wallThermal;
exposure=fig.UserData.aerothermalExposure;
assert(exposure.summary.calculatedFaces==wall.summary.calculatedFaces);
assert(exposure.summary.resolvedAreaFraction>0.9);
assert(exposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(any(contains(string(ui.wallDescription.Value),'热量积分')));
assert(strcmp(ui.wallTimeSlider.Enable,'on'));
ui.wallTimeSlider.Value=1;
ui.wallTimeSlider.ValueChangedFcn(ui.wallTimeSlider,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData,wall.temperature(:,1)));
assert(fig.UserData.wallTimeIndex==1);
id=find(wall.validFace,1);
click=struct('IntersectionPoint',mesh.localGeometry.center(id,:));
ui.wallSurface.ButtonDownFcn(ui.wallSurface,click);
assert(fig.UserData.wallSelectedFace==id);
ui.wallTimeSlider.Value=numel(wall.time);
ui.wallTimeSlider.ValueChangedFcn(ui.wallTimeSlider,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData,wall.temperature(:,end)));
assert(numel(findobj(ui.wallHistoryAxis,'Type','line'))>=2);
% 逐项核对滑块对应的热量与峰值图，不能把全历程峰值误标为中途时刻。
middle=ceil(numel(wall.time)/2);
ui.wallTimeSlider.Value=middle;
ui.wallTimeSlider.ValueChangedFcn(ui.wallTimeSlider,[]);
ui.wallMetric.Value=2;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData, ...
    exposure.cumulativeNetHeat(:,middle)/1e6));
% v37 明确显示材料净入壁热量，区别于尚未扣除内侧散热的外侧净输入。
assert(contains(ui.wallAxis.Title.String{1},'累计净入壁'));
ui.wallMetric.Value=3;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData, ...
    exposure.cumulativeConvectiveHeat(:,middle)/1e6));
ui.wallMetric.Value=5;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData, ...
    exposure.cumulativeRadiativeHeat(:,middle)/1e6));
ui.wallMetric.Value=4;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData, ...
    max(wall.convectiveHeatFlux(:,1:middle),[],2)/1e4));
assert(contains(ui.wallAxis.Title.String{1},'截至当前时刻峰值'));
ui.wallTimeSlider.Value=numel(wall.time);
ui.wallTimeSlider.ValueChangedFcn(ui.wallTimeSlider,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData, ...
    exposure.peakConvectiveHeatFlux/1e4));
ui.tabs.SelectedTab=ui.wallTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'wall_peak_gui.png'));
end
ui.wallMetric.Value=2;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
ui.tabs.SelectedTab=ui.wallTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'wall_exposure_gui.png'));
end
ui.wallMetric.Value=1;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
if nargin>0
    ui.tabs.SelectedTab=ui.wallTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'coupled_wall_temperature_gui.png'));
end
ui.addFlightButton.ButtonPushedFcn(ui.addFlightButton,[]);
assert(size(ui.flightTable.Data,1)==4);
ui.removeFlightButton.ButtonPushedFcn(ui.removeFlightButton,[]);
assert(size(ui.flightTable.Data,1)==3);
ui.flightTable.Data=[0 30 5 0 0;0.2 20 5 10 0;0.4 30 5 0 0];
ui.tabs.SelectedTab=ui.flightTab;
drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'flight_inputs_gui.png'));
end
ui.runFlightButton.ButtonPushedFcn(ui.runFlightButton,[]);
assert(isfield(fig.UserData,'flightHeatingHistory'));
assert(any(contains(string(ui.flightDescription.Value),'各站有解面积')));
assert(any(contains(string(ui.flightDescription.Value), ...
    '第 1 计算站') & contains(string(ui.flightDescription.Value),'已采用')));
assert(any(contains(string(ui.flightDescription.Value), ...
    '第 2 计算站') & contains(string(ui.flightDescription.Value),'未采用')));
assert(fig.UserData.wallThermal.summary.calculatedFaces>1000);
assert(all(ismember([0 .2 .4],fig.UserData.wallThermal.time)));
assert(isfield(fig.UserData.wallThermal,'coupling'));
assert(~fig.UserData.wallThermal.stagnationPoint.available);
assert(isequal(ui.tabs.SelectedTab,ui.wallTab));
assert(any(contains(string(ui.wallDescription.Value),'逐时刻外流热边界')));
flightExposure=fig.UserData.aerothermalExposure;
% 历程也不能为了旧80%比例跨越未知上游；核对真实面积及逐站原因。
assert(flightExposure.summary.resolvedAreaFraction>0 && ...
    flightExposure.summary.resolvedAreaFraction<1);
assert(abs(flightExposure.summary.resolvedArea- ...
    sum(mesh.localGeometry.area(fig.UserData.wallThermal.validFace)))<1e-12);
reasonIndex=find(strcmp(fig.UserData.flightHeatingHistory.reasonNames, ...
    '上游外缘热力状态未解/域外占位影响'),1);
assert(~isempty(reasonIndex) && ...
    fig.UserData.flightHeatingHistory.summary.reasonAreaByStation(reasonIndex,2)>0);
assert(fig.UserData.flightHeatingHistory.unknownThermalFacesByStation(2)>0);
assert(fig.UserData.flightHeatingHistory.upstreamUnknownThermalFacesByStation(2)>0);
assert(fig.UserData.flightHeatingHistory.rarefiedFacesByStation(2)==0, ...
    '该工况域外占位不能计入真实稀薄流面数。');
fprintf('Current flight-history heat area %.9f%%; upstream thermal unknown remains missing.\n', ...
    100*flightExposure.summary.resolvedAreaFraction);
assert(flightExposure.summary.maximumEnergyBalanceResidual<1e-3);
ui.wallMetric.Value=2;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData, ...
    flightExposure.cumulativeNetHeat(:,end)/1e6));
ui.tabs.SelectedTab=ui.flightTab; drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'flight_coverage_gui.png'));
end
ui.tabs.SelectedTab=ui.wallTab;
drawnow;
if nargin>0
    exportapp(fig,fullfile(outputFolder,'flight_exposure_gui.png'));
end
% GUI 的最大计算间隔直接控制内部站点数，而不是仅改变显示文字。
ui.tabs.SelectedTab=ui.flightTab;
ui.maxFlightStep.Value=0.1;
ui.runFlightButton.ButtonPushedFcn(ui.runFlightButton,[]);
assert(fig.UserData.flightHeatingHistory.summary.inputStations==3);
assert(fig.UserData.flightHeatingHistory.summary.evaluationStations==5);
assert(all(ismember(fig.UserData.flightHeatingHistory.time,fig.UserData.wallThermal.time)));
assert(any(contains(string(ui.flightDescription.Value),'输入 3 站 → 内部计算 5 站')));
ui.wallMetric.Value=1;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
ui.runWallButton.ButtonPushedFcn(ui.runWallButton,[]);
assert(~isfield(fig.UserData,'flightHeatingHistory'));
assert(~isfield(evalin('base','mesh'),'flightHeatingHistory'));
fprintf('Heating passed: %d/%d faces, q=[%.3f, %.3f] W/cm2.\n', ...
    h.summary.calculatedFaces,numel(h.heatFlux), ...
    h.summary.minimumHeatFlux/1e4,h.summary.maximumHeatFlux/1e4);
end
