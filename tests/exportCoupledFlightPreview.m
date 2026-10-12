function exportCoupledFlightPreview(outputFolder)
%EXPORTCOUPLEDFLIGHTPREVIEW 独立发布目录中复核新耦合交互并导出界面。
% 此辅助程序仅生成可复现的预览，不改变主程序默认飞行/材料设置。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg)); [raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
cfg.pressureGUI.visible='off';
fig=pressureGUI(mesh,cfg); cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
ui=fig.UserData;
ui.flightTable.Data=[0 30 5 0 0; .2 20 5 10 0; .4 30 5 0 0];
ui.maxFlightStep.Value=.2;
material=[2700 900 .003 180 .8 300];
for j=1:6, ui.wallFields(j).Value=material(j); end
ui.runFlightButton.ButtonPushedFcn(ui.runFlightButton,[]);
data=fig.UserData; w=data.wallThermal;
assert(isfield(w,'coupling') && ~w.stagnationPoint.available);
assert(nnz(w.validFace)>1000);
assert(data.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(any(contains(string(ui.flightDescription.Value),'实际热壁历程')));
% 点击最热的有效侧壁面元，核对当前时间来源，而不是只展示总算法。
valid=find(w.validFace & w.coupling.turbulentSourceIDHistory(:,end)>0);
[~,chosen]=max(w.temperature(valid,end)); id=valid(chosen);
ui.wallSurface.ButtonDownFcn(ui.wallSurface,struct( ...
    'IntersectionPoint',mesh.localGeometry.center(id,:)));
assert(fig.UserData.wallSelectedFace==id);
ui.wallTimeSlider.Value=numel(w.time);
ui.wallTimeSlider.ValueChangedFcn(ui.wallTimeSlider,[]);
assert(any(contains(string(ui.wallHistoryAxis.Title.String),'当前热流')));
assert(isequaln(ui.wallSurface.FaceVertexCData,w.temperature(:,end)));
ui.tabs.SelectedTab=ui.wallTab; drawnow;
exportapp(fig,fullfile(outputFolder,'coupled_flight_temperature_gui.png'));
ui.wallTimeSlider.Value=1; ui.wallTimeSlider.ValueChangedFcn(ui.wallTimeSlider,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData,w.temperature(:,1)));
ui.wallTimeSlider.Value=numel(w.time);
ui.wallMetric.Value=2; ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData, ...
    data.aerothermalExposure.cumulativeNetHeat(:,end)/1e6));
exportapp(fig,fullfile(outputFolder,'coupled_flight_exposure_gui.png'));
ui.tabs.SelectedTab=ui.flightTab; drawnow;
exportapp(fig,fullfile(outputFolder,'coupled_flight_coverage_gui.png'));
fprintf('Coupled flight GUI passed: %d faces, %.2f%% area, %d accepted times, residual %.3g J/m2.\n', ...
    nnz(w.validFace),100*data.aerothermalExposure.summary.resolvedAreaFraction, ...
    numel(w.time),data.aerothermalExposure.summary.maximumEnergyBalanceResidual);
end
