function testLayeredWallGUI(outputFolder)
%TESTLAYEREDWALLGUI 真实箭体的分层编辑、恢复、固定/飞行耦合和三维显示。
% 以下两层及温度表都是人工验证数据，不是实际热防护材料建议。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg)); [raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
cfg.pressureGUI.visible='off';
fig=pressureGUI(mesh,cfg); cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
ui=fig.UserData;
ui.layerEnabled.Value=true; ui.layerEnabled.ValueChangedFcn(ui.layerEnabled,[]);
ui.layerTable.Data=[.001 1000 900 2 4 .002;.003 2600 800 5 5 0];
entries=[0 0 0 0 .8 300 .08 3];
for j=1:8, ui.wallFields(j).Value=entries(j); end
% 各层表数据必须独立保留；增加/删除层不能误删其他层的物性表。
table1=[200 800 1.8;400 1000 2.2;2000 1800 4];
ui.materialEnabled.Value=true; ui.materialTable.Data=table1;
ui.layerSelector.Value=2; ui.layerSelector.ValueChangedFcn(ui.layerSelector,[]);
assert(~ui.materialEnabled.Value && all(isnan(ui.materialTable.Data),'all'));
ui.addLayerButton.ButtonPushedFcn(ui.addLayerButton,[]);
assert(size(ui.layerTable.Data,1)==3);
ui.removeLayerButton.ButtonPushedFcn(ui.removeLayerButton,[]);
assert(size(ui.layerTable.Data,1)==2 && ui.layerTable.Data(end,6)==0);
ui.layerSelector.Value=1; ui.layerSelector.ValueChangedFcn(ui.layerSelector,[]);
assert(ui.materialEnabled.Value && isequal(ui.materialTable.Data,table1));
% 模式切换保留单层空表，同时恢复分层表。
ui.layerEnabled.Value=false; ui.layerEnabled.ValueChangedFcn(ui.layerEnabled,[]);
assert(~ui.materialEnabled.Value && all(isnan(ui.materialTable.Data),'all'));
ui.layerEnabled.Value=true; ui.layerEnabled.ValueChangedFcn(ui.layerEnabled,[]);
assert(ui.materialEnabled.Value && isequal(ui.materialTable.Data,table1));
ui.runWallButton.ButtonPushedFcn(ui.runWallButton,[]);
w=fig.UserData.wallThermal;
assert(nnz(w.validFace)>1000 && all(isnan(w.Biot),'all') && ~any(w.lumpedFace));
assert(numel(w.material.layers)==2 && w.material.contactResistance==.002);
assert(isequal(w.material.layers{1}.propertyTable.specificHeat,table1(:,2)));
assert(~isfield(w.material.layers{2},'propertyTable'));
assert(strcmp(ui.wallFields(1).Enable,'off'));
assert(fig.UserData.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(any(contains(string(ui.wallDescription.Value),'不用单一 Biot')));
ui.wallMetric.Value=6; ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData(w.oneDimensionalFace), ...
    w.oneDimensional.backTemperature(:,end)));
ui.wallDepthView.Value=true; ui.wallDepthView.ValueChangedFcn(ui.wallDepthView,[]);
lines=findobj(ui.wallHistoryAxis,'Type','line');
assert(numel(lines)==1 && nnz(diff(lines.XData)==0)==1);
assert(numel(findall(ui.wallHistoryAxis,'Type','constantline'))==1);
ui.wallTimeSlider.Value=1; ui.wallTimeSlider.ValueChangedFcn(ui.wallTimeSlider,[]);
lines=findobj(ui.wallHistoryAxis,'Type','line');
selected=find(find(w.oneDimensionalFace)==fig.UserData.wallSelectedFace,1);
assert(isequal(lines.YData,w.oneDimensional.finalDepthTemperature(selected,:)));
% 材料编辑没有污染压力共用状态，恢复已有计算时仍保留分层。
ui.selector.Value=2; ui.selector.ValueChangedFcn(ui.selector,[]);
assert(isequaln(ui.surface.FaceVertexCData,mesh.pressure.staticPressure/1000));
saved=evalin('base','mesh'); restored=pressureGUI(saved,cfg);
cleanup2=onCleanup(@()delete(restored)); %#ok<NASGU>
assert(restored.UserData.layerEnabled.Value);
assert(isequal(restored.UserData.materialTable.Data,table1));
assert(restored.UserData.layerTable.Data(1,6)==.002);
if nargin>0
    ui.tabs.SelectedTab=ui.materialTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'layered_material_gui.png'));
    ui.wallMetric.Value=1; ui.wallTimeSlider.Value=numel(w.time);
    ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
    ui.tabs.SelectedTab=ui.wallTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'layered_wall_profile_gui.png'));
end
ui.flightTable.Data=[0 30 5 0 0;.02 29 5 2 0;.04 30 5 0 0];
ui.runFlightButton.ButtonPushedFcn(ui.runFlightButton,[]);
wf=fig.UserData.wallThermal;
assert(isfield(wf,'coupling') && numel(wf.material.layers)==2);
assert(fig.UserData.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(size(wf.oneDimensional.backTemperature,2)==numel(wf.time));
ui.wallDepthView.Value=false; ui.wallMetric.Value=6;
ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(max(abs(ui.wallHistoryAxis.XLim-[wf.time(1) wf.time(end)]))<1e-12);
assert(diff(ui.wallAxis.CLim)>=.1-1e-10);
assert(numel(findall(ui.wallHistoryAxis,'Type','constantline'))==1);
if nargin>0
    ui.tabs.SelectedTab=ui.wallTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'layered_flight_back_temperature_gui.png'));
end
fprintf('Layered wall GUI passed: independent layer tables, restoration, fixed/flight coupling, energy, back map and final depth plot.\n');
end
