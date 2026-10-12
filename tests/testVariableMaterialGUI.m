function testVariableMaterialGUI(outputFolder)
%TESTVARIABLEMATERIALGUI 实际箭体的材料表输入、热耦合、恢复和统一窗口。
% 测试表为人工数据，只验证软件路径，不代表真实材料性质。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg)); [raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
cfg.pressureGUI.visible='off';
fig=pressureGUI(mesh,cfg); cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
ui=fig.UserData;
assert(~ui.materialEnabled.Value && all(isnan(ui.materialTable.Data),'all'));
ui.addMaterialButton.ButtonPushedFcn(ui.addMaterialButton,[]);
assert(size(ui.materialTable.Data,1)==3);
ui.removeMaterialButton.ButtonPushedFcn(ui.removeMaterialButton,[]);
assert(size(ui.materialTable.Data,1)==2);
% 启用表后，壁温页标量 cp/k 可以留空值零；实际值必须来自表。
entries=[1000 0 .01 0 .8 300 .04 3];
for j=1:8, ui.wallFields(j).Value=entries(j); end
ui.materialEnabled.Value=true;
failed=false;
try, ui.runWallButton.ButtonPushedFcn(ui.runWallButton,[]); catch, failed=true; end
assert(failed && ~isfield(fig.UserData,'wallThermal'));
values=[200 900 10;400 1100 12;2000 1900 22];
ui.materialTable.Data=values;
ui.runWallButton.ButtonPushedFcn(ui.runWallButton,[]);
w=fig.UserData.wallThermal;
assert(nnz(w.validFace)>1000 && ~any(w.lumpedFace));
assert(w.material.specificHeat==1000 && w.material.thermalConductivity==11);
assert(strcmp(ui.wallFields(2).Enable,'off') && strcmp(ui.wallFields(4).Enable,'off'));
assert(ui.wallFields(2).Value==1000 && ui.wallFields(4).Value==11);
assert(contains(w.materialMethod,'分段线性'));
assert(any(contains(string(ui.wallDescription.Value),'材料物性')));
assert(fig.UserData.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
% 材料输入不能覆盖压力页的共享状态；重新切换压力显示应保持原数据。
ui.selector.Value=2; ui.selector.ValueChangedFcn(ui.selector,[]);
assert(isequaln(ui.surface.FaceVertexCData,mesh.pressure.staticPressure/1000));
ui.selector.Value=1; ui.selector.ValueChangedFcn(ui.selector,[]);
saved=evalin('base','mesh');
restored=pressureGUI(saved,cfg); cleanup2=onCleanup(@()delete(restored)); %#ok<NASGU>
assert(restored.UserData.materialEnabled.Value);
assert(isequal(restored.UserData.materialTable.Data,values));
assert(restored.UserData.wallFields(1).Value==1000);
if nargin>0
    ui.tabs.SelectedTab=ui.materialTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'material_property_table_gui.png'));
    ui.tabs.SelectedTab=ui.wallTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'variable_material_wall_gui.png'));
end
% 同一材料表进入变工况耦合，仍按比热积分核对材料储热。
ui.flightTable.Data=[0 30 5 0 0;.02 29 5 2 0;.04 30 5 0 0];
ui.runFlightButton.ButtonPushedFcn(ui.runFlightButton,[]);
wf=fig.UserData.wallThermal;
assert(contains(wf.materialMethod,'分段线性') && isfield(wf,'coupling'));
assert(fig.UserData.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(isequal(wf.material.propertyTable.specificHeat,values(:,2)));
if nargin>0
    ui.tabs.SelectedTab=ui.wallTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'variable_material_flight_gui.png'));
end
fprintf('Variable material GUI passed: table validation, fixed/flight coupling, saved inputs and energy.\n');
end
