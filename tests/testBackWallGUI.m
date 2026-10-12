function testBackWallGUI(outputFolder)
%TESTBACKWALLGUI 在真实箭体上验证内侧输入、双面能量云图与固定/飞行耦合。
% 数据仅用于验证，内侧 h/温度不是程序推算的实际舱内状态。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg)); [raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
cfg.pressureGUI.visible='off'; fig=pressureGUI(mesh,cfg);
cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
ui=fig.UserData;
assert(~ui.backBoundaryEnabled.Value && strcmp(ui.backBoundaryFields(1).Enable,'off'));
ui.backBoundaryEnabled.Value=true;
ui.backBoundaryEnabled.ValueChangedFcn(ui.backBoundaryEnabled,[]);
assert(strcmp(ui.backBoundaryFields(1).Enable,'on'));
entries=[1000 1000 .001 10 .8 300 .1 3];
for j=1:8, ui.wallFields(j).Value=entries(j); end
failed=false; try, ui.runWallButton.ButtonPushedFcn(ui.runWallButton,[]); catch, failed=true; end
assert(failed && ~isfield(fig.UserData,'wallThermal'));
inner=[300 270 .8 250];
for j=1:4, ui.backBoundaryFields(j).Value=inner(j); end
ui.runWallButton.ButtonPushedFcn(ui.runWallButton,[]);
w=fig.UserData.wallThermal; e=fig.UserData.aerothermalExposure;
assert(w.backBoundary.active && ~any(w.lumpedFace));
assert(e.summary.maximumEnergyBalanceResidual<1e-3);
assert(e.summary.totalBackEnergyResolved>0 && all(w.backHeatFlux(w.validFace,1)>0));
assert(max(abs(w.netHeatFlux(w.validFace,:)-w.outerNetHeatFlux(w.validFace,:)+w.backHeatFlux(w.validFace,:)),[],'all')<1e-8);
assert(any(contains(string(ui.wallDescription.Value),'内侧累计净散热')));
for metric=[6 7 8 9 10 2]
    ui.wallMetric.Value=metric; ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
    switch metric
        case 6, expected=w.backTemperature(:,end);
        case 7, expected=w.backHeatFlux(:,end)/1e4;
        case 8, expected=e.cumulativeBackHeat(:,end)/1e6;
        case 9, expected=w.backConvectiveHeatFlux(:,end)/1e4;
        case 10, expected=w.backRadiativeHeatFlux(:,end)/1e4;
        case 2, expected=e.cumulativeNetHeat(:,end)/1e6;
    end
    assert(isequaln(ui.wallSurface.FaceVertexCData,expected));
    if metric>=7
        assert(abs(sum(ui.wallAxis.CLim))<1e-12);
        assert(any(contains(string(ui.wallHistoryAxis.Title.String),'背面 Newton')));
    end
end
% 已有结果恢复后仍采用原内侧条件；压力共享变量未受材料/边界编辑污染。
saved=evalin('base','mesh'); restored=pressureGUI(saved,cfg);
cleanup2=onCleanup(@()delete(restored)); %#ok<NASGU>
assert(restored.UserData.backBoundaryEnabled.Value);
assert(isequal(arrayfun(@(x)x.Value,restored.UserData.backBoundaryFields)',inner));
ui.selector.Value=2; ui.selector.ValueChangedFcn(ui.selector,[]);
assert(isequaln(ui.surface.FaceVertexCData,mesh.pressure.staticPressure/1000));
if nargin>0
    ui.tabs.SelectedTab=ui.backBoundaryTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'back_boundary_inputs_gui.png'));
    ui.tabs.SelectedTab=ui.wallTab; ui.wallMetric.Value=7;
    ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]); drawnow;
    exportapp(fig,fullfile(outputFolder,'back_heat_flux_gui.png'));
end
% 多层内侧换热进入实际变攻角飞行历程，继续用逐层储热核对。
ui.layerEnabled.Value=true; ui.layerEnabled.ValueChangedFcn(ui.layerEnabled,[]);
ui.layerTable.Data=[.0005 1000 1000 5 4 .001;.001 2000 800 8 4 0];
ui.flightTable.Data=[0 30 5 0 0;.02 29 5 2 0;.04 30 5 0 0];
ui.runFlightButton.ButtonPushedFcn(ui.runFlightButton,[]);
wf=fig.UserData.wallThermal; ef=fig.UserData.aerothermalExposure;
assert(wf.backBoundary.active && numel(wf.material.layers)==2 && isfield(wf,'coupling'));
assert(ef.summary.maximumEnergyBalanceResidual<1e-3);
assert(all(isfinite(wf.backHeatFlux(wf.validFace,:)),'all'));
assert(all(isnan(wf.backHeatFlux(~wf.validFace,:)),'all'));
ui.wallMetric.Value=8; ui.wallMetric.ValueChangedFcn(ui.wallMetric,[]);
assert(isequaln(ui.wallSurface.FaceVertexCData,ef.cumulativeBackHeat(:,end)/1e6));
if nargin>0
    drawnow; exportapp(fig,fullfile(outputFolder,'back_flight_energy_gui.png'));
end
fprintf('Back wall GUI passed: validation, restoration, flux/energy maps, signed scale and fixed/flight coupling.\n');
end
