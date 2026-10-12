function testPressure(outputFolder)
% 压力模块回归：解析极限、原始入射量、背风近似、真实模型和 GUI 回调。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg = createCaseConfig();'); cfg.output.verbose = false;
fs = calcFreestream(cfg.flight); fs.Mach = 5; fs.gamma = 1.4;
[cp,ref] = calcCpMax(fs);
assert(abs(cp-1.808769960702756) < 1e-10);
assert(abs(ref.postShockPressureRatio-29) < 1e-12);
assert(abs(ref.stagnationPressure-(fs.pressure+ref.dynamicPressure*cp)) < 1e-9*ref.stagnationPressure);
sonic = fs; sonic.Mach = 1;
[~,s] = calcCpMax(sonic);
assert(abs(s.postShockMach-1) < 1e-12);
assert(abs(s.stagnationPressureRatio-1.2^3.5) < 1e-12);
fast = fs; fast.Mach = 1000;
assert(abs(calcCpMax(fast)-1.839371) < 1e-5);
bad = fs; bad.Mach = 0.9;
mustFail(@()calcCpMax(bad));

%% 解析面元：法向分别沿 +X、-X 和略微迎风方向
normals = [1 0 0;-1 0 0;0.001 sqrt(1-0.001^2) 0];
m.V = zeros(9,3); m.F = reshape(1:9,3,[])';
for i=1:3
    n = normals(i,:); a = [0 0 1]; b = cross(n,a);
    m.V(m.F(i,:),:) = [0 0 0;a;b]+[i*3 0 0];
end
g.normal = normals; g.center = zeros(3,3); g.area = 0.5*ones(3,1);
for i=1:3, g.center(i,:) = mean(m.V(m.F(i,:),:),1); end
g.region = repmat({'sidewall'},3,1);
for name = {'isAftCap','isNoseCap','isAftTipAdjacent','isNoseTipAdjacent','poorAngle','featureFace','topologyIssue'}
    g.(name{1}) = false(3,1);
end
g.isAftCap(2) = true; m.localGeometry = g;
fs.unitFlowDirection = [-1 0 0];
m = modifiedNewtonianPressure(m,fs,cfg);
assert(max(abs(m.pressure.Cp-[cp;0;cp*1e-6])) < 1e-12);
assert(m.aerothermalRegion.isGrazing(3) && m.pressure.Cp(3)>0);
custom = cfg; custom.pressure.aftBaseCp = -0.01;
mc = modifiedNewtonianPressure(m,fs,custom);
assert(mc.pressure.Cp(2) == -0.01);
custom.pressure.aftBaseCp = -100;
mustFail(@()modifiedNewtonianPressure(m,fs,custom));

%% 默认 STL 全链路与 GUI：切换字段后核对实际着色数组
evalc('raw = loadSTLMesh(cfg); raw = calcFaceGeometry(raw);');
[raw,~] = alignGeometryAxis(raw,'pca');
evalc('mesh = autoRemeshSurface(raw,cfg); mesh = fixNormalDirection(mesh,cfg);');
evalc('mesh = buildSurfaceCoordinates(mesh,cfg);');
mesh = analyzeLocalGeometry(mesh,cfg);
fs = calcFreestream(cfg.flight);
before = mesh; mesh = modifiedNewtonianPressure(mesh,fs,cfg);
assert(isequal(mesh.V,before.V) && isequal(mesh.F,before.F));
assert(all(mesh.pressure.Cp >= 0 & mesh.pressure.Cp <= mesh.pressure.CpMax));
cfg.pressureGUI.visible = 'off';
fig = pressureGUI(mesh,cfg); cleanup = onCleanup(@()delete(fig)); %#ok<NASGU>
ui = fig.UserData;
% 径向放大只改变三维显示比例，四个页面保持一致，截面取样仍按原坐标。
assert(ui.displayScale.Value==3);
originalVertices=ui.surface.Vertices;
originalPressure=ui.surface.FaceVertexCData;
originalHeat=ui.heatSurface.FaceVertexCData;
originalProfile=ui.selectedProfile;
for scale=[1 5 3]
    ui.displayScale.Value=scale;
    ui.displayScale.ValueChangedFcn(ui.displayScale,[]);
    axesList=[ui.axes ui.sectionAxes ui.heatAxis ui.wallAxis];
    for axisIndex=1:numel(axesList)
        assert(isequal(axesList(axisIndex).DataAspectRatio,[scale 1 1]));
    end
    assert(isequal(ui.surface.Vertices,originalVertices));
    assert(isequal(ui.surface.FaceVertexCData,originalPressure));
    assert(isequaln(ui.heatSurface.FaceVertexCData,originalHeat));
    assert(all(isfinite([ui.axes.CameraViewAngle ui.heatAxis.CameraViewAngle])));
end
ui.sectionSurface.ButtonDownFcn(ui.sectionSurface, ...
    struct('IntersectionPoint',[originalProfile.stationX 0 0]));
selected=fig.UserData.selectedProfile;
assert(selected.stationX==originalProfile.stationX);
assert(isequaln(selected.Cp,originalProfile.Cp));
% 修改参考面积只改变系数，不改变已计算压力与实际气动力。
originalTable=ui.loadsTable.Data;
ui.referenceFields(1).Value=2*ui.referenceFields(1).Value;
ui.applyReference.ButtonPushedFcn(ui.applyReference,[]);
newTable=ui.loadsTable.Data;
assert(abs(newTable{13,2}-originalTable{13,2}/2)<1e-12);
assert(isequal(newTable(1:9,2),originalTable(1:9,2)));
ui.autoReference.ButtonPushedFcn(ui.autoReference,[]);
assert(isequal(ui.loadsTable.Data,originalTable));
% 固定侧视角按钮应分别显示上下两侧，攻角数值必须出现在标题中。
ui.viewSelector.Value = 3; feval(ui.viewSelector.ValueChangedFcn,ui.viewSelector,[]);
[~,el] = view(ui.axes); assert(el < 0);
ui.viewSelector.Value = 2; feval(ui.viewSelector.ValueChangedFcn,ui.viewSelector,[]);
[~,el] = view(ui.axes); assert(el > 0);
ui.viewSelector.Value = 1; feval(ui.viewSelector.ValueChangedFcn,ui.viewSelector,[]);
for choice = 1:3
    ui.selector.Value = choice;
    feval(ui.selector.ValueChangedFcn,ui.selector,[]);
    expected = {mesh.pressure.Cp,mesh.pressure.staticPressure/1000,mesh.pressure.pressureIncrement/1000};
    assert(isequal(ui.surface.FaceVertexCData,expected{choice}));
end
ui.selector.Value=5; feval(ui.selector.ValueChangedFcn,ui.selector,[]);
assert(isequaln(ui.surface.FaceVertexCData,ui.boundaryLayerEdge.pressure/1000));
ui.wire.Value = 1; feval(ui.wire.ValueChangedFcn,ui.wire,[]);
assert(isnumeric(ui.surface.EdgeColor));
ui.wire.Value = 0; feval(ui.wire.ValueChangedFcn,ui.wire,[]);
ui.selector.Value = 1; feval(ui.selector.ValueChangedFcn,ui.selector,[]); drawnow;
if nargin > 0
    exportgraphics(ui.axes,fullfile(outputFolder,'pressure_preview.png'),'Resolution',120);
    % 完整窗口截图包含控件，供布局检查。
    exportapp(fig,fullfile(outputFolder,'pressure_gui.png'));
end
fprintf('Pressure and GUI tests passed: CpMax=%.8f, p=[%.3f, %.3f] Pa, faces=%d.\n', ...
    mesh.pressure.CpMax,min(mesh.pressure.staticPressure),max(mesh.pressure.staticPressure),size(mesh.F,1));
% 正攻角默认观察下侧迎风面，且圆柱压力随攻角增大。
side = abs(mesh.localGeometry.normal(:,1)) < 0.01 & mesh.localGeometry.isSidewall;
lastMaximum = -Inf;
for angle = [0 10 20 30]
    flight = cfg.flight; flight.alpha = angle; flight.beta = 0;
    tilted = modifiedNewtonianPressure(mesh,calcFreestream(flight),cfg);
    maximum = max(tilted.pressure.Cp(side));
    assert(maximum > lastMaximum); lastMaximum = maximum;
end
tiltedFigure = pressureGUI(tilted,cfg);
cleanupTilted = onCleanup(@()delete(tiltedFigure)); %#ok<NASGU>
[~,el] = view(tiltedFigure.UserData.axes); assert(el < 0);
if nargin > 0
    exportapp(tiltedFigure,fullfile(outputFolder,'pressure_alpha30_gui.png'));
    tiltedUI=tiltedFigure.UserData;
    tiltedUI.selector.Value=5;
    feval(tiltedUI.selector.ValueChangedFcn,tiltedUI.selector,[]);
    exportapp(tiltedFigure,fullfile(outputFolder,'external_alpha30_gui.png'));
end
end

function mustFail(f)
% 非法物理输入必须报错，不得产生可用外观的结果。
try, f(); catch, return; end
error('预期错误未发生。');
end
