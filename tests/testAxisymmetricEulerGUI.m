function testAxisymmetricEulerGUI(outputFolder)
%TESTAXISYMMETRICEULERGUI 验证实际 STL 提取、内部流场和统一页绘图。
% 同时检查新场不会悄悄覆盖原有外流/气动热的算法来源。
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
raw=calcFaceGeometry(loadSTLMesh(cfg));[raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg);mesh=fixNormalDirection(mesh,cfg);');
[contour,geo]=extractAxisymmetricEulerContour(mesh);
assert(contour(end,2)==0 && geo.relativeRoundnessError<.02);
bad=mesh;bad.V(:,2)=2*bad.V(:,2);
try,extractAxisymmetricEulerContour(bad);error('Test:NotRejected','未拒绝非圆截面');
catch err,assert(strcmp(err.identifier,'RocketAeroHeat:Euler:Axisymmetry'));end
mesh=buildSurfaceCoordinates(mesh,cfg);mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
before=mesh.heating;cfg.pressureGUI.visible='off';fig=pressureGUI(mesh,cfg);
cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
ui=fig.UserData.eulerUI;fig.UserData.tabs.SelectedTab=ui.tab;
fig.UserData.tabs.SelectionChangedFcn(fig.UserData.tabs,[]);
assert(strcmp(fig.UserData.selector.Enable,'off'));
% 明确产生一次未收敛结果，再通过同网格续算完成；不继承旧状态标记。
ui.maxIterations.Value=1;
ui.runButton.ButtonPushedFcn(ui.runButton,[]);
assert(~fig.UserData.eulerResult.converged && ~fig.UserData.eulerResult.available);
ui.maxIterations.Value=50000;
ui.runButton.ButtonPushedFcn(ui.runButton,[]);
sol=fig.UserData.eulerResult;
fprintf('Actual Euler: available=%d residual=%g outflow=%d time=%g\n', ...
    sol.available,sol.residual,sol.boundary.supersonicOutflow,sol.elapsedSeconds);
assert(sol.converged && sol.boundary.clear && ~sol.edgeStateAvailable);
assert(sol.restart.used && ~sol.restart.sourceConverged);
assert(isfield(sol,'edgeCandidate') && ~sol.edgeCandidate.available);
assert(isequal(sol.edgeCandidate.wallIndex,flipud(find(sol.grid.bodyWall & sol.grid.wallR>0))));
assert(strcmp(sol.fluxScheme,'HLLenthalpy') && sol.diagnostics.maximumTotalEnthalpyDeviation<1e-6);
assert(contains(ui.description.Text,'上风总焓') && contains(ui.description.Text,'初场'));
% 该实际箭体尚未通过网格精度验证，界面必须公开总焓偏差。
% 不把一次粗网格恰好接近驻点关联值当作独立精度证明。
assert(~sol.diagnostics.gridAccuracyVerified);
assert(contains(ui.status.Text,'总焓最大偏差') && contains(ui.status.Text,'未接入热流'));
assert(all(isfinite(sol.diagnostics.totalEnthalpyRelativeError),'all'));
assert(isequaln(before,fig.UserData.heating));
assert(contains(ui.description.Text,'未') || contains(ui.guide.Text,'既有外流'));
for j=1:numel(ui.metric.Items)
    ui.metric.Value=ui.metric.Items{j};ui.metric.ValueChangedFcn(ui.metric,[]);drawnow;
    assert(~isempty(ui.axes.Children));
    switch j
        case 1
            plotted=findobj(ui.axes,'Type','patch');
            assert(isequal(plotted.FaceVertexCData,sol.pressure(:)/1000));
        case 2
            plotted=findobj(ui.axes,'Type','patch');
            assert(isequal(plotted.FaceVertexCData,sol.Mach(:)));
        case 3
            plotted=findobj(ui.axes,'Type','surface');
            assert(isequal(plotted.CData(:,1),sol.wall.pressure(sol.grid.bodyWall)/1000));
        case 4
            plotted=findobj(ui.axes,'Type','line');
            assert(isequal(plotted.YData(:),sol.wall.pressure(sol.grid.bodyWall)/1000));
        case 5
            plotted=findobj(ui.axes,'Type','line');
            assert(isequal(plotted.YData(:),max(sol.residualHistory(:,2),eps)));
        case 6
            plotted=findobj(ui.axes,'Type','surface');
            assert(isequal(plotted.CData(:,1),sol.edgeCandidate.signedSpeed));
        case 7
            plotted=findobj(ui.axes,'Type','surface');
            assert(isequal(plotted.CData(:,1),sol.edgeCandidate.temperature));
        case 8
            plotted=findobj(ui.axes,'Type','line','DisplayName','投影切向速度');
            assert(isequal(plotted.YData(:),sol.edgeCandidate.signedSpeed));
        case 9
            plotted=findobj(ui.axes,'Type','line');
            assert(isequaln(plotted.YData(:),sol.edgeCandidate.velocityGradient));
        case 10
            plotted=findobj(ui.axes,'Type','line');
            assert(isequal(plotted.YData(:),sol.edgeCandidate.diagnostics.totalPressureRatio));
        case 11
            plotted=findobj(ui.axes,'Type','line','DisplayName','第一层单元压力');
            assert(isequal(plotted.YData(:),sol.edgeCandidate.pressure/1000));
    end
    if nargin>0
        if ~isfolder(outputFolder),mkdir(outputFolder);end
        exportapp(fig,fullfile(outputFolder,sprintf('euler_view_%d.png',j)));
    end
end
ui.metric.Value=ui.metric.Items{1};ui.metric.ValueChangedFcn(ui.metric,[]);
assert(strcmp(ui.axes.ColorScale,'log'));
ui.fullField.Value=true;ui.fullField.ValueChangedFcn(ui.fullField,[]);
assert(abs(ui.axes.YLim(2)-sol.grid.outerRadius)<1e-10);
ui.logPressure.Value=false;ui.logPressure.ValueChangedFcn(ui.logPressure,[]);
assert(strcmp(ui.axes.ColorScale,'linear'));
% 实际箭体选择鼻部法向过渡：无需换窗口，也不能覆盖旧工程热流。
ui.gridType.Value='noseBlend';ui.runButton.ButtonPushedFcn(ui.runButton,[]);
blended=fig.UserData.eulerResult;
assert(blended.available && strcmp(blended.boundary.topology,'axisAtNose'));
assert(contains(ui.description.Text,'过渡网格') && ~blended.edgeStateAvailable);
assert(isequaln(before,fig.UserData.heating));
ui.fullField.Value=false;ui.logPressure.Value=true;
for j=[1 6]
    ui.metric.Value=ui.metric.Items{j};ui.metric.ValueChangedFcn(ui.metric,[]);drawnow;
    if j==1
        plotted=findobj(ui.axes,'Type','patch');
        assert(isequal(plotted.FaceVertexCData,blended.pressure(:)/1000));
    else
        plotted=findobj(ui.axes,'Type','surface');
        assert(isequal(plotted.CData(:,1),blended.edgeCandidate.signedSpeed));
    end
    if nargin>0,exportapp(fig,fullfile(outputFolder,sprintf('euler_blended_%d.png',j)));end
end
fig.UserData.tabs.SelectedTab=fig.UserData.mapTab;
fig.UserData.tabs.SelectionChangedFcn(fig.UserData.tabs,[]);
assert(strcmp(fig.UserData.selector.Enable,'on'));
oblique=sol.freestream;oblique.alpha=5;
offUI=createAxisymmetricEulerTab(fig.UserData.tabs,mesh,oblique);
assert(strcmp(offUI.runButton.Enable,'off'));
% 原始 STL 轮廓和当前计算场供可复现诊断，不作为真实精度校核。
if nargin>0,save(fullfile(outputFolder,'euler_actual.mat'),'sol','blended','contour','geo');end
fprintf('PASS actual model Euler and unified GUI\n');
end
