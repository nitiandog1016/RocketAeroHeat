function testEulerLinearGUI(outputFolder)
%TESTEULERLINEARGUI 独立解析圆锥在统一页计算，核对算法与取样来源。
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);
c=[linspace(-1,0,9)',-linspace(-1,0,9)'*tand(10)];
theta=(0:31)'*2*pi/32;mesh.V=zeros(0,3);ids=cell(size(c,1),1);
for k=1:size(c,1)
 start=size(mesh.V,1);
 if c(k,2)==0,points=[c(k,1),0,0];
 else,points=[repmat(c(k,1),32,1),c(k,2)*cos(theta),c(k,2)*sin(theta)];end
 mesh.V=[mesh.V;points];ids{k}=(start+1:size(mesh.V,1))'; %#ok<AGROW>
end
mesh.param.sectionX=c(:,1);mesh.param.sectionNodeIDs=ids;
fig=uifigure('Visible','off','Position',[30 30 1580 940]);cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
fig.UserData=struct;tabs=uitabgroup(fig,'Position',[0 0 1580 940]);
ui=createAxisymmetricEulerTab(tabs,mesh,fs);assert(ui.spatialOrder.Value==1);
ui.axialCells.Value=48;ui.radialCells.Value=36;ui.spatialOrder.Value=2;
ui.fluxScheme.Value='HLL';ui.maxIterations.Value=16000;
ui.runButton.ButtonPushedFcn(ui.runButton,[]);sol=fig.UserData.eulerResult;
assert(sol.available && sol.spatialOrder==2 && sol.edgeCandidate.sample.fromReconstructedTrace);
assert(contains(ui.description.Text,'Venkatakrishnan') && ~sol.edgeStateAvailable);
ui.spatialOrder.Value=1; % 尚未重算：方法标签必须仍是既有解的来源。
for index=[3 6 8 11]
 ui.metric.Value=ui.metric.Items{index};ui.metric.ValueChangedFcn(ui.metric,[]);
 assert(contains(ui.description.Text,'重构'));
 if index==8
  lines=findobj(ui.axes,'Type','line');assert(any(strcmp({lines.DisplayName},'重构壁面迹法向速度')));
 elseif index==11
  lines=findobj(ui.axes,'Type','line');assert(any(strcmp({lines.DisplayName},'重构壁面迹压力')));
 end
end
ui.spatialOrder.Value=2;ui.metric.Value=ui.metric.Items{3};ui.metric.ValueChangedFcn(ui.metric,[]);
if nargin>0
 if ~isfolder(outputFolder),mkdir(outputFolder);end
 exportapp(fig,fullfile(outputFolder,'linear_euler_cone_pressure.png'));
end
%% 镜像与总焓组合必须使用真实回调，图中文字记录实际算法
ui.reconstructionMode.Value=4;ui.fluxScheme.Value='HLLenthalpy';
ui.runButton.ButtonPushedFcn(ui.runButton,[]);thermo=fig.UserData.eulerResult;
assert(thermo.available && strcmp(thermo.boundaryReconstruction,'slipGhost'));
assert(strcmp(thermo.reconstructionVariables,'enthalpy') && thermo.diagnostics.maximumTotalEnthalpyDeviation<1e-6);
assert(contains(ui.description.Text,'镜像') && contains(ui.description.Text,'rho/u/v/H'));
assert(thermo.restart.used && strcmp(thermo.restart.sourceReconstructionVariables,'primitive'));
ui.metric.Value=ui.metric.Items{7};ui.metric.ValueChangedFcn(ui.metric,[]);
plotted=findobj(ui.axes,'Type','surface');assert(isequal(plotted.CData(:,1),thermo.edgeCandidate.temperature));
ui.reconstructionMode.Value=1;ui.metric.ValueChangedFcn(ui.metric,[]);
assert(contains(ui.description.Text,'rho/u/v/H')); % 不把尚未应用的选择写成结果来源。
ui.reconstructionMode.Value=4;
if nargin>0,exportapp(fig,fullfile(outputFolder,'thermodynamic_euler_cone_temperature.png'));end
fprintf('Limited linear GUI passed: unified selector, actual solver, provenance and trace labels.\n');
end
