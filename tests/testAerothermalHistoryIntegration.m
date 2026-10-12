function testAerothermalHistoryIntegration(out)
%TESTAEROTHERMALHISTORYINTEGRATION 经典锥柱变工况热历史与内嵌 WebGL 联调。
% 模块一：运行真实内部外流/沿程热耦合；不将合成显示样本当成物理验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(tempdir,'RocketAeroHeat_history');end
if ~exist(out,'dir'),mkdir(out);end
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=0;cfg.flight.beta=0;cfg.flight.Twall=300;
mesh=buildCanonicalBody('coneCylinder',struct('resolution',24));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
mesh=calcFlightStationState(mesh,cfg.flight,cfg);
t=[0 .02 .06];
flight=struct('time',t,'Mach',[4 4 4],'altitude',[20000 19900 19800], ...
    'alpha',[0 5 10],'beta',[0 0 0]);
material=struct('density',2700,'specificHeat',900,'thickness',.003, ...
    'thermalConductivity',180,'emissivity',.8,'radiationTemperature',300);
h=calcFlightHeatingHistory(mesh,flight,cfg);
mesh=calcWallThermalHistory(mesh,h,material,cfg);mesh=calcAerothermalExposure(mesh);
w=mesh.wallThermal;assert(any(w.validFace) && isfield(w,'coupling'));
data=buildAerothermalViewerData(mesh);field=data.fields(strcmp({data.fields.key},'wallSolved'));
assert(isequal(field.history.time,w.time));
actual=reshape(field.history.values,size(w.temperature));
assert(isequaln(actual,w.temperature));
assert(contains(field.history.timeNote{end},'攻角 10'));
assert(mesh.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
% 外流显示与每个实际接受时间站的求解快照逐值一致，不能借用末站。
for key={'pressure','cp','edgePressure','edgeTemperature'}
    flowField=data.fields(strcmp({data.fields.key},key{1}));
    assert(isscalar(flowField) && isequal(flowField.history.time,w.time));
    for j=1:numel(w.time)
        meta=w.coupling.boundaryMetadata{j};s=meta.flowSnapshot;
        expected=s.fields(strcmp({s.fields.key},key{1}));
        ids=(j-1)*size(mesh.F,1)+(1:size(mesh.F,1));
        assert(isequaln(flowField.history.values(ids),expected.values));
        assert(isequal(flowField.history.status(ids),expected.status));
    end
end
p=data.fields(strcmp({data.fields.key},'pressure')).history;
values=reshape(p.values,size(mesh.F,1),[]);
assert(max(abs(values(:,end)-values(:,1)))>1);
assert(any(p.status==2),'给定背风或底压必须与计算值区分');
exportAerothermalWeb(mesh,fullfile(out,'cone_history'));save(fullfile(out,'cone_history.mat'),'mesh','cfg');
% 模块二：真实 uihtml 握手及像素回执；不是仅检查生成的 HTML 文本。
fig=uifigure('Visible','off','Position',[100 100 1200 850]);clean=onCleanup(@()delete(fig)); %#ok<NASGU>
tabs=uitabgroup(fig,'Position',[0 0 1200 850]);ui=createAerothermalWebTab(tabs,mesh);
started=tic;
while (isempty(ui.html.UserData.rendered) || ~ui.html.UserData.ready) && toc(started)<40
    drawnow;pause(.1);
end
gpu=ui.html.UserData.rendered;
assert(ui.html.UserData.ready && ~isempty(gpu) && gpu.webglError==0 && gpu.visiblePixelSamples>0);
assert(gpu.faces==size(mesh.F,1));
fprintf('History integration passed: %d faces, %d accepted stations, %.2f%% thermal area, energy residual %.3g J/m2; embedded GL error %d, visible samples %d.\n', ...
    size(mesh.F,1),numel(w.time),100*mesh.aerothermalExposure.summary.resolvedAreaFraction, ...
    mesh.aerothermalExposure.summary.maximumEnergyBalanceResidual,gpu.webglError,gpu.visiblePixelSamples);
end
