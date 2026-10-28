function testBasePressureWorkbench(out)
%TESTBASEPRESSUREWORKBENCH 当前压力/气动力/网页和同窗口配置的底压集成校核。
% 模块用途：真实尖锥底面的关联压力进入积分和导出；不支持区域没有
% 给定Cp回退，整体载荷保持未解，底热/外缘温度不因底压接入而补值。
% 独立区域公式及薄层审核另由testMeshBasePressure校核；不是实验验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(root,'outputs','v54_validation');end
if ~exist(out,'dir'),mkdir(out);end
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.altitude=40000;cfg.flight.Mach=3;cfg.flight.alpha=0;cfg.flight.beta=0;
cfg.flight.Twall=300;cfg.heating.airThermalModel='caloricallyPerfect';
[mesh,meta]=buildCanonicalBody('cone',struct('radius',.2,'coneAngleDeg',9,'resolution',24));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);fs=calcFreestream(cfg.flight);
default=calcFlightStationState(mesh,cfg.flight,cfg);
default=integrateAerodynamicLoads(default,cfg);
base=strcmp(meta.faceRegion,'base');side=~base;

%% 1. 与报告区域结果一致，真实面积/法向独立积分底压力的载荷增量
selected=cfg;selected.basePressure=struct('model','sand93Laminar', ...
    'hasRearSupport',false,'hasGasInjection',false, ...
    'operationSource','经典孤立尖锥算例；工况设定无尾撑、无喷流');
current=calcFlightStationState(mesh,selected.flight,selected);
current=integrateAerodynamicLoads(current,selected);
p=current.pressure;r=p.baseRegion;
assert(r.pressureAvailable,r.status);
assert(all(p.fieldStatus(base)==1) && all(p.isComputedBase(base)));
assert(~any(p.isPrescribedBase(base)) && ~any(p.isPrescribedLeeward(base)));
assert(isequaln(p.staticPressure(side),default.pressure.staticPressure(side)));
assert(max(abs(p.staticPressure(base)-r.pressure(base)))<1e-12);
assert(max(abs(p.pressureIncrement-(p.staticPressure-fs.pressure)))<1e-10);
assert(current.aerodynamicLoads.available);
deltaP=p.staticPressure(base)-default.pressure.staticPressure(base);
expected=-sum((deltaP.*current.localGeometry.area(base)).*current.localGeometry.normal(base,:),1);
actual=current.aerodynamicLoads.forceBody-default.aerodynamicLoads.forceBody;
assert(norm(actual-expected)<1e-9);
assert(isequaln(current.heating.heatFlux,default.heating.heatFlux));
assert(all(isnan(current.heating.heatFlux(base))));
assert(all(isnan(current.boundaryLayerEdge.temperature(base))));
data=buildAerothermalViewerData(current);pressure=field(data,'pressure');cp=field(data,'cp');
for f={pressure,cp}
    e=f{1};assert(all(e.status(base)==1));
    assert(all(contains(e.methodNames(e.methodIndex(base)+1),'Eq.(8b)')));
    assert(all(contains(e.scopeNames(e.scopeIndex(base)+1),'无喷流/尾撑')));
end
heat=field(data,'heatFlux');edgeT=field(data,'edgeTemperature');
assert(all(heat.status(base)==0) && all(edgeT.status(base)==0));
defaultData=buildAerothermalViewerData(default);defaultP=field(defaultData,'pressure');
assert(abs(pressure.computedAreaFraction-defaultP.computedAreaFraction- ...
    sum(data.faceArea(base))/sum(data.faceArea))<1e-12);
exportAerothermalWeb(current,fullfile(out,'base_pressure'));

%% 2. 非零攻角与缺少操作声明时，底压缺失而非给定边界；其他区域仍计算
blockedCfg=selected;blockedCfg.flight.alpha=5;
blocked=calcFlightStationState(mesh,blockedCfg.flight,blockedCfg);
blocked=integrateAerodynamicLoads(blocked,blockedCfg);
assert(all(isnan(blocked.pressure.staticPressure(base))));
assert(all(blocked.pressure.fieldStatus(base)==0));
assert(~any(blocked.pressure.isPrescribedBase(base)) && ~blocked.aerodynamicLoads.available);
assert(all(isnan(blocked.aerodynamicLoads.forceBody)) && isnan(blocked.aerodynamicLoads.CD));
assert(abs(blocked.aerodynamicLoads.missingPressureArea-sum(data.faceArea(base)))<1e-12);
blockedData=buildAerothermalViewerData(blocked);blockedP=field(blockedData,'pressure');
assert(all(blockedP.status(base)==0) && all(isnan(blockedP.values(base))));
assert(all(contains(blockedP.methodNames(blockedP.methodIndex(base)+1),'零攻角')));
undeclared=selected;undeclared.basePressure=rmfield(undeclared.basePressure,'hasGasInjection');
u=modifiedNewtonianPressure(mesh,fs,undeclared);
assert(all(isnan(u.pressure.Cp(base))) && ~u.pressure.baseRegion.available);
bad=default;bad.pressure.pressureIncrement(find(side,1))=NaN;
mustFail(@()integrateAerodynamicLoads(bad,cfg),'RocketAeroHeat:Loads:MissingPressure');
placeholder=default;placeholder.pressure.fieldStatus(find(side,1))=0;
mustFail(@()integrateAerodynamicLoads(placeholder,cfg),'RocketAeroHeat:Loads:MissingPressure');
% 尾盖识别不完整时必须明确拒绝应用；不能让剩余底面偷偷保留给定Cp。
tilted=mesh;tilted.V(:,1)=tilted.V(:,1)+.01*tilted.V(:,2);
tilted=prepare(tilted,cfg);
mustFail(@()modifiedNewtonianPressure(tilted,fs,selected),'RocketAeroHeat:Pressure:BaseTarget');
warped=mesh;baseVertex=unique(mesh.F(base,:));
warped.V(baseVertex(1),1)=warped.V(baseVertex(1),1)+.001;
warped=prepare(warped,cfg);
mustFail(@()modifiedNewtonianPressure(warped,fs,selected),'RocketAeroHeat:Pressure:BaseTarget');

%% 3. 实际GUI在同一宿主显示设置，取消/应用与网页逐字段来源均一致
fig=pressureGUI(default,cfg);cleanup=onCleanup(@()closeFigure(fig)); %#ok<NASGU>
button=fig.UserData.basePressureButton;button.ButtonPushedFcn(button,[]);
editor=fig.UserData.basePressureEditor;
assert(ancestor(editor.page,'figure')==fig && strcmp(editor.page.Visible,'on'));
editor.model.Value='sand93Laminar';editor.thermalModel.Value='caloricallyPerfect';
editor.support.Value='no';editor.injection.Value='no';editor.source.Value='GUI经典孤立尖锥：无尾撑、无喷流';
editor.apply.ButtonPushedFcn(editor.apply,[]);drawnow;
shown=fig.UserData.currentMesh;
assert(strcmp(fig.UserData.currentCfg.basePressure.model,'sand93Laminar'));
assert(max(abs(shown.pressure.staticPressure(base)-p.staticPressure(base)))<1e-10);
assert(shown.aerodynamicLoads.available && strcmp(fig.UserData.basePressureButton.Visible,'on'));
ui=fig.UserData.webUI;start=tic;
while (~ui.html.UserData.ready || isempty(ui.html.UserData.rendered)) && toc(start)<40
    drawnow;pause(.1);
end
assert(ui.html.UserData.ready && ui.html.UserData.rendered.webglError==0);
assert(ui.html.UserData.rendered.visiblePixelSamples>0);
button=fig.UserData.basePressureButton;button.ButtonPushedFcn(button,[]);
editor=fig.UserData.basePressureEditor;editor.cancel.ButtonPushedFcn(editor.cancel,[]);drawnow;
assert(strcmp(fig.UserData.currentMesh.pressure.baseModel,'sand93Laminar'));
button=fig.UserData.basePressureButton;button.ButtonPushedFcn(button,[]);
editor=fig.UserData.basePressureEditor;editor.model.Value='prescribed';
editor.apply.ButtonPushedFcn(editor.apply,[]);drawnow;
assert(isequaln(fig.UserData.currentMesh.pressure.Cp,default.pressure.Cp));
assert(all(fig.UserData.currentMesh.pressure.fieldStatus(base)==2));
% 超域工况也必须实际进入工作台，不因载荷缺失而崩溃。
pressureGUI(blocked,blockedCfg,ancestor(fig.UserData.basePressureButton,'uipanel'));
assert(~fig.UserData.currentMesh.aerodynamicLoads.available);
save(fullfile(out,'base_pressure_workbench.mat'),'current','default','blocked','selected','cfg','meta');
fprintf('Base pressure workbench passed: %d actual base faces, computed pressure/loads/provenance; rejected base remains NaN, no heat fill; same-window GUI apply/cancel and WebGL.\n',nnz(base));
end
function f=field(data,key)
f=data.fields(strcmp({data.fields.key},key));
end
function mustFail(fun,id)
try,fun();catch ex,assert(strcmp(ex.identifier,id),ex.identifier);return;end
error('Expected error not raised: %s',id);
end
function closeFigure(fig)
if isgraphics(fig),delete(fig);end
end
function mesh=prepare(mesh,cfg)
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
end
