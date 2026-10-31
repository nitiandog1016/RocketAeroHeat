function report=reviewSU2ThermalView(out,capture)
%REVIEWSU2THERMALVIEW 复用唯一经典数值结果核对最终温度文案、绝热字段及旧等温兼容。
% 只打开同一生产结果组件的临时测试窗口，不重算网格/流场，不建立永久研究页。
if nargin<2,capture=false;end
assert(~isfolder(out));mkdir(out);d=load('outputs/v114_validation/actual01/classic_result.mat','r');r=d.r;
fig=uifigure('Name','RocketAeroHeat · 热边界结果复核','Position',[20 30 1460 850]);guard=onCleanup(@()close(fig));
host=uipanel(fig,'Units','normalized','Position',[0 0 1 1]);ui=createSU2ResearchView(host,r);
ui.select.Value='heatFluxRaw';ui.select.ValueChangedFcn(ui.select,[]);ui.setSectionX(.158496);waitRender(ui);
ui.exportTo(fullfile(out,'fields.csv'));ui.exportSectionTo(fullfile(out,'section.csv'));ui.inspectFace(1454,[]);
assert(contains(strjoin(ui.probe.Value,newline),'未指定（绝热壁）') && ~contains(strjoin(ui.probe.Value,newline),'NaN'));
ui.select.Value='pressureRaw';ui.select.ValueChangedFcn(ui.select,[]);waitRender(ui);
assert(isequal(ui.geometry.html.UserData.payload.faceField.values,r.pressureRaw));
report=struct('productionComponentReused',true,'meshOrSolverRerun',false,'givenTwExplicitlyUnspecified',true, ...
    'rendered',ui.geometry.html.UserData.rendered,'physicsQualified',false,'newQualifiedArea',0);
write(fullfile(out,'report.json'),report);
if capture
    write(fullfile(out,'ready_for_capture.json'),struct('name',fig.Name));fprintf('V114_REVIEW_READY\n');
    t=tic;while ~isfile(fullfile(out,'capture_done.json')) && toc(t)<60,drawnow;pause(.1);end
end
% 原v107等温快照没有thermalBoundary字段，给定Tw仍为300K，原面/场不变。
old=load('docs/validation/v107/replay_seed/mapped_research.mat','result');delete(host.Children);
ui=createSU2ResearchView(host,old.result);ui.select.Value='wallInput';ui.select.ValueChangedFcn(ui.select,[]);waitRender(ui);
assert(all(ui.geometry.html.UserData.payload.faceField.values==300) && ui.geometry.html.UserData.payload.sectionView.valueStatus==2);
ui.exportTo(fullfile(out,'legacy_isothermal.csv'));
write(fullfile(out,'checks.json'),struct('legacyScalarWallTemperature',true,'givenAdiabaticHeatFlux',true,'unknownAdiabaticTw',true));
fprintf('Final production thermal view and legacy scalar isothermal snapshot passed without solver/mesher.\n');
end
function waitRender(ui)
t=tic;while toc(t)<40
    drawnow;s=ui.geometry.html.UserData;
    if ~isempty(s.rendered)&&s.rendered.revision==s.revision&&~isempty(ui.getSection()),break;end
    pause(.1);
end
assert(~isempty(s.rendered)&&s.rendered.revision==s.revision&&s.rendered.webglError==0&&~isempty(ui.getSection()));
end
function write(file,data)
fid=fopen(file,'w','n','UTF-8');assert(fid>=0);guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(data));
end
