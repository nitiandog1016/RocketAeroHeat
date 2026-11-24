function check_ui(out,job)
% 模块一：在真实五步工作台导入米制 STL，重开本次从原生帧复读的 14D 两半步结果。
fig=rocketAeroHeatWorkbench();closeFig=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(out,'supported_m.stl');c.scale.Value=1;press(c.load);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'shortSupport65');
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
choose(fig.UserData.controls.shortBranch,'half');press(fig.UserData.controls.confirm);
c=fig.UserData.controls;c.shortDirectory.Value=job;press(c.reopenShort);
state=fig.UserData.workflowState;assert(state.currentStep==5 && strcmp(state.status,'solved'),'%s',fig.UserData.status.Text);
m=state.result;t=fig.UserData.caseReport.stations;c=fig.UserData.controls;
assert(height(t)==7 && all(t.gradientStatus=="已重建") && isequal(t,m.shortSupportStationTable));
assert(all(t.gradientFrameSHA256==string(m.evidenceIdentity.sourceFrame.sha256)));
% 模块二：实际切换三类列，核对守恒值、梯度值、面号及取样信息。
assert(isequal(c.stations.Data,t{:,{'angleDeg','fileFaceID','reportedNu','computedNu'}}));
choose(c.stationMethod,'nu');assert(isequal(c.stations.Data,t{:,{'angleDeg','reportedNu','computedNu','gradientNu'}}));
choose(c.stationMethod,'flux');assert(isequal(c.stations.Data,t{:,{'angleDeg','fluxInWm2','gradientInWm2'}}));
choose(c.stationMethod,'profile');assert(isequal(c.stations.Data,t{:,{'angleDeg','layerCount','firstDistanceM','linearNu','wideNu'}}));
choose(c.stationMethod,'main');
assert(~isempty(c.exportStations.ButtonPushedFcn));
ui=fig.UserData.webUI;timer=tic;
while isempty(ui.html.UserData.rendered) && toc(timer)<50,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.webglError==0 && gl.faces==32480 && gl.visiblePixelSamples>0);
s=jsondecode(fileread(fullfile(job,'export','half','summary.json')));s.outW=s.energy.wallW;
v=addShortSupportSurfaceTab(ui.evidenceTab.Parent,t,s);assert(numel(v.line.Parent.Children)==3);
v.exportTo(fullfile(out,'tab.csv'));exportAerothermalWeb(m,fullfile(out,'workbench'));
% 模块三：单独核对旧归档、稳态、20D 工况和受损记录的来源边界。
cases=fullfile(fileparts(which('loadShortSupportSurface')),'cases','short_sphere');
for name={'supportedSurface65.zip','supportedFine65.zip','supportedAngularSurface65.zip','supportedTime65Full.zip','supportedTime65Half.zip'}
    [station,summary]=fromArchive(cases,name{1},out);
    a=attachShortStationGradients(station,summary,cases,'');
    assert(height(a)==7 && all(isfinite(a.gradientNu)) && all(a.gradientFrameSHA256==string(summary.frameSHA256)));
end
[station,other]=fromArchive(fullfile(cases,'outer20'),'supportedTime65Half.zip',out);
assert(isequal(attachShortStationGradients(station,other,cases,''),station));
bad=fullfile(out,'damaged');mkdir(fullfile(bad,'export'));copyfile(fullfile(job,'completed.json'),bad);
file=fopen(fullfile(bad,'export','stationGradients.json'),'w');fprintf(file,'{}');fclose(file);
rejected=false;
try,attachShortStationGradients(t,s,cases,bad);catch ex,rejected=strcmp(ex.identifier,'RAH:Short:Gradient');end
assert(rejected);
proof=struct('frame',m.evidenceIdentity.sourceFrame.sha256,'faces',gl.faces,'visiblePixels',gl.visiblePixelSamples, ...
    'stationsRead',height(t),'referenceFaceIDs',t.referenceFaceID,'fileFaceIDs',t.fileFaceID,'gradientNu',t.gradientNu, ...
    'allMainNu',t.computedNu,'status',t.gradientStatus,'damagedRejected',rejected,'archiveSets',5);
f=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(f>=0);guard=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function [t,s]=fromArchive(folder,name,out)
% 测站原表和帧摘要取自随版本交付的表面结果包。
fixture=fullfile(out,'station-fixture');unzip(fullfile(folder,name),fixture);
t=readtable(fullfile(fixture,'stations.csv'));s=jsondecode(fileread(fullfile(fixture,'summary.json')));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,val)
c.Value=val;c.ValueChangedFcn(c,[]);drawnow;
end
