function check_ui(out,job,archived)
% 模块一：按普通 STL 流程导入反序文件，在新窗口从真实半步任务打开曲线。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(out,'input.stl');c.scale.Value=1;press(c.load);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'shortSupport65');press(fig.UserData.controls.generate);
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.shortBranch,'half');press(fig.UserData.controls.confirm);
c=fig.UserData.controls;
if archived,press(c.readShort);else,c.shortDirectory.Value=job;press(c.reopenShort);end
state=fig.UserData.workflowState;assert(state.currentStep==5 && strcmp(state.status,'solved'),'%s',fig.UserData.status.Text);
m=state.result;ui=fig.UserData.webUI;c=fig.UserData.controls;press(c.openAzimuth);
tab=ui.tab.Parent.SelectedTab;assert(strcmp(tab.Tag,'Short65Azimuth'));html=tab.UserData.html;
rear=waitView(tab,[119.8;166.7],'nu',0);assert(all(rear.gaps==0) && numel(rear.values)==2);
t=m.shortSupportAzimuth.points;
if archived,assert(isequal(t.fileFaceID,t.referenceFaceID));else,assert(all(t.fileFaceID==32481-t.referenceFaceID));end
assert(isequal(rear.faces(:),t.fileFaceID(ismember(t.angleDeg,[119.8 166.7]) & t.phiDeg==0)));
% 模块二：使用同一组件的方位与变量选择，确认管遮盖断线、两个背风值和极点读数。
html.Data=struct('kind','azimuth-view','angle',90,'field','flux','phi',90);
gap=waitView(tab,90,'flux',90);assert(gap.gaps==2);
html.Data=struct('kind','azimuth-view','angle','rear','field','nu','phi',280);
pair=waitView(tab,[119.8;166.7],'nu',280);assert(pair.values(1)>pair.values(2));
html.Data=struct('kind','azimuth-view','angle',0,'field','nu','phi',0);
pole=waitView(tab,0,'nu',0);assert(pole.rows==1);
% 将曲线反馈的面号送入同一个探针回调，核对反序 STL 的实际面号。
html.Data=struct('kind','azimuth-locate','frame',rear.frame,'face',pair.faces(1));html.DataChangedFcn(html,[]);drawnow;
probe=ui.evidenceView.getState();assert(ui.tab.Parent.SelectedTab==ui.evidenceTab && probe.selectedFaceID==pair.faces(1));
exportAerothermalWeb(m,fullfile(out,'workbench'));
% 模块三：验证归档与任务的源绑定、拒绝受损包，并读整步确认没有半步曲线。
sf=jsondecode(fileread(fullfile(job,'export','half','summary.json')));
cases=fullfile(fileparts(which('loadShortSupportSurface')),'cases','short_sphere');
arch=readShortSupportAzimuth(cases,'',sf,fig.UserData.caseReport.stations,out);
assert(isequaln(arch.points.nu,t.nu));
bad=fullfile(out,'bad');mkdir(fullfile(bad,'export'));copyfile(fullfile(job,'completed.json'),bad);
f=fopen(fullfile(bad,'export','stationAzimuth.zip'),'w');fprintf(f,'x');fclose(f);
rejected=false;
try,readShortSupportAzimuth(cases,bad,sf,fig.UserData.caseReport.stations,out);catch ex,rejected=strcmp(ex.identifier,'RAH:Short:Azimuth');end
assert(rejected);
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.shortBranch,'full');press(fig.UserData.controls.confirm);
press(fig.UserData.controls.readShort);assert(strcmp(fig.UserData.workflowState.status,'solved'),'%s',fig.UserData.status.Text);
assert(~isfield(fig.UserData.workflowState.result,'shortSupportAzimuth') && isempty(findobj(fig.UserData.webUI.tab.Parent,'Tag','Short65Azimuth')));
proof=struct('archived',archived,'rear',rear,'gap',gap,'otherAzimuth',pair,'pole',pole,'rows',height(t),'selectedFileFace',probe.selectedFaceID, ...
    'archiveMatches',true,'fullStepIsolated',true,'damagedRejected',rejected);
f=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(f>=0);close=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function msg=waitView(tab,angles,field,phi)
timer=tic;
while toc(timer)<35
    drawnow;msg=tab.UserData.rendered;
    if isstruct(msg) && isequal(msg.angles(:),angles(:)) && strcmp(msg.field,field) && abs(msg.phi-phi)<1e-6,return;end
    pause(.1);
end
error('未收到正确的周向曲线渲染回执。');
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,val)
c.Value=val;c.ValueChangedFcn(c,[]);drawnow;
end
