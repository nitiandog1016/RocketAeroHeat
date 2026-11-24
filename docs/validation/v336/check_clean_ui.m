function check_clean_ui(out,mode,job)
% 模块一：在独立安装的真实工作台加载米制 STL，使用主工况 14D 和两半步。
if nargin<3,job='';end
fig=rocketAeroHeatWorkbench();clean=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(fileparts(out),'supported_m.stl');c.scale.Value=1;press(c.load);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'shortSupport65');
assert(fig.UserData.controls.shortRadius.Value==14);press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;assert(~isempty(s.mesh),'%s',fig.UserData.status.Text);
assert(isequal(s.mesh.shortSupportImport.referenceFaceID,(1:32480)'));
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.shortBranch,'half');press(fig.UserData.controls.confirm);
c=fig.UserData.controls;assert(fig.UserData.workflowState.currentStep==4);
bad={};
if strcmp(mode,'solve')
    fprintf('14D new calculation selected in clean installation\n');press(c.run);
else
    % 重开前检查缺失回执、已改动面表及已改动初场记录。只改试验副本。
    ready=fig.UserData.workflowState;bad=damagedJobs(job,out,ready);
    c.shortDirectory.Value=job;fprintf('Reopening completed calculation in a new window\n');press(c.reopenShort);
end
% 模块二：核对实际完成的终帧、反力、来源指纹、七站及全表面读数。
s=fig.UserData.workflowState;
assert(s.currentStep==5 && strcmp(s.status,'solved'),'%s',fig.UserData.status.Text);
m=s.result;r=fig.UserData.caseReport;j=m.shortSupportJob;
assert(strcmp(m.shortSupportAction,mode) && ~m.evidenceIdentity.reusedArchive);
f=shortSupportJobFrame(j,m.shortSupportImport,s.freestream);
assert(strcmp(f.terminal,m.evidenceIdentity.sourceFrame.sha256) && strcmp(f.initial,m.evidenceIdentity.sourceFrame.initialSHA256));
assert(strcmp(s.freestream.initialFrame,f.initial));
ui=fig.UserData.webUI;t=tic;
while isempty(ui.html.UserData.rendered) && toc(t)<45,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.faces==32480 && gl.webglError==0 && gl.visiblePixelSamples>0);
assert(r.fullSurfaceComputed && r.jointComputedFaceCount==32480 && r.jointExperimentalAreaM2==0);
ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;ev=ui.evidenceView;ev.selectField('heatFlux');
ids=[r.nativePeakFileFaceID;r.stations.fileFaceID];
for i=ids',ev.inspectFace(i,[]);v=ev.getState();assert(v.selectedFaceID==i);end
wall=ev.exportWallPairTo(fullfile(out,'wall.csv'));assert(wall.rows==32480 && wall.joint.pressureHeatFluxQualifiedArea==0);
% 未参与求解的三个温度不能在查看或重开时成为有限数。
data=ev.getData();
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    field=data.fields(strcmp({data.fields.key},key{1}));assert(isscalar(field) && all(isnan(field.values)));
end
% 模块三：记录真实窗口的终态。新算和重开可以独立比对，位姿与面序均来自实际 STL。
proof=struct('action',mode,'directory',j,'source',m.evidenceIdentity.sourceFrame,'faces',wall.rows,'elapsedS',s.freestream.elapsedS, ...
    'dtS',s.freestream.dtS,'areaM2',r.jointComputedAreaM2,'outW',-r.wholeSurfaceHeatPowerW,'peak',r.nativePeakReferenceFaceID, ...
    'reverseFaces',numel(r.reverseFaces),'regionalOutW',[r.regions.outW],'stations',fig.UserData.controls.stations.Data, ...
    'stationRMSE',r.stationRMSE,'equivalentRelativeError',r.equivalentRelativeError,'pixels',gl.visiblePixelSamples,'webGLError',gl.webglError);
proof.rejectedDamagedJobs=bad;
if strcmp(mode,'reopen')
    % 将本地重开的同一快照导出为离线三维页，可在浏览器检查原面和测量单位。
    exportAerothermalWeb(m,fullfile(out,'workbench'));
end
file=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(file>=0);guard=onCleanup(@()fclose(file)); %#ok<NASGU>
fprintf(file,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function results=damagedJobs(job,out,state)
% 模块用途：用与当前面网、边界相同的任务副本执行失败路径，并证实原任务可继续打开。
folder=fullfile(out,'damaged');mkdir(folder);
request=@()aeroHeatWorkflowState('beginSolve',state,struct('shortSupportAction','reopen','runDirectory',folder));
results={reject(request,'RAH:Short:Incomplete')};
mkdir(fullfile(folder,'export','full'));mkdir(fullfile(folder,'export','half'));
copyfile(fullfile(job,'completed.json'),folder);copyfile(fullfile(job,'construction.json'),folder);
copyfile(fullfile(job,'export','pair.json'),fullfile(folder,'export'));
copyfile(fullfile(job,'export','half','surface.zip'),fullfile(folder,'export','half'));
p=fullfile(folder,'export','full','surface.zip');file=fopen(p,'wb');assert(file>=0);fwrite(file,0,'uint8');fclose(file);
results{end+1}=reject(request,'RAH:Short:Job');
copyfile(fullfile(job,'export','full','surface.zip'),p,'f');
file=fopen(fullfile(folder,'construction.json'),'a');assert(file>=0);fprintf(file,' ');fclose(file);
results{end+1}=reject(request,'RAH:Short:Job');
end
function actual=reject(run,expected)
try,run();catch ex,actual=ex.identifier;assert(strcmp(actual,expected),'%s',ex.message);return;end
error('应拒绝的任务通过了来源检查。');
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
