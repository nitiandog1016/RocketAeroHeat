function check_domain_ui(out,branch)
% 模块一：从普通 STL 先绑定 14D，再切到 20D，验证旧网格、边界与结果一起失效。
fig=rocketAeroHeatWorkbench();closeFig=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(out,'specimen.stl');c.scale.Value=1;
if strcmp(branch,'half'),c.scale.Value=.001;end
press(c.load);
if strcmp(branch,'half'),c=fig.UserData.controls;c.rz.Value=-90;press(c.rotate);end
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'shortSupport65');press(fig.UserData.controls.generate);
fprintf('14D STL bound\n');
s=fig.UserData.workflowState;assert(~isempty(s.mesh),'%s',fig.UserData.status.Text);binding14=s.mesh.shortSupportImport;
choose(fig.UserData.controls.shortRadius,20);s=fig.UserData.workflowState;
assert(isempty(s.mesh) && isempty(s.result) && ~s.boundaryConfirmed);
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;b=s.mesh.shortSupportImport;
fprintf('20D STL bound\n');
assert(~strcmp(b.meshSHA256,binding14.meshSHA256) && isequal(b.referenceFaceID,binding14.referenceFaceID));
order=readmatrix(fullfile(out,'expected_order.txt'));assert(isequal(order,b.referenceFaceID));
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.shortBranch,branch);press(fig.UserData.controls.confirm);
ready=fig.UserData.workflowState;assert(ready.currentStep==4 && ready.cfg.shortSupportRadiusD==20);
if strcmp(branch,'full')
    % 在进程边界截获一次新算请求，确认使用 20D 的程序和新目录；实际三步另有原生回执。
    folder=tempname;mkdir(folder);path=fullfile(folder,'runCM011Python.m');
    template=fullfile(fileparts(fileparts(mfilename('fullpath'))),'v334','dispatch_probe.m.in');
    fid=fopen(path,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',strrep(fileread(template),'shortSupportTime65.py','shortSupportOuter65.py'));fclose(fid);
    addpath(folder,'-begin');cleanPath=onCleanup(@()removeProbe(folder,path));
    press(fig.UserData.controls.run);failed=fig.UserData.workflowState;
    assert(strcmp(failed.error.identifier,'RAH:ShortTest:DispatchObserved') && isempty(failed.result));
    clear cleanPath;
end
press(fig.UserData.controls.readShort);s=fig.UserData.workflowState;
assert(s.currentStep==5 && strcmp(s.status,'solved'),'%s',fig.UserData.status.Text);
fprintf('20D archive accepted\n');
% 模块二：两个网的七站、每面的差量、参考帧、探针与截线属于同一快照。
report=fig.UserData.caseReport;result=s.result;c=fig.UserData.controls;ui=fig.UserData.webUI;
assert(isequal(c.stations.Data(:,4),report.stations.baseNu) && isequal(c.stations.Data(:,5),report.stations.computedNu));
timer=tic;while isempty(ui.html.UserData.rendered) && toc(timer)<45,drawnow;pause(.12);end
gl=ui.html.UserData.rendered;assert(gl.faces==32480 && gl.webglError==0 && gl.visiblePixelSamples>0);
fprintf('WebGL confirmed\n');
ui.evidenceTab.Parent.SelectedTab=ui.evidenceTab;ev=ui.evidenceView;data=ev.getData();
assert(contains(data.resultCondition.caption,'20D') && data.resultCondition.physicalStepSizeS==ready.freestream.dtS);
for key={'outerGaugeChange','outerFluxChange'}
    ev.selectField(key{1});f=data.fields(strcmp({data.fields.key},key{1}));
    [~,face]=max(abs(f.values));ev.inspectFace(face,[]);view=ev.getState();assert(view.selectedFaceID==face && strcmp(view.fieldKey,key{1}));
end
x=.00137+b.translationM(1);ev.setSectionX(x);timer=tic;
while toc(timer)<45
    drawnow;cut=ev.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'outerFluxChange') && abs(cut.x-x)<1e-13,break;end
    pause(.12);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'outerFluxChange') && abs(cut.x-x)<1e-13);
ev.exportSectionTo(fullfile(out,'section.csv'));wall=ev.exportWallPairTo(fullfile(out,'wall.csv'));
fprintf('Section and wall exported\n');
assert(wall.rows==32480 && wall.joint.pressureHeatFluxQualifiedArea==0);
% 模块三：保存两种测度的身份；离线页只导出一份，原生计算与界面单列。
proof=struct('branch',branch,'baseFrame',result.shortSupportOuter.baseFrame,'frame',result.evidenceIdentity.sourceFrame.sha256, ...
    'baseMesh',binding14.meshSHA256,'mesh',b.meshSHA256,'stlSHA256',b.source.sha256,'elapsedS',s.freestream.elapsedS, ...
    'dtS',s.freestream.dtS,'faceCount',wall.rows,'powerW',-report.wholeSurfaceHeatPowerW, ...
    'basePowerW',report.baseOutW,'peakFace',report.peakOutFace,'referencePeakFace',report.nativePeakReferenceFaceID, ...
    'stations',c.stations.Data,'pixels',gl.visiblePixelSamples,'glError',gl.webglError);
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.shortRadius,14);s=fig.UserData.workflowState;
assert(isempty(s.mesh) && isempty(s.result) && ~s.steps(5).available && ~s.boundaryConfirmed);
if strcmp(branch,'half'),exportAerothermalWeb(result,fullfile(out,'workbench'));end
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
function removeProbe(folder,path)
rmpath(folder);delete(path);rmdir(folder);
end
