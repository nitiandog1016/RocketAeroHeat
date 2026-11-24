function check_ui(out,mode,job)
% 模块一：从米制普通 STL 确认网格、时间与步长，在真实窗口加载所选端点。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
c=fig.UserData.controls;c.file.Value=fullfile(out,'input.stl');c.scale.Value=1;press(c.load);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'shortSupport65');press(fig.UserData.controls.generate);
press(fig.UserData.controls.nextBoundary);early=strcmp(mode,'early');local=strcmp(mode,'local');
if early,branch='half';target=.25;elseif local,branch='full';target=1;else,branch='half';target=1;end
choose(fig.UserData.controls.shortTarget,target);choose(fig.UserData.controls.shortBranch,branch);
press(fig.UserData.controls.confirm);c=fig.UserData.controls;
if local
    % 新算请求另行验证其冻结的起点字段。实际数值运行通过相同后台入口完成。
    ready=fig.UserData.workflowState;origin=fullfile(fileparts(fileparts(job)),'v338_validation','reread01');
    new=aeroHeatWorkflowState('beginSolve',ready,struct('shortSupportAction','solve','shortStartDirectory',origin));
    assert(strcmp(new.frozenInputs.cfg.shortSupportSourceJob,char(java.io.File(origin).getCanonicalPath())));
    c.shortDirectory.Value=job;press(c.reopenShort);
else
    press(c.readShort);
end
state=fig.UserData.workflowState;assert(strcmp(state.status,'solved'),'%s',fig.UserData.status.Text);
m=state.result;c=fig.UserData.controls;ui=fig.UserData.webUI;r=fig.UserData.caseReport;
assert(r.jointComputedFaceCount==32480 && r.jointExperimentalAreaM2==0);
if ~early
    choose(c.stationMethod,'time');assert(isequal(c.stations.Data,r.stations{:,{'angleDeg','reportedNu','earlyNu','computedNu'}}));
    assert(isfield(m,'shortSupportTemporal') && m.shortSupportTemporal.laterTimeS==state.freestream.elapsedS);
    if local
        assert(isequal(m.shortSupportImport.referenceFaceID,(32480:-1:1)') && strcmp(m.shortSupportAction,'reopen'));
    else
        assert(isequal(m.shortSupportImport.referenceFaceID,(1:32480)') && m.evidenceIdentity.reusedArchive);
    end
    press(c.openAzimuth);tab=ui.tab.Parent.SelectedTab;timer=tic;
    while isempty(tab.UserData.rendered) && toc(timer)<35,drawnow;pause(.1);end
    msg=tab.UserData.rendered;assert(strcmp(msg.frame,m.evidenceIdentity.sourceFrame.sha256));
    assert(max(abs(msg.values(:)-r.stations.computedNu(end-1:end)))<1e-8);
else
    assert(~isfield(m,'shortSupportTemporal') && state.freestream.elapsedS<.003);
end
% 模块二：读取实际渲染、热峰探针和逐面导出，验证温度身份及冻结后的旧时刻拒收。
ui.tab.Parent.SelectedTab=ui.tab;timer=tic;
while isempty(ui.html.UserData.rendered) && toc(timer)<45,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(gl.faces==32480 && gl.webglError==0 && gl.visiblePixelSamples>0);
ui.tab.Parent.SelectedTab=ui.evidenceTab;ev=ui.evidenceView;d=ev.getData();
for key={'wallSolved','edgeTemperature','recoveryTemperature'}
    f=d.fields(strcmp({d.fields.key},key{1}));assert(all(isnan(f.values)));
end
if ~early
    for name={'timeGaugeChange','timeFirstChange','timeFluxChange'}
        key=name{1};ev.selectField(key);ev.inspectFace(r.peakOutFace,[]);probe=ev.getState();
        assert(probe.selectedFaceID==r.peakOutFace && strcmp(probe.fieldKey,key));
    end
    x=m.shortSupportImport.translationM(1)+.00137;ev.setSectionX(x);timer=tic;
    while toc(timer)<45
        drawnow;cut=ev.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'timeFluxChange') && abs(cut.x-x)<1e-13,break;end
        pause(.12);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'timeFluxChange') && abs(cut.x-x)<1e-13);
    ev.exportSectionTo(fullfile(out,'section.csv'));
    rec=ev.exportWallPairTo(fullfile(out,'wall.csv'));assert(rec.rows==32480);
    writetable(r.stations,fullfile(out,'stations.csv'),'Encoding','UTF-8');
    exportAerothermalWeb(m,fullfile(out,'workbench'));
end
proof=struct('mode',mode,'branch',branch,'timeS',state.freestream.elapsedS,'dtS',state.freestream.dtS,'faces',gl.faces, ...
    'pixels',gl.visiblePixelSamples,'terminal',m.evidenceIdentity.sourceFrame.sha256,'previous',m.evidenceIdentity.sourceFrame.previousSHA256, ...
    'peakFile',r.nativePeakFileFaceID,'peakReference',r.nativePeakReferenceFaceID,'outW',-r.wholeSurfaceHeatPowerW,'regionsW',[r.regions.outW]);
if ~early
    proof.powerDeltaW=m.shortSupportTemporal.powerChangeW;proof.earlierFrame=m.shortSupportTemporal.baseFrame;proof.azimuth=msg;
    press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.shortTarget,.25);
    assert(isempty(fig.UserData.workflowState.result));press(fig.UserData.controls.confirm);ready=fig.UserData.workflowState;
    rejected=false;try,aeroHeatWorkflowState('beginSolve',ready,struct('shortSupportAction','reopen','runDirectory',job));
    catch ex,rejected=strcmp(ex.identifier,'RAH:Short:Job');end
    assert(rejected);proof.wrongTimeRejected=rejected;
end
f=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(f>=0);close=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,val)
c.Value=val;c.ValueChangedFcn(c,[]);drawnow;
end
