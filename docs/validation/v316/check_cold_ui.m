function check_cold_ui(out,mode,job)
% 模块一：从模型和网格选择进入 12 K 计算页，核对原有来源和新增操作。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);c=fig.UserData.controls.leadGrid;c.Value=4;c.ValueChangedFcn(c,[]);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);
assert(fig.UserData.workflowState.boundaryConfirmed && fig.UserData.workflowState.mesh.thermistorPrepared.cells==558048);
assert(isfield(fig.UserData.controls,'coldMeasured') && isfield(fig.UserData.controls,'resumeBead'));
if strcmp(mode,'controls')
    % 将暂停请求交给实际按钮，在读入旧热态或开始重复推进之前保存一个新生成的冷初场。
    t=timer('ExecutionMode','fixedSpacing','Period',.2,'TimerFcn',@requestPause);
    closeTimer=onCleanup(@()delete(t));start(t);press(fig.UserData.controls.coldMeasured);stop(t); %#ok<NASGU>
    s=fig.UserData.workflowState;assert(strcmp(s.status,'paused') && isempty(s.result),fig.UserData.status.Text);
    folder=fig.UserData.controls.beadDirectory.Value;
    assert(isfile(fullfile(folder,'cold','cold.json')) && ~isfile(fullfile(folder,'completed.json')));
    proof=struct('pausedDirectory',folder,'paused',jsondecode(fileread(fullfile(folder,'progress.json'))));
    for origin={'archive','solve'}
        if strcmp(origin{1},'archive'),press(fig.UserData.controls.readBead);else,press(fig.UserData.controls.run);end
        s=fig.UserData.workflowState;assert(s.currentStep==5 && strcmp(s.result.thermistorMeasuredPoint.action,origin{1}),fig.UserData.status.Text);
        proof.(origin{1})=s.result.evidenceIdentity.sourceFrame.sha256;
        press(fig.UserData.workflowButtons(4));
    end
else
    % 模块二：用同一侧栏接回本次冷启动，检验终帧、全表面积、峰值、截面及导出的同面来源。
    fig.UserData.controls.beadDirectory.Value=job;press(fig.UserData.controls.resumeBead);s=fig.UserData.workflowState;
    assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
    m=s.result.thermistorMeasuredPoint;id=s.result.evidenceIdentity;
    canonicalJob=char(java.io.File(job).getCanonicalPath());
    assert(strcmp(m.action,'resume') && ~id.reusedArchive && strcmpi(id.runDirectory,canonicalJob) && m.jointVerifiedM2==0, ...
        '冷态动作 %s，归档标记 %d，实际目录 %s，预期目录 %s。',m.action,id.reusedArchive,id.runDirectory,canonicalJob);
    v=fig.UserData.webUI.evidenceView;d=v.getData();c=fig.UserData.controls.measuredPoint;
    assert(strcmp(c.Data{2,2},sprintf('%.2f',m.beadHWm2K)) && strcmp(d.measuredPoint.coldSource.frame,id.sourceFrame.sha256));
    html=fig.UserData.webUI.html;clock=tic;
    while isempty(html.UserData.rendered) && toc(clock)<40,drawnow;pause(.1);end
    gl=html.UserData.rendered;assert(gl.webglError==0 && gl.faces==10272 && gl.visiblePixelSamples>0);
    v.selectField('heatFlux');v.inspectFace(m.peakFaceID,[]);u=v.getState();assert(u.selectedFaceID==m.peakFaceID);
    w=v.exportWallPairTo(fullfile(out,'wall.csv'));q=v.exportTo(fullfile(out,'heat.csv'));
    assert(w.rows==10272 && q.rows==10272 && w.joint.pressureHeatFluxQualifiedArea==0);
    x=.000172;v.setSectionX(x);clock=tic;
    while toc(clock)<40
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-12);
    v.exportSectionTo(fullfile(out,'section.csv'));page=exportAerothermalWeb(s.result,fullfile(out,'workbench'));
    proof=struct('page',page,'frame',id.sourceFrame.sha256,'wallSHA256',id.sourceFrame.wallSHA256,'peakFace',u.selectedFaceID, ...
        'faces',gl.faces,'pixels',gl.visiblePixelSamples,'webglError',gl.webglError,'sectionXM',x,'job',job,'areaM2',m.areaM2,'h',m.beadHWm2K);
end
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
    function requestPause(timerObj,~)
        if ~fig.UserData.workflowState.busy,return;end
        controls=fig.UserData.controls;
        if isfield(controls,'pauseMeasured') && strcmp(controls.pauseMeasured.Enable,'on')
            stop(timerObj);press(controls.pauseMeasured);
        end
    end
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
