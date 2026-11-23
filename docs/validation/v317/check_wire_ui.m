function check_wire_ui(out,job)
% 模块一：在同一侧栏选择长引线边界，冻结本次热条件并读取实际原生任务。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,4);press(fig.UserData.controls.generate);
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.wireBoundary,'longWire');
press(fig.UserData.controls.confirm);st=fig.UserData.workflowState;
assert(st.boundaryConfirmed && strcmp(st.cfg.thermistorBoundary.leadThermalModel,'longWire'),fig.UserData.status.Text);
% 实际新算入口先冻结目录与均温来源，在首个阶段保存，避免为界面测试重复推进相同全场。
profile=jsondecode(fileread(fullfile(job,'prepared','profile.json')));
fig.UserData.controls.wireSeed.Value=fileparts(profile.base);
timerObj=timer('ExecutionMode','fixedSpacing','Period',.25,'TimerFcn',@pauseNew);
timerGuard=onCleanup(@()delete(timerObj));start(timerObj);press(fig.UserData.controls.run);stop(timerObj); %#ok<NASGU>
st=fig.UserData.workflowState;assert(strcmp(st.status,'paused') && isempty(st.result),fig.UserData.status.Text);
pausedDirectory=fig.UserData.controls.beadDirectory.Value;
saved=jsondecode(fileread(fullfile(pausedDirectory,'wire_job.json')));assert(strcmpi(saved.externalSeed,fileparts(profile.base)));
assert(~isfile(fullfile(pausedDirectory,'run','completed.json')));
fig.UserData.controls.beadDirectory.Value=job;press(fig.UserData.controls.resumeWire);
st=fig.UserData.workflowState;assert(st.currentStep==5 && ~isempty(st.result),fig.UserData.status.Text);
record=st.result.thermistorWire;assert(~st.result.evidenceIdentity.reusedArchive && strcmp(record.action,'resumeWire'));
v=fig.UserData.webUI.evidenceView;data=v.getData();tab=fig.UserData.controls.wireComparison;
assert(strcmp(tab.Data{3,2},sprintf('%.2f',record.beadHWm2K)) && strcmp(data.resultCondition.wallInputType,'perFace'));
powerTable=fig.UserData.controls.wireBudget;assert(size(powerTable.Data,1)==5 && strcmp(powerTable.Data{1,2},sprintf('%.5f',1000*record.profile.oneLeadRootW)));
choose(fig.UserData.controls.wireBudgetMode,'bands');
assert(size(powerTable.Data,1)==12 && strcmp(powerTable.Data{12,4},sprintf('%.3f',1e6*record.bands(12).frontSideW)));
assert(~isfield(st.result,'hybridSolvedWallTemperatureK') && strcmp(record.sourceFrameSHA256,data.thermistorWire.sourceFrameSHA256));
% 模块二：检查原面探针、非均温输入和热流差截线，再从归档重开同一终帧。
v.selectField('wallInput');v.inspectFace(record.peakFaceID,[]);u=v.getState();assert(u.selectedFaceID==record.peakFaceID);
html=fig.UserData.webUI.html;t=tic;
while isempty(html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
gl=html.UserData.rendered;assert(gl.webglError==0 && gl.visiblePixelSamples>0 && gl.faces==10272);
wall=v.exportWallPairTo(fullfile(out,'wall.csv'));assert(wall.rows==10272 && wall.joint.pressureHeatFluxQualifiedArea==0);
single=v.exportTo(fullfile(out,'temperature.csv'));assert(single.rows==10272);
v.selectField('wireHeatFluxChange');x=.00073;v.setSectionX(x);t=tic;
while toc(t)<40
    drawnow;cut=v.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'wireHeatFluxChange') && abs(cut.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'wireHeatFluxChange') && abs(cut.x-x)<1e-12);
v.exportSectionTo(fullfile(out,'section.csv'));page=exportAerothermalWeb(st.result,fullfile(out,'workbench'));
press(fig.UserData.workflowButtons(4));press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && s.result.evidenceIdentity.reusedArchive && isequaln(s.result.hybridWall,st.result.hybridWall) && ...
    isequal(s.result.hybridWallTemperatureK,st.result.hybridWallTemperatureK),fig.UserData.status.Text);
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.wireBoundary,'uniform');s=fig.UserData.workflowState;
assert(isempty(s.result) && ~s.boundaryConfirmed && ~isfield(s.cfg.thermistorBoundary,'leadThermalModel'));
proof=struct('page',page,'frame',record.sourceFrameSHA256,'faces',gl.faces,'pixels',gl.visiblePixelSamples,'webglError',gl.webglError, ...
    'peakFace',u.selectedFaceID,'sectionXM',x,'pausedNewJob',pausedDirectory,'bands',numel(record.bands), ...
    'archiveSameWall',true,'changedBoundaryInvalidated',true);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
    function pauseNew(timerHandle,~)
        if ~fig.UserData.workflowState.busy,return;end
        cs=fig.UserData.controls;
        if isfield(cs,'pauseMeasured') && strcmp(cs.pauseMeasured.Enable,'on')
            stop(timerHandle);press(cs.pauseMeasured);
        end
    end
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
