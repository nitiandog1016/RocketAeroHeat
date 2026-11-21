function check_workbench(out,mode,folder)
% 模块一：经实际五步侧栏选用公开试件，核对全周网格预览、冻结边界与步骤门禁。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.beadModel);c=fig.UserData.controls;
c.beadCells.Value=48;c.beadCells.ValueChangedFcn(c.beadCells,[]);press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;assert(~isempty(s.mesh),fig.UserData.status.Text);
p=s.mesh.thermistorPrepared;assert(p.nodes==172482 && p.cells==172032 && p.wallFaces==3072);
d=fig.UserData.geometryView.html.UserData.payload;
assert(isequal(d.vertices,p.preview.V) && isequal(d.faces,p.preview.F-1) && ~s.steps(4).available);
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);s=fig.UserData.workflowState;
assert(s.currentStep==4 && s.freestream.speed==6 && s.cfg.thermistorBoundary.wallTemperatureK==303.15);
c=fig.UserData.controls;
switch mode
    case 'run'
        press(c.run);
    case 'resume'
        c.beadDirectory.Value=folder;press(c.resumeBead);
    otherwise
        c.beadDirectory.Value=folder;press(c.readBead);
end
s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
r=fig.UserData.caseReport;job=r.outputDirectory;
folder=char(java.io.File(folder).getCanonicalPath());
if strcmp(mode,'run'),assert(~strcmpi(job,folder));else,assert(strcmpi(job,folder));end

% 模块二：逐面与本次目录核对，检查实际三维页、物理时间、探针、符号与积分功率。
wall=readtable(fullfile(job,'second_energy','wall.csv'));w=s.result.hybridWall;
assert(isequal(w.faceID,wall.faceID) && all(abs(w.heatFluxWm2-wall.heatFluxWm2)<1e-9) && ...
    all(abs(w.pressurePa-wall.pressurePa)<1e-9) && all(abs(w.nearWallCellTemperatureK-wall.nearWallCellTemperatureK)<1e-9));
assert(abs(sum(wall.areaM2.*w.heatFluxWm2)-r.wholeSurfaceHeatPowerW)<1e-13 && isfile(r.workbenchPage));
ui=fig.UserData.webUI;t=tic;
while isempty(ui.html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
gl=ui.html.UserData.rendered;assert(~isempty(gl) && gl.faces==5216 && gl.webglError==0 && gl.visiblePixelSamples>0);
v=ui.evidenceView;v.selectField('heatFlux');v.inspectFace(5211,[]);d=v.getData();
assert(v.getState().selectedFaceID==5211 && d.resultCondition.elapsedPhysicalTimeS==7.80936086544745e-6);
assert(abs(d.evidenceLedger.joint.pressureHeatFluxAvailableArea-sum(wall.areaM2))<1e-15);
e=s.result.evidenceIdentity;
assert(strcmp(e.runDirectory,job) && strcmp(r.startupAction,mode) && ...
    strcmp(e.sourceFrame.wallSHA256,cm011FileSHA256(fullfile(job,'second_energy','wall.csv'))));
proof=struct('mode',mode,'directory',job,'physicalTimeS',d.resultCondition.elapsedPhysicalTimeS,'faces',height(w), ...
    'areaM2',sum(wall.areaM2),'intoWallW',r.wholeSurfaceHeatPowerW,'sourceFrame',e.sourceFrame.sha256, ...
    'wallSHA256',e.sourceFrame.wallSHA256,'webglError',gl.webglError,'visiblePixelSamples',gl.visiblePixelSamples, ...
    'pickedFace',v.getState().selectedFaceID,'page',r.workbenchPage,'prepared',p.manifest,'resultReceipt',r.resultReceipt);
if strcmp(mode,'read')
    % 以真实已保存端点验证目录身份和热边界拒绝路径，两项均不启动求解器。
    trial=aeroHeatWorkflowState('beginSolve',s,struct('startupAction','read','runDirectory',out));
    try
        aeroHeatWorkflowState('completeSolve',trial,struct('result',s.result,'sourceRevision',trial.revision,'solveId',trial.solveId));
        error('RAH:Test:Accepted','错误目录的终场被接收。');
    catch ex
        assert(strcmp(ex.identifier,'RAH:Bead:Snapshot'),ex.message);
    end
    wrong=s.cfg.thermistorBoundary;wrong.wallTemperatureK=304;
    try
        thermistorWorkbenchBoundary(wrong,48);error('RAH:Test:Accepted','错误壁温被接收。');
    catch ex
        assert(strcmp(ex.identifier,'RAH:Bead:Boundary'),ex.message);
    end
    proof.rejectedWrongDirectoryAndWallTemperature=true;
end
% 模块三：改变网格后旧的边界、结果和结果页同时失效。
press(fig.UserData.workflowButtons(2));c=fig.UserData.controls;c.beadCells.Value=2;c.beadCells.ValueChangedFcn(c.beadCells,[]);
s=fig.UserData.workflowState;assert(isempty(s.mesh) && isempty(s.result) && isempty(s.freestream) && ~s.steps(5).available);
proof.invalidatedOnGridChange=true;
fid=fopen(fullfile(out,'workflow.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end

function press(button)
% 触发生产界面的按钮回调，并处理计算完成时的实际网页消息。
button.ButtonPushedFcn(button,[]);drawnow;
end
