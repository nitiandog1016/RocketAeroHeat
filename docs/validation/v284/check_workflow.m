function check_workflow(out,fresh)
% 走真实主侧栏，在没有前置输入的步骤核对门禁，比较两个网格的预览及逐面结果。
fig=rocketAeroHeatWorkbench();cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.beadModel);s=fig.UserData.workflowState;
assert(s.currentStep==2 && s.steps(1).completed && ~s.steps(3).available && isempty(s.result));
assert(strcmp(s.modelSource.kind,'thermistorB') && size(s.orientedMesh.F,1)==5216);
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh),fig.UserData.status.Text);
assert(s.mesh.thermistorPrepared.cells==6144 && ~s.steps(4).available);
proof.grid15=s.mesh.thermistorPrepared.directory;showBoundary();
if fresh
    press(fig.UserData.controls.run);
else
    press(fig.UserData.controls.readBead);
end
[result15,rep]=checkResult(15);proof.realSolve=~rep.reusedArchive;
assert(proof.realSolve==fresh);
if fresh
    proof.output=rep.outputDirectory;
    proof.frame=result15.evidenceIdentity.sourceFrame;
    proof.powerW=rep.wholeSurfaceHeatPowerW;
    proof.view=verifyDisplay();
    saveReport();return;
end
proof.view15=verifyDisplay();

% 撤销边界后原网格仍在，旧云图不可进入；再次确认不会自动重算。
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.beadCondition,'none');
s=fig.UserData.workflowState;
assert(~s.steps(4).available && ~s.steps(5).available && isempty(s.result) && isempty(fig.UserData.webUI) && ~isempty(s.mesh));
choose(fig.UserData.controls.beadCondition,'published');press(fig.UserData.controls.confirm);
assert(~fig.UserData.workflowState.steps(5).available);
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.beadCells,4);
s=fig.UserData.workflowState;assert(isempty(s.mesh) && isempty(s.freestream) && ~any([s.steps(3:5).available]));
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh),fig.UserData.status.Text);p=s.mesh.thermistorPrepared;
assert(p.cells==12288 && p.wallFaces==256 && size(p.preview.F,1)>256);
payload=fig.UserData.geometryView.html.UserData.payload;
assert(isequal(payload.vertices,p.preview.V) && isequal(payload.faces,p.preview.F-1) && payload.showWire);
proof.grid7p5=p.directory;showBoundary();
% 当前版本与编号正确时也要拒绝来自另一周向网格的结果。
trial=aeroHeatWorkflowState('beginSolve',fig.UserData.workflowState);
wrong=trial.cfg;wrong.thermistorAngularCells=2;
try
    aeroHeatWorkflowState('setBoundary',fig.UserData.workflowState,struct('cfg',wrong));
    error('RAH:Test:Accepted','边界页跨级更换了离散网格。');
catch ex
    assert(strcmp(ex.identifier,'RocketAeroHeat:WorkflowInputHierarchy'),ex.message);
end
try
    aeroHeatWorkflowState('completeSolve',trial,struct('result',result15,'solveId',trial.solveId,'sourceRevision',trial.solveInputRevision));
    error('RAH:Test:Accepted','错误的周向网格通过了核对。');
catch ex
    assert(strcmp(ex.identifier,'RAH:Bead:Snapshot'),ex.message);
end
press(fig.UserData.controls.readBead);[fine,rep]=checkResult(7.5);assert(rep.reusedArchive);
proof.view7p5=verifyDisplay();proof.deltaPeak=max(abs(fine.beadAngularDelta.heatFluxWm2));
press(fig.UserData.workflowButtons(1));press(fig.UserData.controls.changeModel);
s=fig.UserData.workflowState;assert(isempty(s.rawMesh) && isempty(s.result) && ~any([s.steps(2:5).available]));
assert(isfield(fig.UserData.controls,'thermistor') && isfield(fig.UserData.controls,'load'));
choose(fig.UserData.controls.thermistorMode,'archive7p5');press(fig.UserData.controls.thermistor);
quick=fig.UserData.thermistorReview;
assert(~isempty(quick) && quick.html.UserData.payload.thermistorBead.cellAngleDeg==7.5 && isempty(fig.UserData.workflowState.result));
proof.legacyArchiveOpened=true;
proof.invalidatedBoundary=true;proof.invalidatedGrid=true;proof.rejectedWrongGrid=true;proof.returnedToModelPicker=true;
saveReport();

    function showBoundary()
        press(fig.UserData.workflowButtons(3));
        assert(strcmp(fig.UserData.controls.confirm.Enable,'off'));
        choose(fig.UserData.controls.beadCondition,'published');press(fig.UserData.controls.confirm);
        st=fig.UserData.workflowState;
        assert(st.currentStep==4 && st.boundaryConfirmed && ~st.steps(5).available);
        assert(isequal(st.freestream.velocityMps,[-3 0 0]) && st.freestream.temperature==293.15);
    end
    function [r,rp]=checkResult(angle)
        st=fig.UserData.workflowState;
        assert(st.currentStep==5 && ~isempty(st.result),fig.UserData.status.Text);
        r=st.result;rp=fig.UserData.caseReport;
        assert(r.thermistorBead.cellAngleDeg==angle && rp.fullSurfaceComputed && size(r.F,1)==5216);
        assert(strcmp(r.evidenceIdentity.meshSHA256,st.mesh.thermistorPrepared.meshSHA256));
        assert(isequal(r.hybridWall.faceID,(0:5215)') && all(r.hybridWallTemperatureK==305.15));
        assert(~r.evidenceIdentity.productionQualified);
    end
    function stats=verifyDisplay()
        ui=fig.UserData.webUI;t=tic;
        while isempty(ui.html.UserData.rendered) && toc(t)<35,drawnow;pause(.1);end
        gl=ui.html.UserData.rendered;assert(~isempty(gl) && gl.faces==5216 && gl.webglError==0 && gl.visiblePixelSamples>0);
        r=fig.UserData.workflowState.result;ui.evidenceView.inspectFace(5204,[]);
        path=fullfile(out,sprintf('wall_%g.csv',r.thermistorBead.cellAngleDeg));
        receipt=ui.evidenceView.exportWallPairTo(path);rows=readtable(path);
        assert(isequal(rows.face_id,r.hybridWall.faceID+1) && max(abs(rows.pressure*1000-r.hybridWall.pressurePa))<1e-8 ...
            && max(abs(rows.heat_flux*1000-r.hybridWall.heatFluxWm2))<1e-9);
        assert(receipt.joint.pressureHeatFluxQualifiedArea==0 && strcmp(receipt.identity.origin,'research'));
        stats=struct('render',gl,'csv',path,'rows',height(rows),'selectedFace',5204,'powerW',fig.UserData.caseReport.wholeSurfaceHeatPowerW);
    end
    function saveReport()
        fid=fopen(fullfile(out,'workflow.json'),'w','n','UTF-8');assert(fid>=0);
        closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
        fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
    end
end
function press(button)
% 使用实际按钮回调，让异步网页消息可进入 MATLAB 事件循环。
button.ButtonPushedFcn(button,[]);drawnow;
end
function choose(control,value)
control.Value=value;control.ValueChangedFcn(control,[]);drawnow;
end
