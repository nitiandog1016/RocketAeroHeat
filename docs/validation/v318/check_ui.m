function check_ui(out)
% 模块一：选择 B 珠的内部电阻测温边界，确认中间求解或浏览不会改写该温度身份。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,4);press(fig.UserData.controls.generate);
press(fig.UserData.controls.nextBoundary);choose(fig.UserData.controls.wireBoundary,'resistance');press(fig.UserData.controls.confirm);
s=fig.UserData.workflowState;assert(s.boundaryConfirmed && s.cfg.thermistorBoundary.resistanceTemperatureK==305.15 && ~isfield(s.cfg.thermistorBoundary,'wallTemperatureK'));
for name={'solve','archive'}
    action=name{1};folder=fullfile(out,action);mkdir(folder);
    if strcmp(action,'solve'),press(fig.UserData.controls.run);else,press(fig.UserData.controls.readBead);end
    s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
    w=s.result.thermistorWire;m=w.thermometry;v=fig.UserData.webUI.evidenceView;d=v.getData();
    assert(strcmp(w.action,action) && w.jointExperimentVerifiedM2==0 && m.internalK==305.15 && m.surfaceK<m.internalK);
    assert(strcmp(d.resultCondition.wallInputType,'perFace') && abs(d.resultCondition.wallMaxK-m.surfaceK)<1e-9 && ...
        ~isfield(s.result,'hybridSolvedWallTemperatureK') && strcmp(s.result.evidenceIdentity.sourceFrame.temperatureMethod,m.fieldMethod));
    rows=fig.UserData.controls.wireComparison.Data;
    assert(strcmp(rows{3,2},sprintf('%.2f',m.apparentHWm2K)) && strcmp(rows{4,2},sprintf('%.2f',m.surfaceHWm2K)));
    b=fig.UserData.controls.wireBudget;assert(strcmp(b.Data{2,2},sprintf('%.4f',m.surfaceK)));
    choose(fig.UserData.controls.wireBudgetMode,'bands');assert(size(b.Data,1)==12);
    % 模块二：用实际 WebGL 结果页检查两个不同的来源以及截线、探针和逐面导出。
    html=fig.UserData.webUI.html;t=tic;
    while isempty(html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
    gl=html.UserData.rendered;assert(gl.webglError==0 && gl.faces==10272 && gl.visiblePixelSamples>0);
    v.selectField('wallInput');v.inspectFace(w.peakFaceID,[]);u=v.getState();assert(u.selectedFaceID==w.peakFaceID);
    a=v.exportWallPairTo(fullfile(folder,'wall.csv'));assert(a.rows==10272 && a.joint.pressureHeatFluxQualifiedArea==0);
    v.selectField('thermometryHeatFluxChange');result=v.exportTo(fullfile(folder,'change.csv'));assert(result.rows==10272);
    x=.000172;v.setSectionX(x);t=tic;
    while toc(t)<40
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'thermometryHeatFluxChange') && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'thermometryHeatFluxChange') && abs(cut.x-x)<1e-12);
    v.exportSectionTo(fullfile(folder,'section.csv'));page=exportAerothermalWeb(s.result,fullfile(folder,'workbench'));
    proof.(action)=struct('page',page,'sourceFrame',w.sourceFrameSHA256,'sourceFieldFactor',m.fieldFactor,'faces',gl.faces, ...
        'pixels',gl.visiblePixelSamples,'webglError',gl.webglError,'sectionXM',x,'peakFace',u.selectedFaceID,'model',m);
    if strcmp(action,'solve'),press(fig.UserData.workflowButtons(4));end
end
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.wireBoundary,'longWire');s=fig.UserData.workflowState;
assert(isempty(s.result) && ~s.boundaryConfirmed && ~isfield(s.cfg.thermistorBoundary,'resistanceTemperatureK'));
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && ~isfield(s.result.thermistorWire,'thermometry') && s.result.hybridCondition.beadWallK==305.15);
proof.originalFrame=s.result.evidenceIdentity.sourceFrame.sha256;
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);cl=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
