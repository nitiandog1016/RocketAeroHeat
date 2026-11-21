function check_point_ui(out)
% 模块一：由模型、网格、边界进入计算页，核验工作台实际使用 96 档体网和 12 K 温差。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,4);press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;assert(~isempty(s.mesh),fig.UserData.status.Text);p=s.mesh.thermistorPrepared;
assert(p.angularCells==96 && p.cells==558048 && size(p.preview.F,1)==20544 && size(s.mesh.F,1)==10272);
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);s=fig.UserData.workflowState;
assert(s.boundaryConfirmed && s.cfg.thermistorBoundary.wallTemperatureK==305.15 && isequal(s.freestream.velocityMps,[0 -3 0]));
% 模块二：分别操作换算与直接读回按钮，核验面元、三维页、联合表和单变量来源回执。
for mode={'solve','archive'}
    action=mode{1};folder=fullfile(out,action);mkdir(folder);
    if strcmp(action,'solve'),press(fig.UserData.controls.run);else,press(fig.UserData.controls.readBead);end
    s=fig.UserData.workflowState;assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
    m=s.result.thermistorMeasuredPoint;assert(strcmp(m.action,action) && m.jointVerifiedM2==0 && numel(m.regions)==7);
    v=fig.UserData.webUI.evidenceView;d=v.getData();c=fig.UserData.controls.measuredPoint;
    assert(size(c.Data,1)==4 && strcmp(c.Data{2,2},sprintf('%.2f',m.beadHWm2K)) && d.resultCondition.wallUniformValueK==305.15);
    assert(~isfield(s.result,'hybridSolvedWallTemperatureK') && all(s.result.hybridWallTemperatureK==305.15));
    html=fig.UserData.webUI.html;t=tic;
    while isempty(html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
    gl=html.UserData.rendered;assert(gl.webglError==0 && gl.faces==10272 && gl.visiblePixelSamples>0);
    v.selectField('heatFlux');v.inspectFace(m.peakFaceID,[]);u=v.getState();assert(u.selectedFaceID==m.peakFaceID);
    w=v.exportWallPairTo(fullfile(folder,'wall.csv'));assert(w.rows==10272 && w.joint.pressureHeatFluxQualifiedArea==0);
    q=v.exportTo(fullfile(folder,'heat.csv'));assert(q.rows==10272);
    x=.000172;v.setSectionX(x);t=tic;
    while toc(t)<40
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'heatFlux') && abs(cut.x-x)<1e-12);
    v.exportSectionTo(fullfile(folder,'section.csv'));
    page=exportAerothermalWeb(s.result,fullfile(folder,'workbench'));
    proof.(action)=struct('page',page,'peakFace',m.peakFaceID,'frame',s.result.evidenceIdentity.sourceFrame.sha256, ...
        'faces',gl.faces,'pixels',gl.visiblePixelSamples,'webglError',gl.webglError,'sectionXM',x,'h',m.beadHWm2K);
    if strcmp(action,'solve'),press(fig.UserData.workflowButtons(4));end
end
% 模块三：切回原 6 K 网格，旧测点边界和两种来源的结果均失效。
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.leadGrid,3);s=fig.UserData.workflowState;
assert(isempty(s.mesh) && isempty(s.result) && isempty(s.cfg.thermistorBoundary) && s.cfg.thermistorAngularCells==48);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);guardFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
