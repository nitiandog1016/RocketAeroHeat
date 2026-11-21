function check_three_ui(out)
% 模块一：在生产侧栏选择第三档网、生成体网，分别确认两种热边界并打开原面结果。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,2);
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(s.mesh.thermistorPrepared.cells==42224 && s.cfg.thermistorRefineLeads==2);
press(fig.UserData.controls.nextBoundary);
for k=1:2
    mode='givenWall';
    if k==2
        mode='givenPower';press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.beadCondition,'power');
        assert(isempty(fig.UserData.workflowState.result) && fig.UserData.controls.beadPower.Value==.004);
    end
    press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
    assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
    r=s.result;web=fig.UserData.webUI;v=web.evidenceView;d=v.getData();g=r.thermistorThreeGrids;summary=g.summary;
    assert(summary.available && strcmp(summary.mode,mode) && numel(summary.grids)==3 && numel(summary.regions)==7);
    if k==1,assert(all([summary.grids.wallK]==305.15));else,assert(all(abs([summary.grids.totalOutwardW]-.004)<1e-12));end
    assert(isfield(r,'hybridSolvedWallTemperatureK')==(k==2));
    t=tic;while isempty(web.html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
    gl=web.html.UserData.rendered;assert(gl.faces==10272 && gl.webglError==0 && gl.visiblePixelSamples>0);
    ranges=d.comparisonColorRanges;assert(isequal(ranges.gridHeat1,ranges.gridHeat2,ranges.heatFlux));
    assert(isequal(ranges.gridHeat21,ranges.gridHeat32));
    % 模块二：拾取最大热流差的原面，导出同一字段的整表、截面及离线三维页。
    press(fig.UserData.controls.compareGrids);q=g.heat(:,3)-g.heat(:,2);
    [~,peak]=max(abs(q));state=v.getState();assert(state.selectedFaceID==peak && strcmp(state.fieldKey,'gridHeat32'));
    folder=fullfile(out,mode);mkdir(folder);v.exportTo(fullfile(folder,'faces.csv'));
    x=.000172;v.setSectionX(x);t=tic;
    while toc(t)<40
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'gridHeat32') && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'gridHeat32') && abs(cut.x-x)<1e-12);
    path=fullfile(folder,'section.csv');v.exportSectionTo(path);rows=readtable(path,'Encoding','UTF-8');
    assert(height(rows)>0 && max(abs(rows.value-q(rows.face_id)))<1e-9);
    page=exportAerothermalWeb(r,fullfile(folder,'workbench'));
    proof.(mode)=struct('page',page,'faces',gl.faces,'visibleSamples',gl.visiblePixelSamples,'webglError',gl.webglError, ...
        'probeFace',peak,'sectionX',x,'cutRows',height(rows),'summary',summary);
end
% 模块三：零功率时三档均返回环境温度和物理零热流，切回旧网后必须重新建网。
cfg=s.cfg;cfg.thermistorResultSource='archive';cfg.thermistorBoundary.netHeatingPowerW=0;
r=runAeroHeatSingleCase(s.mesh,cfg);g=r.thermistorThreeGrids;
assert(all(g.heat==0,'all') && all([g.summary.grids.wallK]==293.15));
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.leadGrid,1);
s=fig.UserData.workflowState;assert(isempty(s.mesh) && isempty(s.result) && isempty(s.cfg.thermistorBoundary));
proof.zeroPowerAndInvalidation=true;
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && isfield(s.result,'thermistorOuterDomain') && isfield(s.result,'thermistorLeadGrid'));
assert(strcmp(s.mesh.thermistorPrepared.meshSHA256,summary.grids(2).meshSHA256) && ...
    strcmp(s.result.evidenceIdentity.sourceFrame.sha256,summary.grids(2).frameSHA256));
proof.secondGridFrame=s.result.evidenceIdentity.sourceFrame.sha256;
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
