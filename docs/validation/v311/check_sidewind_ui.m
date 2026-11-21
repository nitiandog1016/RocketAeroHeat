function check_sidewind_ui(out)
% 模块一：从模型页进入全周侧吹，实际生成、确认边界并打开归档。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,3);
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh),fig.UserData.status.Text);
assert(s.mesh.thermistorPrepared.cells==278688 && s.cfg.thermistorAngularCells==48);
press(fig.UserData.controls.nextBoundary);press(fig.UserData.controls.confirm);
assert(fig.UserData.workflowState.boundaryConfirmed,fig.UserData.status.Text);
press(fig.UserData.controls.readBead);s=fig.UserData.workflowState;
assert(s.currentStep==5 && ~isempty(s.result),fig.UserData.status.Text);
r=s.result;web=fig.UserData.webUI;v=web.evidenceView;d=v.getData();stats=r.thermistorSidewind;
assert(isequal(s.freestream.velocityMps,[0 -3 0]) && isequal(r.hybridCondition.velocityMps,[0 -3 0]));
assert(contains(d.coordinateNote,'来流 −Y') && contains(d.caseNote,'侧吹') && d.resultCondition.wallUniformValueK==299.15);
t=tic;while isempty(web.html.UserData.rendered) && toc(t)<40,drawnow;pause(.1);end
gl=web.html.UserData.rendered;assert(gl.faces==10272 && gl.webglError==0 && gl.visiblePixelSamples>0);
% 模块二：检查同一截面的两个变量；导出面向科研核对的有符号数值。
x=.000172;v.inspectFace(stats.peakFaceID,[]);assert(v.getState().selectedFaceID==stats.peakFaceID);
pair=v.exportWallPairTo(fullfile(out,'wall.csv'));assert(pair.joint.pressureHeatFluxQualifiedArea==0);
for key={'pressureGauge','heatFlux'}
    field=key{1};v.selectField(field);v.setSectionX(x);t=tic;
    while toc(t)<40
        drawnow;cut=v.getSection();
        if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,field) && abs(cut.x-x)<1e-12,break;end
        pause(.1);
    end
    assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,field) && abs(cut.x-x)<1e-12);
    path=fullfile(out,[field '_section.csv']);v.exportSectionTo(path);rows=readtable(path,'Encoding','UTF-8');
    values=d.fields(strcmp({d.fields.key},field)).values;
    assert(height(rows)>0 && max(rows.value)>min(rows.value) && max(abs(rows.value-values(rows.face_id)))<1e-9);
end
page=exportAerothermalWeb(r,fullfile(out,'workbench'));
% 模块三：切回轴向网格后使旧的侧吹结果及热边界失效。
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.leadGrid,1);s=fig.UserData.workflowState;
assert(isempty(s.mesh) && isempty(s.result) && isempty(s.cfg.thermistorBoundary) && s.cfg.thermistorAngularCells==4);
proof=struct('page',page,'faces',gl.faces,'visibleSamples',gl.visiblePixelSamples,'webglError',gl.webglError, ...
    'peakFaceID',stats.peakFaceID,'meshSHA256',r.evidenceIdentity.meshSHA256,'sourceFrame',stats.sourceFrameSHA256, ...
    'wallPowerW',stats.totalOutwardW,'sectionX',x,'jointExperimentVerifiedM2',pair.joint.pressureHeatFluxQualifiedArea);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof,'PrettyPrint',true));fclose(fid);
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
