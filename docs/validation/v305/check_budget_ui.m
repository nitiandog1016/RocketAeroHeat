function check_budget_ui(out)
% 模块一：经实际侧栏打开第三档给定壁温结果，验证分带云图、逐带探针及同一色标。
fig=rocketAeroHeatWorkbench();guard=onCleanup(@()delete(fig)); %#ok<NASGU>
press(fig.UserData.controls.wiredBeadModel);choose(fig.UserData.controls.leadGrid,2);
press(fig.UserData.controls.generate);press(fig.UserData.controls.nextBoundary);
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);
s=fig.UserData.workflowState;assert(s.currentStep==5 && isfield(s.result,'thermistorLayerBudget'),fig.UserData.status.Text);
web=fig.UserData.webUI;v=web.evidenceView;d=v.getData();b=s.result.thermistorLayerBudget;
t=tic;while isempty(web.html.UserData.rendered) && toc(t)<35,drawnow;pause(.1);end
gl=web.html.UserData.rendered;assert(gl.faces==10272 && gl.visiblePixelSamples>0 && gl.webglError==0);
press(fig.UserData.controls.layerBudget);view=v.getState();peakFace=view.selectedFaceID;
assert(strcmp(view.fieldKey,'bandRelease') && b.faceBand(peakFace)==6);
delta=[b.bands.delta];keys={'bandRelease','bandConvection','bandDiffusion','bandSide'};
watts={[delta.faceOutW],[delta.outerConvW],[delta.outerDiffW],[delta.sideConvW]+[delta.sideDiffW]};
for j=1:numel(keys)
    f=d.fields(strcmp({d.fields.key},keys{j}));assert(all(f.status==1) && strcmp(f.unit,'W/m²'));
    assert(isequal(d.comparisonColorRanges.bandRelease,d.comparisonColorRanges.(keys{j})));
    for k=1:numel(b.bands)
        rows=b.faceBand==k;face=find(rows,1);v.selectField(keys{j});v.inspectFace(face,[]);
        ss=v.getState();assert(ss.selectedFaceID==face && strcmp(ss.fieldKey,keys{j}));
        assert(all(abs(f.values(rows)*b.bands(k).areaM2-watts{j}(k))<1e-13));
        assert(contains(f.scopeNames{f.scopeIndex(face)+1},b.bands(k).name));
    end
end
% 模块二：导出最大变化所含的前中段截面、可追溯的全表面表及独立离线页。
press(fig.UserData.controls.layerBudget);v.exportTo(fullfile(out,'faces.csv'));
x=.000172;v.setSectionX(x);t=tic;
while toc(t)<35
    drawnow;cut=v.getSection();
    if isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'bandRelease') && abs(cut.x-x)<1e-12,break;end
    pause(.1);
end
assert(isstruct(cut) && isfield(cut,'csv') && strcmp(cut.fieldKey,'bandRelease') && abs(cut.x-x)<1e-12);
v.exportSectionTo(fullfile(out,'section.csv'));rows=readtable(fullfile(out,'section.csv'),'Encoding','UTF-8');
assert(height(rows)>0 && all(b.faceBand(rows.face_id)==6));
exportAerothermalWeb(s.result,fullfile(out,'workbench'));
% 模块三：改用 4 mW 时原页面失效，复查功率工况仅显示自身已验收的结果。
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.beadCondition,'power');
assert(isempty(fig.UserData.workflowState.result));
press(fig.UserData.controls.confirm);press(fig.UserData.controls.readBead);p=fig.UserData.workflowState.result;
assert(isfield(p,'hybridSolvedWallTemperatureK') && ~isfield(p,'thermistorLayerBudget') && ...
    ~isfield(fig.UserData.controls,'layerBudget') && ~isfield(fig.UserData.webUI.evidenceView.getData(),'layerBudget'));
assert(abs(p.hybridCondition.netHeatingPowerW-.004)<1e-14 && abs(p.hybridSolvedWallTemperatureK(1)-297.694198579349)<1e-9);
proof=struct('faces',gl.faces,'samples',gl.visiblePixelSamples,'webglError',gl.webglError,'bands',numel(b.bands), ...
    'fields',{keys},'faceSelectedByButton',peakFace,'sectionX',x,'sectionRows',height(rows),'powerCaseSeparated',true);
fid=fopen(fullfile(out,'ui.json'),'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s\n',jsonencode(proof));fclose(fid);
end
function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
function choose(control,value)
control.Value=value;control.ValueChangedFcn(control,[]);drawnow;
end
