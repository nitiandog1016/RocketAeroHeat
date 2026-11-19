function check_workflow(out)
%CHECK_WORKFLOW 用第二份 STL 的真实五步控件计算，并从同一结果页保留绘制与选面回执。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
cfg=createCaseConfig();cfg.output.verbose=false;cfg.geometry.scale=1;
cfg.geometry.stlFile=fullfile(root,'cases','CM011-alternate-SI.stl');
fig=rocketAeroHeatWorkbench(cfg);cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
fig.Visible='off';press(fig.UserData.controls.load);
s=fig.UserData.workflowState;assert(size(s.rawMesh.F,1)==10080,fig.UserData.status.Text);
press(fig.UserData.controls.confirm);choose(fig.UserData.controls.meshBackend,'cm011Imported');
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh)&&isfield(s.mesh,'cm011Prepared'),fig.UserData.status.Text);
p=s.mesh.cm011Prepared;saveText(out,'mesh_folder.txt',p.directory);
assert(isequal(s.mesh.F,s.rawMesh.F)&&isequal(s.mesh.V,s.rawMesh.V)&&p.maxFacetSampleDistanceM<5e-6);
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.cm011Condition,'inclined');
press(fig.UserData.controls.confirm);assert(fig.UserData.workflowState.boundaryConfirmed);
press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(strcmp(s.status,'solved')&&s.currentStep==5,fig.UserData.status.Text);
assert(isequal(s.result.V,s.rawMesh.V)&&isequal(s.result.F,s.rawMesh.F)&&s.resultSourceRevision==s.revision);
report=fig.UserData.caseReport;assert(report.fullSurfaceComputed&&report.native.acceptedNumericalWall&&~report.reusedNativeRun);
saveText(out,'run_folder.txt',report.outputDirectory);
ui=fig.UserData.webUI;d=ui.html.UserData.payload;n=size(s.mesh.F,1);
assert(isequal(d.faces,s.rawMesh.F-1)&&isequal(s.result.hybridWall.faceID,(0:n-1)')&&d.resultCondition.alphaDeg==10);
assert(~d.researchIdentity.productionQualified && strcmp(d.researchIdentity.stlSHA256,s.modelSource.sha256));
page=exportAerothermalWeb(s.result,fullfile(out,'workbench'));assert(isfile(page));
started=tic;while toc(started)<16 && isempty(ui.html.UserData.rendered),drawnow;pause(.1);end
rendered=ui.html.UserData.rendered;
ids=[1,1701,n];q=d.fields(strcmp({d.fields.key},'heatFlux'));
for id=ids
    ui.html.Data=struct('kind','selection','face',id);ui.html.DataChangedFcn(ui.html,[]);
    assert(ui.html.UserData.selectedFace==id && abs(q.values(id)*1000-s.result.hybridWall.heatFluxWm2(id))<1e-14);
end
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.cm011Condition,'none');
assert(isempty(fig.UserData.webUI)&&isempty(fig.UserData.workflowState.result));
saveText(out,'workflow.json',jsonencode(struct('originalFaces',n,'sourceRevision',double(s.revision),'solveId',double(s.solveId), ...
    'areaM2',report.totalArea,'powerW',report.wholeSurfaceHeatPowerW,'rendered',rendered,'pickedByEvent',ids,'mouseDriven',false)));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function choose(c,v)
c.Value=v;c.ValueChangedFcn(c,[]);drawnow;
end
function saveText(folder,name,value)
fid=fopen(fullfile(folder,name),'w','n','UTF-8');assert(fid>=0);clean=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',value);
end
