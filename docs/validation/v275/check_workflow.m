function check_workflow(out)
%CHECK_WORKFLOW 从实际五步控件发起一次低速计算，核对快照、导出及上游变更。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
cfg=createCaseConfig();cfg.geometry.stlFile=fullfile(root,'cases','CM011-mm-scale-SI.stl');cfg.geometry.scale=1;
cfg.output.verbose=false;fig=rocketAeroHeatWorkbench(cfg);
cleanup=onCleanup(@()deleteIfOpen(fig)); %#ok<NASGU>
fig.Visible='off';press(fig.UserData.controls.load);
s=fig.UserData.workflowState;assert(size(s.rawMesh.F,1)==14140,fig.UserData.status.Text);
press(fig.UserData.controls.confirm);setValue(fig.UserData.controls.meshBackend,'cm011Imported');
press(fig.UserData.controls.generate);s=fig.UserData.workflowState;
assert(~isempty(s.mesh)&&isfield(s.mesh,'cm011Prepared'),fig.UserData.status.Text);
p=s.mesh.cm011Prepared;
writeText(fullfile(out,'mesh_folder.txt'),p.directory);
assert(isequal(s.mesh.F,s.rawMesh.F)&&isequal(s.mesh.V,s.rawMesh.V));
assert(strcmp(cm011FileSHA256(fullfile(p.directory,'mesh.su2')),'adfd41231c7fba5a5d7f9469c69c0080f1e465143d528768cfb2092bcc7c19c0'));
press(fig.UserData.workflowButtons(3));assert(strcmp(fig.UserData.controls.confirm.Enable,'off'));
setValue(fig.UserData.controls.cm011Condition,'inclined');press(fig.UserData.controls.confirm);
s=fig.UserData.workflowState;assert(s.currentStep==4 && s.boundaryConfirmed && isempty(s.result));
press(fig.UserData.controls.run);s=fig.UserData.workflowState;
assert(strcmp(s.status,'solved')&&s.currentStep==5,fig.UserData.status.Text);
assert(isequal(s.result.V,s.rawMesh.V)&&isequal(s.result.F,s.rawMesh.F)&&s.resultSourceRevision==s.revision);
n=size(s.mesh.F,1);assert(isequal(s.result.hybridWall.faceID,(0:n-1)'));
report=fig.UserData.caseReport;assert(report.fullSurfaceComputed && report.native.acceptedNumericalWall && ~report.reusedNativeRun);
writeText(fullfile(out,'run_folder.txt'),report.outputDirectory);
ui=fig.UserData.webUI;assert(~isempty(ui));
d=ui.html.UserData.payload;assert(isequal(d.faces,s.rawMesh.F-1)&&d.resultCondition.alphaDeg==10);
assert(~d.researchIdentity.productionQualified&&strcmp(d.researchIdentity.stlSHA256,s.modelSource.sha256));
page=exportAerothermalWeb(s.result,fullfile(out,'workbench'));assert(isfile(page));
rendered=ui.html.UserData.rendered;
press(fig.UserData.workflowButtons(2));press(fig.UserData.workflowButtons(5));
assert(fig.UserData.workflowState.solveId==s.solveId && fig.UserData.workflowState.resultSourceRevision==s.revision);
press(fig.UserData.workflowButtons(3));setValue(fig.UserData.controls.cm011Condition,'none');
x=fig.UserData.workflowState;assert(isempty(x.result)&&~x.boundaryConfirmed&&~x.steps(5).available&&x.revision>s.revision);
assert(isempty(fig.UserData.webUI));
press(fig.UserData.workflowButtons(2));setValue(fig.UserData.controls.meshBackend,'engineering');
x=fig.UserData.workflowState;assert(isempty(x.mesh)&&~x.steps(3).available&&~isempty(x.orientedMesh));
press(fig.UserData.workflowButtons(1));press(fig.UserData.controls.flip);
assert(isempty(fig.UserData.workflowState.orientedMesh));setValue(fig.UserData.controls.scale,.001);
assert(isempty(fig.UserData.workflowState.rawMesh)&&~fig.UserData.workflowState.steps(2).available);
result=struct('originalFaces',n,'allStepsCompleted',true,'revision',double(s.revision),'solveId',double(s.solveId), ...
    'originalNumberingPreserved',true,'boundaryInvalidated',true,'meshInvalidated',true,'orientationInvalidated',true, ...
    'scaleInvalidated',true,'powerW',report.wholeSurfaceHeatPowerW,'renderedReceipt',rendered);
writeText(fullfile(out,'workflow.json'),jsonencode(result));
end
function press(c)
c.ButtonPushedFcn(c,[]);drawnow;
end
function setValue(c,value)
c.Value=value;c.ValueChangedFcn(c,[]);drawnow;
end
function writeText(path,text)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',text);
end
function deleteIfOpen(fig)
if isgraphics(fig),delete(fig);end
end
