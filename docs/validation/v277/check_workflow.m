function check_workflow(out)
% 在五步主窗口导入同一 STL，由新网格从冷态求解，再验证切换分辨率的失效行为。
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
cfg=createCaseConfig();cfg.output.verbose=false;cfg.geometry.scale=1;
cfg.geometry.stlFile=fullfile(root,'cases','CM011-alternate-SI.stl');
fig=rocketAeroHeatWorkbench(cfg);cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
fig.Visible='off';press(fig.UserData.controls.load);press(fig.UserData.controls.confirm);
choose(fig.UserData.controls.meshBackend,'cm011Imported');
assert(fig.UserData.controls.cm011HalfSectors.Value==12);
choose(fig.UserData.controls.cm011HalfSectors,24);press(fig.UserData.controls.generate);
s=fig.UserData.workflowState;
assert(~isempty(s.mesh)&&s.mesh.cm011Prepared.halfSectors==24,fig.UserData.status.Text);
p=s.mesh.cm011Prepared;saveText(out,'mesh_folder.txt',p.directory);
runCM011Python({fullfile(root,'docs','validation','v277','check_inputs.py'),p.directory,fullfile(out,'inputs.json')},fullfile(out,'preflight.log'));
press(fig.UserData.workflowButtons(3));choose(fig.UserData.controls.cm011Condition,'inclined');
press(fig.UserData.controls.confirm);press(fig.UserData.controls.run);
s=fig.UserData.workflowState;assert(strcmp(s.status,'solved')&&s.currentStep==5,fig.UserData.status.Text);
report=fig.UserData.caseReport;
assert(report.fullSurfaceComputed&&report.native.acceptedNumericalWall&&~report.reusedNativeRun);
saveText(out,'run_folder.txt',report.outputDirectory);
ui=fig.UserData.webUI;d=ui.html.UserData.payload;
assert(isequal(d.faces,s.rawMesh.F-1)&&d.researchIdentity.halfSectors==24);
started=tic;while toc(started)<16&&isempty(ui.html.UserData.rendered),drawnow;pause(.1);end
rendered=ui.html.UserData.rendered;
assert(rendered.faces==10080&&rendered.visiblePixelSamples>0&&rendered.webglError==0);
id=1701;ui.html.Data=struct('kind','selection','face',id);ui.html.DataChangedFcn(ui.html,[]);
assert(ui.html.UserData.selectedFace==id);
press(fig.UserData.workflowButtons(2));choose(fig.UserData.controls.cm011HalfSectors,12);
s=fig.UserData.workflowState;
assert(isempty(s.mesh)&&isempty(s.result)&&~s.boundaryConfirmed&&isempty(fig.UserData.webUI));
saveText(out,'workflow.json',jsonencode(struct('rendered',rendered,'pickedByEvent',id,'mouseDriven',false, ...
    'halfSectors',24,'sourceSHA256',p.sourceSHA256,'resetOnGridChange',true,'powerW',report.wholeSurfaceHeatPowerW)));
end
function press(control)
control.ButtonPushedFcn(control,[]);drawnow;
end
function choose(control,value)
control.Value=value;control.ValueChangedFcn(control,[]);drawnow;
end
function saveText(out,name,value)
fid=fopen(fullfile(out,name),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',value);
end
