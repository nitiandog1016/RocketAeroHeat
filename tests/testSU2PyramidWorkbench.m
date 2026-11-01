function report=testSU2PyramidWorkbench(out)
%TESTSU2PYRAMIDWORKBENCH 在主五步骤核对实际混合单元说明，失败候选不能接入。
% 模块一：只读小型已完成接口网格，原壁面/质量状态与真实计数保持一致。
% 模块二：拒绝真实失败网格和不相容初猜；不启动流场，不重复旧全量UI验收。
assert(~isfolder(out));mkdir(out);d=load('outputs/v135_validation/fixture/mesh_input.mat');
s=struct('V',d.V,'F',d.F);c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;cfg.output.verbose=false;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',2,'outputEvery',2,'maximumSeconds',120,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));
fig.UserData.actions.useSurfaceForSU2(s);fig.UserData.actions.readHybridMesh('outputs/v135_validation/fixture');
info=fig.UserData.volumeInfo;method=fig.UserData.controls.meshMethod.Text;
assert(info.nodes==38 && info.cells==92 && info.construction.pyramids==8 && info.construction.noseRetainedLayers==1);
assert(contains(method,'16棱柱 / 8金字塔 / 68四面体') && contains(method,'鼻部保留1层') && contains(method,'未获物理资格'));
assert(strcmp(fig.UserData.controls.generate.Text,'重新生成常规体网格'));
waitGeometry(fig);payload=fig.UserData.geometryView.html.UserData.payload;
assert(isequal(fig.UserData.workflowState.mesh.V,s.V) && isequal(fig.UserData.workflowState.mesh.F,s.F));
rejected=false;try,fig.UserData.actions.readHybridMesh('outputs/v135_validation/actual');catch,rejected=true;end
assert(rejected && fig.UserData.volumeInfo.cells==92 && isempty(fig.UserData.workflowState.researchResult));
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed);
fig.UserData.actions.selectStep(3);
controls=fig.UserData.controls;controls.initialMode.Value='wallEnthalpy';controls.initialMode.ValueChangedFcn(controls.initialMode,[]);
assert(~fig.UserData.workflowState.boundaryConfirmed && contains(fig.UserData.status.Text,'尚未适配该初猜'));
assert(isempty(fig.UserData.engineJob));
save(fullfile(out,'preview.mat'),'info','payload','method','-v7');
report=struct('originalFaces',8,'nodes',38,'cells',92,'pyramids',8,'noseRetainedLayers',1, ...
    'mainWorkbenchWebglPassed',true,'exactCountsAndLayerText',true,'failedActualCandidateRejected',rejected, ...
    'incompatibleInitializerRejected',true,'sourceSurfaceUnchanged',true,'nativeMouseOrScreenshotTest',false, ...
    'nativeFlowRuns',0,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
function waitGeometry(fig)
t=tic;while toc(t)<40
    drawnow;data=fig.UserData.geometryView.html.UserData;
    if ~isempty(data.rendered) && data.rendered.revision==data.revision,break;end
    pause(.1);
end
assert(~isempty(data.rendered) && data.rendered.revision==data.revision && data.rendered.faces==8 && data.rendered.webglError==0);
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
