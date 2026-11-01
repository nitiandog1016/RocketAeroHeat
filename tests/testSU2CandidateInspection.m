function report=testSU2CandidateInspection(out)
%TESTSU2CANDIDATEINSPECTION 在主五步骤查看真实失败候选，验证不授予计算资格。
% 模块一：Delaunay/HXT同原面切换、WebGL数值/算法和探针保持真实来源。
% 模块二：拒绝原拓扑失败候选；返回原表面和上游变化不遗留候选身份。
assert(~isfolder(out));mkdir(out);d=load('outputs/v118_validation/actual01/surface.mat');
c=load('outputs/v56_validation/state_validity/state_validity_actual_rocket.mat','cfg');cfg=c.cfg;cfg.output.verbose=false;
cfg.su2=struct('volumeMode','hybridTriangles','iterations',2,'outputEvery',2,'maximumSeconds',120,'pythonExecutable','D:/anacondapython/python.exe');
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));fig.UserData.actions.useSurfaceForSU2(d.s);
before=fig.UserData.workflowState;folders={'outputs/v136_validation/explicit','outputs/v136_validation/hxt'};
counts=[205 475];algorithms={'Delaunay','HXT'};infos=cell(2,1);payloads=cell(2,1);picked=cell(2,1);
for j=1:2
    fig.UserData.actions.inspectMeshCandidate(folders{j});waitGeometry(fig);
    info=fig.UserData.meshQuality;payload=fig.UserData.geometryView.html.UserData.payload;
    assert(info.lowQualityCells==counts(j) && strcmp(info.exteriorAlgorithm,algorithms{j}) && ~info.acceptedForCalculation);
    assert(contains(payload.caption,'失败候选') && contains(payload.caption,'未接入计算') && isequal(payload.faceField.range,[0 .1]));
    assert(isequaln(before,fig.UserData.workflowState) && isempty(fig.UserData.volumeInfo) && isempty(fig.UserData.engineJob));
    assert(isempty(fig.UserData.workflowState.mesh) && ~fig.UserData.workflowState.steps(4).available);
    message=struct('kind','geometry-face-probe','revision',payload.revision,'faceIndex',0);
    html=fig.UserData.geometryView.html;html.Data=message;html.DataChangedFcn(html,[]);
    picked{j}=fig.UserData.meshQualityProbe;assert(picked{j}.cellID==info.displaySourceCellIDs(1));
    infos{j}=info;payloads{j}=payload;
end
rejected=false;try,fig.UserData.actions.inspectMeshCandidate('outputs/v135_validation/actual');catch,rejected=true;end
assert(rejected && fig.UserData.meshQuality.lowQualityCells==475 && isequaln(before,fig.UserData.workflowState));
button=fig.UserData.controls.meshQuality;button.ButtonPushedFcn(button,[]);waitGeometry(fig);
assert(isempty(fig.UserData.meshQuality) && size(fig.UserData.geometryView.html.UserData.payload.faces,1)==19752);
% 重新显示已缓存原始候选，再改上游；无需重生成或原生流场。
fig.UserData.actions.inspectMeshCandidate(folders{1});waitGeometry(fig);
ctrl=fig.UserData.controls.volumeMode;ctrl.Value='normalLayers';ctrl.ValueChangedFcn(ctrl,[]);waitGeometry(fig);
assert(isempty(fig.UserData.meshQuality) && ~isfield(fig.UserData.geometryView.html.UserData.payload,'faceField'));
save(fullfile(out,'candidate_previews.mat'),'infos','payloads','picked','-v7');
report=struct('originalFaces',19752,'candidateBadTets',counts,'actualAlgorithms',{algorithms}, ...
    'mainWorkbenchWebglPassed',true,'calculationInputsUnchanged',true,'noSolverMeshAccepted',true,'bridgeProbePassed',true, ...
    'topologyFailureRejected',rejected,'returnToOriginalPassed',true,'upstreamInvalidationPassed',true, ...
    'nativeMouseOrScreenshotTest',false,'nativeFlowRuns',0,'newQualifiedArea',0);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);g=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
end
function waitGeometry(fig)
t=tic;while toc(t)<40
    drawnow;s=fig.UserData.geometryView.html.UserData;
    if ~isempty(s.rendered) && s.rendered.revision==s.revision,break;end
    pause(.1);
end
assert(~isempty(s.rendered) && s.rendered.revision==s.revision && s.rendered.faces==size(s.payload.faces,1) && s.rendered.webglError==0);
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
