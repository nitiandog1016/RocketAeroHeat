function testSU2MeshModeUI(out,geometrySource)
%TESTSU2MESHM0DEUI 核对主工作台网格选项、真实生成与下游失效，不重跑流场。
% 模块一：新工况默认法向层；缺字段的历史工况保留射线层且不偷偷补字段。
% 模块二：生成结果与独立研究网格参数一致，切换模式清除体网格和边界确认。
assert(~isfolder(out));mkdir(out);evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
[options,mode]=su2WorkbenchMeshOptions(cfg);assert(strcmp(mode,'normalLayers') && options.radialLayers==40);
legacy=cfg;legacy.su2=rmfield(legacy.su2,'volumeMode');
[old,mode]=su2WorkbenchMeshOptions(legacy);assert(strcmp(mode,'starRay') && old.radialLayers==8);
bad=cfg;bad.su2.volumeMode='unknown';refused=false;
try,su2WorkbenchMeshOptions(bad);catch ex,refused=strcmp(ex.identifier,'RocketAeroHeat:SU2:VolumeMode');end
assert(refused);
fig=rocketAeroHeatWorkbench(cfg);guard=onCleanup(@()closeIfPresent(fig));fig.UserData.actions.loadNative();
assert(strcmp(fig.UserData.controls.volumeMode.Value,'normalLayers'));
fig.UserData.actions.generateMesh();s=fig.UserData.workflowState;v=fig.UserData.volumeInfo;
assert(s.steps(2).completed && v.cells==58560 && v.nodes==59122);
saved=load(geometrySource,'volume','surface','options');
assert(isequaln(v.construction,saved.volume.meshConstruction) && isequaln(options,saved.options));
assert(isequal(s.mesh.V,saved.surface.V) && isequal(s.mesh.F,saved.surface.F));
assert(contains(fig.UserData.controls.volumeSummary.Text,'58560'));
fig.UserData.actions.selectStep(3);fig.UserData.actions.confirmBoundary();assert(fig.UserData.workflowState.boundaryConfirmed);
fig.UserData.actions.selectStep(2);before=fig.UserData.workflowState;
c=fig.UserData.controls.volumeMode;c.Value='starRay';c.ValueChangedFcn(c,[]);
after=fig.UserData.workflowState;
assert(after.revision>before.revision && isempty(after.mesh) && isempty(fig.UserData.volumeInfo));
assert(~after.boundaryConfirmed && isempty(after.researchResult) && ~after.steps(4).available && ~after.steps(5).available);
assert(isequaln(after.rawMesh,before.rawMesh));
% 仅切页不更改输入，也不生成新网格。
revision=after.revision;fig.UserData.actions.selectStep(1);fig.UserData.actions.selectStep(2);
assert(fig.UserData.workflowState.revision==revision && isempty(fig.UserData.volumeInfo));close(fig);clear guard;
fig=rocketAeroHeatWorkbench(legacy);guard=onCleanup(@()closeIfPresent(fig));fig.UserData.actions.loadNative();
assert(strcmp(fig.UserData.controls.volumeMode.Value,'starRay') && ~isfield(fig.UserData.workflowState.cfg.su2,'volumeMode'));
report=struct('newDefaultNormal',true,'oldSnapshotStarRayPreserved',true,'actualGeneratedCells',v.cells, ...
    'meshMatchesIndependentStudy',true,'modeChangeInvalidatesDownstream',true,'noSolverRun',true,'physicsQualified',false);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);fileGuard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Main workbench mesh-mode generation/legacy fallback/invalidation passed; no flow recomputation.\n');
end
function closeIfPresent(fig)
if isgraphics(fig),close(fig);end
end
