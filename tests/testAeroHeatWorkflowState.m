function testAeroHeatWorkflowState(out)
%TESTAEROHEATWORKFLOWSTATE 验证五步骤集中状态、真实上游修改及过期快照保护。
% 模块一：按输入/朝向/网格/边界/显式计算/结果执行完整路径。
% 模块二：缩放、旋转、同面数新网格及边界编辑的下游失效；确认输入
% 壁温不冒充求解值，失败保留输入，浏览与幂等提交不自动重算。
% 模块三：通过实际经典锥柱内部外流/热流后端确认完整快照的场量、
% NaN 与算法来源逐字保留。本测试是软件流程回归，不是实验精度验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=5;
cfg.flight.beta=0;cfg.flight.Twall=300;
raw=struct('V',[0 0 0;1 0 0;0 1 0;0 0 1], ...
    'F',[1 3 2;1 2 4;1 4 3;2 3 4]);
s=aeroHeatWorkflowState('init',[],struct('cfg',cfg));
assert(s.currentStep==1 && isequal([s.steps.available],[true false false false false]));
assert(~any([s.steps.completed]) && isempty(s.result) && ~s.busy);
mustThrow(@()aeroHeatWorkflowState('selectStep',s,struct('step',5)), ...
    'RocketAeroHeat:WorkflowPrerequisite');
mustThrow(@()aeroHeatWorkflowState('beginSolve',s), ...
    'RocketAeroHeat:WorkflowPrerequisite');

%% 模块一：真实几何和已确认朝向是网格页的前置条件。
s=aeroHeatWorkflowState('setModel',s,struct('rawMesh',raw));
assert(~s.steps(1).completed && ~s.steps(2).available);
s=aeroHeatWorkflowState('setOrientation',s,struct('orientedMesh',raw));
assert(s.steps(1).completed && s.steps(2).available);
s=aeroHeatWorkflowState('setMesh',s,struct('mesh',raw));
assert(s.steps(2).completed && s.steps(3).available && ~s.steps(4).available);
s=aeroHeatWorkflowState('confirmBoundary',s);
assert(s.boundaryConfirmed && s.steps(4).available && ~s.steps(4).completed);
ready=s;
s=aeroHeatWorkflowState('beginSolve',s);
token=solvePayload(s);snapshot=s.frozenInputs;
assert(s.busy && strcmp(s.status,'solving') && s.currentStep==4);
assert(isequaln(snapshot.cfg,cfg) && isequal(snapshot.mesh,raw));
assert(~isfield(snapshot.mesh,'pressure') && ~isfield(snapshot.mesh,'heating'));
mustThrow(@()aeroHeatWorkflowState('beginSolve',s), ...
    'RocketAeroHeat:WorkflowPrerequisite');
result=raw;result.pressure=struct('staticPressure',[120;NaN;80;0], ...
    'Cp',[1;NaN;0;-1],'pressureIncrement',[20;NaN;0;-100], ...
    'fieldStatus',uint8([1;0;2;1]),'missingPressureFace',[false;true;false;false], ...
    'model','软件合同样本','isPrescribedBase',[false;false;true;false]);
result.heating=struct('heatFlux',[10;NaN;-3;0], ...
    'eligibleFace',[true;false;true;true], ...
    'faceMethodID',uint8([1;0;2;1]),'wallTemperatureByFace',ones(4,1)*300);
token.result=result;s=aeroHeatWorkflowState('completeSolve',s,token);
assert(~s.busy && s.steps(5).available && all([s.steps.completed]));
assert(isequaln(s.result,result) && isequal(s.resultSourceRevision,s.revision));
assert(~isfield(s.result,'wallThermal') && s.cfg.flight.Twall==300);
completed=s;

%% 模块二：浏览页面和未改变输入的应用操作不重算、不失效。
for i=1:5,s=aeroHeatWorkflowState('selectStep',s,struct('step',i));end
assert(isequaln(s.result,completed.result) && s.revision==completed.revision);
assert(s.solveId==completed.solveId && ~s.busy);
s=aeroHeatWorkflowState('setBoundary',s,struct('cfg',s.cfg));
s=aeroHeatWorkflowState('setOrientation',s,struct('orientedMesh',raw));
s=aeroHeatWorkflowState('setMesh',s,struct('mesh',raw));
s=aeroHeatWorkflowState('setMeshOptions',s,struct('cfg',s.cfg));
assert(isequaln(s.result,completed.result) && s.revision==completed.revision);

%% 模块三：实际朝向和尺度改变都清除网格、来流、确认状态及快照。
rotated=raw;rotated.V=raw.V*[0 -1 0;1 0 0;0 0 1];
changed=aeroHeatWorkflowState('setOrientation',completed,struct('orientedMesh',rotated));
assert(changed.revision>completed.revision && isempty(changed.mesh));
assertCleared(changed);assert(~changed.steps(3).available);
scaled=raw;scaled.V=2*raw.V;
changed=aeroHeatWorkflowState('setModel',completed, ...
    struct('rawMesh',scaled,'orientedMesh',scaled));
assert(isempty(changed.mesh) && changed.steps(2).available);assertCleared(changed);
newCfg=cfg;newCfg.geometry.scale=2*cfg.geometry.scale;
changed=aeroHeatWorkflowState('invalidateModel',completed,struct('cfg',newCfg));
assert(isempty(changed.rawMesh) && isempty(changed.orientedMesh) && isempty(changed.mesh));
assert(changed.cfg.geometry.scale==newCfg.geometry.scale);assertCleared(changed);
changed=aeroHeatWorkflowState('invalidateOrientation',completed);
assert(isequal(changed.rawMesh,raw) && isempty(changed.orientedMesh) && isempty(changed.mesh));
assert(~changed.steps(2).available);assertCleared(changed);

%% 模块四：同面数新网格也必须清除旧逐面元壁温绑定，保留统一壁温。
byFaceCfg=cfg;byFaceCfg.heating.wallTemperatureByFace=[300;310;320;330];
faceState=aeroHeatWorkflowState('setBoundary',completed,struct('cfg',byFaceCfg));
assert(isequal(faceState.cfg.heating.wallTemperatureByFace,[300;310;320;330]));
assertCleared(faceState);
faceState=aeroHeatWorkflowState('confirmBoundary',faceState);
faceState=aeroHeatWorkflowState('beginSolve',faceState);
faceToken=solvePayload(faceState);faceToken.result=result;
faceState=aeroHeatWorkflowState('completeSolve',faceState,faceToken);
newMesh=raw;newMesh.V(2,1)=1.25;
changed=aeroHeatWorkflowState('setMesh',faceState,struct('mesh',newMesh));
assert(size(changed.mesh.F,1)==size(raw.F,1) && isequal(changed.mesh.V,newMesh.V));
assert(~isfield(changed.cfg.heating,'wallTemperatureByFace'));
assert(changed.cfg.flight.Twall==300 && contains(changed.lastInvalidation,'逐面元壁温'));
assertCleared(changed);assert(changed.steps(3).available && ~changed.steps(4).available);
meshCfg=faceState.cfg;meshCfg.mesh.workflowTestParameter=12;
changed=aeroHeatWorkflowState('setMeshOptions',faceState,struct('cfg',meshCfg));
assert(isempty(changed.mesh) && isequal(changed.orientedMesh,raw));
assert(~isfield(changed.cfg.heating,'wallTemperatureByFace'));assertCleared(changed);
changed=aeroHeatWorkflowState('invalidateOrientation',faceState);
assert(~isfield(changed.cfg.heating,'wallTemperatureByFace'));assertCleared(changed);
changed=aeroHeatWorkflowState('invalidateModel',faceState);
assert(~isfield(changed.cfg.heating,'wallTemperatureByFace'));assertCleared(changed);

%% 模块五：不同攻角或壁温使旧结果失效，单纯查看不将新文字配旧云图。
newCfg=cfg;newCfg.flight.alpha=20;newCfg.flight.Twall=350;
changed=aeroHeatWorkflowState('setBoundary',completed,struct('cfg',newCfg));
assert(changed.cfg.flight.alpha==20 && changed.cfg.flight.Twall==350);
assert(~changed.steps(5).available && ~changed.boundaryConfirmed);assertCleared(changed);
mustThrow(@()aeroHeatWorkflowState('selectStep',changed,struct('step',5)), ...
    'RocketAeroHeat:WorkflowPrerequisite');
bad=changed;bad.cfg.flight.Twall=NaN;
mustThrow(@()aeroHeatWorkflowState('confirmBoundary',bad), ...
    'RocketAeroHeat:WorkflowBoundary');

%% 模块六：求解失败保留当前输入、明确报错，并可同版本重试。
solving=aeroHeatWorkflowState('beginSolve',ready);oldToken=solvePayload(solving);
failedToken=oldToken;failedToken.error=MException('RocketAeroHeat:InjectedFailure','后端测试失败');
failed=aeroHeatWorkflowState('failSolve',solving,failedToken);
assert(strcmp(failed.status,'failed') && strcmp(failed.steps(4).status,'failed'));
assert(~failed.busy && isempty(failed.result) && failed.boundaryConfirmed);
assert(isequal(failed.mesh,ready.mesh) && isequaln(failed.cfg,ready.cfg));
assert(strcmp(failed.error.identifier,'RocketAeroHeat:InjectedFailure'));
retry=aeroHeatWorkflowState('beginSolve',failed);
assert(retry.revision==solving.revision && retry.solveId>solving.solveId);
oldToken.result=result;
mustThrow(@()aeroHeatWorkflowState('completeSolve',retry,oldToken), ...
    'RocketAeroHeat:WorkflowStaleResult');
mustThrow(@()aeroHeatWorkflowState('failSolve',retry,failedToken), ...
    'RocketAeroHeat:WorkflowStaleResult');
assert(retry.busy && isempty(retry.result));
freshToken=solvePayload(retry);freshToken.result=result;
done=aeroHeatWorkflowState('completeSolve',retry,freshToken);assert(done.steps(5).available);

%% 模块七：计算中真实输入修改取消旧版本，旧返回不能覆盖新工况。
solving=aeroHeatWorkflowState('beginSolve',ready);oldToken=solvePayload(solving);
changed=aeroHeatWorkflowState('setBoundary',solving,struct('cfg',newCfg));
assertCleared(changed);assert(~changed.busy);
assert(solving.frozenInputs.cfg.flight.alpha==5 && changed.cfg.flight.alpha==20);
oldToken.result=result;
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,oldToken), ...
    'RocketAeroHeat:WorkflowStaleResult');
changed=aeroHeatWorkflowState('confirmBoundary',changed);
changed=aeroHeatWorkflowState('beginSolve',changed);
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,oldToken), ...
    'RocketAeroHeat:WorkflowStaleResult');
wrong=solvePayload(changed);wrong.result=scaled;
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultGeometry');
wrong.result=raw;
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultFields');
wrong.result=result;wrong.result.heating.heatFlux=wrong.result.heating.heatFlux(1:3);
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultFields');
wrong.result=result;wrong.result.pressure.staticPressure(1)=Inf;
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultFields');
wrong.result=result;wrong.result.pressure.fieldStatus(2)=uint8(1);
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultFields');
wrong.result=result;wrong.result.pressure.staticPressure(2)=0;
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultFields');
wrong.result=result;wrong.result.heating.eligibleFace(2)=true;
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultFields');
wrong.result=result;wrong.result.heating.eligibleFace(1)=false;
mustThrow(@()aeroHeatWorkflowState('completeSolve',changed,wrong), ...
    'RocketAeroHeat:WorkflowResultFields');
assert(changed.busy && isempty(changed.result));
wrongCfg=ready.cfg;wrongCfg.geometry.scale=2*wrongCfg.geometry.scale;
mustThrow(@()aeroHeatWorkflowState('setBoundary',ready,struct('cfg',wrongCfg)), ...
    'RocketAeroHeat:WorkflowInputHierarchy');
mustThrow(@()aeroHeatWorkflowState('setMeshOptions',ready,struct('cfg',wrongCfg)), ...
    'RocketAeroHeat:WorkflowInputHierarchy');
wrongCfg=ready.cfg;wrongCfg.mesh.workflowTestParameter=20;
mustThrow(@()aeroHeatWorkflowState('setBoundary',ready,struct('cfg',wrongCfg)), ...
    'RocketAeroHeat:WorkflowInputHierarchy');
wrongCfg=ready.cfg;wrongCfg.remesh.workflowTestParameter=20;
mustThrow(@()aeroHeatWorkflowState('setBoundary',ready,struct('cfg',wrongCfg)), ...
    'RocketAeroHeat:WorkflowInputHierarchy');

% 完成的是后端执行，不表示全面积已解；全NaN快照不变成物理零值。
allMissing=result;allMissing.pressure.staticPressure(:)=NaN;
allMissing.pressure.Cp(:)=NaN;allMissing.pressure.pressureIncrement(:)=NaN;
allMissing.pressure.fieldStatus(:)=uint8(0);allMissing.pressure.missingPressureFace(:)=true;
allMissing.heating.heatFlux(:)=NaN;allMissing.heating.eligibleFace(:)=false;
missingState=aeroHeatWorkflowState('beginSolve',ready);missingToken=solvePayload(missingState);
missingToken.result=allMissing;
missingState=aeroHeatWorkflowState('completeSolve',missingState,missingToken);
assert(all(isnan(missingState.result.pressure.staticPressure)));
assert(all(isnan(missingState.result.heating.heatFlux)) && missingState.steps(5).available);

%% 模块八：几何输入严格保护；从旧计算网格重新导入不带旧场量。
bad=raw;bad.F(1,1)=0;
mustThrow(@()aeroHeatWorkflowState('setModel',completed,struct('rawMesh',bad)), ...
    'RocketAeroHeat:WorkflowMesh');
bad=raw;bad.V(2,:)=bad.V(1,:);
mustThrow(@()aeroHeatWorkflowState('setModel',completed,struct('rawMesh',bad)), ...
    'RocketAeroHeat:WorkflowMesh');
mustThrow(@()aeroHeatWorkflowState('setModel',completed, ...
    struct('rawMesh',raw,'modelSource',struct('units','mm'))), ...
    'RocketAeroHeat:WorkflowUnits');
imported=aeroHeatWorkflowState('init',[],struct('cfg',cfg));
imported=aeroHeatWorkflowState('setModel',imported, ...
    struct('rawMesh',result,'orientedMesh',result));
imported=aeroHeatWorkflowState('setMesh',imported,struct('mesh',result));
assert(~isfield(imported.mesh,'pressure') && ~isfield(imported.mesh,'heating'));
assert(isempty(imported.result) && ~imported.boundaryConfirmed);

%% 模块九：真实内部后端全快照沿同一状态路径接收，未解 NaN 不补零。
mesh=buildCanonicalBody('coneCylinder',struct('resolution',12));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
realState=aeroHeatWorkflowState('init',[],struct('cfg',cfg));
realState=aeroHeatWorkflowState('setModel',realState, ...
    struct('rawMesh',mesh,'orientedMesh',mesh));
realState=aeroHeatWorkflowState('setMesh',realState,struct('mesh',mesh));
fs=calcFreestream(cfg.flight);
realState=aeroHeatWorkflowState('confirmBoundary',realState,struct('freestream',fs));
realState=aeroHeatWorkflowState('beginSolve',realState);
frozen=realState.frozenInputs;
realResult=calcFlightStationState(frozen.mesh,frozen.cfg.flight,frozen.cfg);
token=solvePayload(realState);token.result=realResult;token.freestream=fs;
realState=aeroHeatWorkflowState('completeSolve',realState,token);
assert(isequaln(realState.result,realResult) && realState.steps(5).available);
assert(any(isfinite(realState.result.heating.heatFlux)));
assert(any(isnan(realState.result.heating.heatFlux)));
assert(isequaln(realState.result.heating.faceMethodID,realResult.heating.faceMethodID));
assert(~isfield(realState.result,'wallThermal'));
if nargin>0
    if ~exist(out,'dir'),mkdir(out);end
    save(fullfile(out,'workflow_state_regression.mat'),'completed','failed','changed', ...
        'realState','realResult');
end
fprintf(['Workflow state regression passed: five actual steps; unchanged viewing; ' ...
    'rotation/scale/same-face-count mesh/boundary invalidation; face-wall binding; ' ...
    'failure/retry/stale result, result fields and configuration hierarchy guards; ' ...
    '%d-face real internal chain ' ...
    'preserved NaN and provenance.\n'],size(realResult.F,1));
end
function token=solvePayload(s)
token=struct('sourceRevision',s.solveInputRevision,'solveId',s.solveId);
end
function assertCleared(s)
assert(isempty(s.result) && isnan(s.resultSourceRevision) && isempty(s.freestream));
assert(~s.boundaryConfirmed && ~s.steps(5).available);
end
function mustThrow(fn,id)
try
    fn();
catch e
    assert(strcmp(e.identifier,id), ...
        '预期错误 %s，实际 %s：%s',id,e.identifier,e.message);
    return
end
error('RocketAeroHeat:TestMissingGuard','缺少预期守卫：%s。',id);
end
