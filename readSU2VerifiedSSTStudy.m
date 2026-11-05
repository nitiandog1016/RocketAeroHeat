function result=readSU2VerifiedSSTStudy(folder,state,volume)
%READSU2VERIFIEDSSTSTUDY 在主工作台只读查看已审核的SST完整域结果。
% 模块一。按定比热研究/变比热同网格合同核对实际终态、来源和物理边界。
% 模块二。原文件节点及墙面与当前原面投影一致，显示实际计算网格。
% 模块三。核对具名迭代、独立EOS及温度身份，单独记录只读查看条件。
assert(~state.busy && state.boundaryConfirmed && ~isempty(volume),'请先完成体网格并确认边界。');
assert(isfield(state.cfg.su2,'flowModel') && strcmp(state.cfg.su2.flowModel,'SST-2003m'),'请先选择SST-2003m。');
model='idealGas';if isfield(state.cfg.su2,'thermodynamicModel'),model=state.cfg.su2.thermodynamicModel;end
assert(ismember(model,{'idealGas','idealAir'}),'本入口支持定比热和理想空气变比热SST。');
variableAir=strcmp(model,'idealAir');
folder=char(java.io.File(folder).getCanonicalPath());root=fileparts(mfilename('fullpath'));
assert(startsWith(lower(folder),lower([fullfile(root,'outputs') filesep])),'请选择本项目输出目录。');
raw=fileread(fullfile(folder,'engine_process.json'));receipt=jsondecode(raw);
assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut && ~receipt.cancelled,'原引擎未成功结束。');
assert(strcmpi(folder,char(java.io.File(receipt.request.workDirectory).getCanonicalPath())),'必须从原运行目录查看。');
isNumerical=isfield(receipt.physicsModel,'numericalStudy');isInitial=isfield(receipt.physicsModel,'initializationStudy');
isTrace=isfield(receipt.physicsModel,'updateTraceStudy');
inputs={'case.cfg','mesh.su2','case_manifest.json','sst_contract.json'};
if variableAir
    assert(strcmp(receipt.runtimeFlavor,'rah-ideal-air-phase-v2') && isfield(receipt,'thermodynamicModel'), ...
        '变比热结果需要实际同阶段运行合同。');
    inputs=[inputs,{'ideal_air_restart.json','ideal_air_contract.json','ideal_air_restart_preparation.json','initial_flow.csv'}];
else
    assert(ismember(receipt.runtimeFlavor,{'official','rah-bounded-conservative-v1','rah-update-trace-v2'}) && ...
        sum([isNumerical isInitial isTrace])==1,'缺少唯一实际数值/初场/更新观测身份。');
    inputs{end+1}='sst_initialization.json';
end
for name=inputs
    file=fullfile(folder,name{1});token=regexp(raw,[regexptranslate('escape',jsonencode(file)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
    assert(numel(token)==1 && strcmp(sha(file),token{1}{1}),'原输入身份缺失或改变：%s',name{1});
end
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));f=manifest.flight;init=manifest.initialization;
assert(strcmp(manifest.thermalBoundary.mode,'isothermal'));
if variableAir
    assert(strcmp(init.kind,'same-mesh-ideal-air-state-v1') && init.samePhysicsAndNumerics && ...
        init.localIterationStartsAt==0 && ~init.physicsQualified);
    assert(isequaln(manifest.idealAir,receipt.thermodynamicModel));
    assert(strcmp(manifest.idealAir.thermodynamicModel,'Lemmon2000-ideal-Air-CoolProp7.2.0-OTH-v1'));
    assert(norm(manifest.expectedVelocity(:)-state.freestream.flowVelocityVector(:))<2e-10, ...
        '当前变比热来流速度与实际原场不同。');
elseif isTrace
    assert(strcmp(init.kind,'same-mesh-SST-checkpoint-research') && init.diagnosticOnly && init.runtimeObservationOnly);
    assert(strcmp(receipt.runtimeFlavor,'rah-update-trace-v2'));
elseif isNumerical
    assert(ismember(init.kind,{'same-mesh-SST-CFL-cap-comparison','same-mesh-SST-bounded-conservative-comparison', ...
        'same-mesh-SST-gradient-comparison','same-mesh-SST-reconstruction-startup-comparison', ...
        'same-mesh-SST-flow-flux-comparison'}));
    isBounded=strcmp(init.kind,'same-mesh-SST-bounded-conservative-comparison');
    assert(isBounded==strcmp(receipt.runtimeFlavor,'rah-bounded-conservative-v1'),'方案与实际运行包身份不同。');
else,assert(ismember(init.kind,{'same-mesh-SST-omega-pocket-initialization-study', ...
        'same-mesh-SST-wall-column-initialization-study'}));assert(strcmp(receipt.runtimeFlavor,'official'));end
assert(isequaln(manifest.sst,receipt.physicsModel));
assert(~isfield(state.cfg.su2,'wallThermalMode') || strcmp(state.cfg.su2.wallThermalMode,'isothermal'),'当前热边界不同。');
actual=[f.Mach,f.alpha,f.beta,f.temperature,f.pressure,f.wallTemperature];
expected=[state.cfg.flight.Mach,state.cfg.flight.alpha,state.cfg.flight.beta,state.freestream.temperature,state.freestream.pressure,state.cfg.flight.Twall];
assert(all(abs(actual-expected)<=1e-11*max(abs(expected),1)),'当前来流/攻角/给定Tw与原场不一致。');
for field={'turbulenceIntensity','turbulentViscosityRatio','turbulentPrandtl'}
    assert(isfield(state.cfg.su2,field{1}),'请先确认显式湍流参数。');
end
assert(isequal([manifest.sst.intensity manifest.sst.viscosityRatio manifest.sst.Prt], ...
    [state.cfg.su2.turbulenceIntensity state.cfg.su2.turbulentViscosityRatio state.cfg.su2.turbulentPrandtl]));
% 从受实际进程绑定的网格读取节点与真实墙面；不依赖未锁定的projection.mat。
mesh=fileread(fullfile(folder,'mesh.su2'));node=regexp(mesh,'(?m)^NPOIN=\s*(\d+)\s*\r?\n','end','once');
marker=regexp(mesh,'(?m)^NMARK=','start','once');assert(~isempty(node) && marker>node);
points=sscanf(mesh(node+1:marker-1),'%f',[4 inf])';
assert(size(points,1)==manifest.nodes && isequal(points(:,4),(0:size(points,1)-1)'));
assert(isequal(points(:,1:3),volume.vertices),'原场节点与当前体网格节点不同；不能套当前投影。');
wallStart=regexp(mesh,'(?m)^MARKER_TAG=\s*wall\s*\r?\nMARKER_ELEMS=\s*\d+\s*\r?\n','end','once');
farStart=regexp(mesh,'(?m)^MARKER_TAG=\s*farfield','start','once');assert(~isempty(wallStart) && farStart>wallStart);
wall=sscanf(mesh(wallStart+1:farStart-1),'%d',[4 inf])';
assert(size(wall,1)==size(state.mesh.F,1) && all(wall(:,1)==5),'当前只读入口需要完整原三角墙面。');
expectedWall=vertcat(volume.wall{:});expectedWall=expectedWall(:,[1 3 2]);
assert(isequal(wall(:,2:4)+1,double(expectedWall)) && isequal(volume.sourceVertices,state.mesh.V) && ...
    isequal(volume.sourceFaces,state.mesh.F),'原墙面/顺序与当前原三角面不同。');
projection=buildSU2WallProjection(volume,state.mesh);
result=mapSU2WallResult(readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk')),projection);
result.V=state.mesh.V;result.F=state.mesh.F;result.receipt=receipt;result.flight=f;result.resultDirectory=folder;
result=attachSU2ThermalBoundary(result,manifest.thermalBoundary);
result.convergenceAudit=auditSU2ResultConvergence(folder,result);
assert(result.convergenceAudit.checkpointIdentityConfirmed,result.convergenceAudit.reason);
result.thermalStateAudit=result.convergenceAudit.thermalStateAudit;
if variableAir
    assert(result.thermalStateAudit.interfacePassed && ...
        result.thermalStateAudit.seriesSuccessfulUpdates==init.priorSuccessfulUpdates+result.thermalStateAudit.successfulUpdates);
    result.actualFreestreamSpeed=receipt.thermodynamicModel.freestreamSpeed;
else
    assert(isfield(result.convergenceAudit,'numericalStudy') || isfield(result.convergenceAudit,'initializationStudy') || isfield(result.convergenceAudit,'updateTrace'), ...
        '未审核的数值/初场方案不能进入查看器。');
end
% 旧源清单可能保留源任务输出设置；实际步数/频率仅从已锁定的当前配置读取。
configuration=fileread(fullfile(folder,'case.cfg'));runControls=struct;
for item={'ITER','OUTPUT_WRT_FREQ'}
    token=regexp(configuration,['(?m)^' item{1} '=\s*(\d+)'],'tokens','once');assert(~isempty(token));
    runControls.(item{1})=str2double(token{1});
end
result.inheritedManifestRunControls=struct('iterations',f.iterations,'outputEvery',f.outputEvery);
result.flight.iterations=runControls.ITER;result.flight.outputEvery=runControls.OUTPUT_WRT_FREQ;
result.actualRunControls=runControls;
result.verifiedCaseReview=struct('sourceDirectory',folder,'manifestSHA256',sha(fullfile(folder,'case_manifest.json')), ...
    'meshSHA256',sha(fullfile(folder,'mesh.su2')),'sourceNodeCount',manifest.nodes,'sourceCellCount',manifest.cells, ...
    'outputSettingsIgnoredForReadOnlyView',true,'currentVolumeUsedOnlyForSurfaceProjection',true, ...
    'historicalUIInputSnapshotCreated',false,'solverRerun',false,'physicsQualified',false);

end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
