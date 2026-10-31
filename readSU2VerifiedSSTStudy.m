function result=readSU2VerifiedSSTStudy(folder,state,volume)
%READSU2VERIFIEDSSTSTUDY 在主工作台只读查看SST数值对照，不伪造旧UI输入快照。
% 模块一：核对真实成功终态、原配置/网格/清单散列及当前物理边界。
% 模块二：原文件中的节点/墙面与当前原面投影一致；显示源场自己的实际网格。
% 模块三：审核具名迭代和独立EOS；预算可不同，绝不据此建立续算身份。
assert(~state.busy && state.boundaryConfirmed && ~isempty(volume),'请先完成体网格并确认边界。');
assert(isfield(state.cfg.su2,'flowModel') && strcmp(state.cfg.su2.flowModel,'SST-2003m'),'请先选择SST-2003m。');
assert(~isfield(state.cfg.su2,'thermodynamicModel') || strcmp(state.cfg.su2.thermodynamicModel,'idealGas'),'本数值对照仅支持定比热。');
folder=char(java.io.File(folder).getCanonicalPath());root=fileparts(mfilename('fullpath'));
assert(startsWith(lower(folder),lower([fullfile(root,'outputs') filesep])),'请选择本项目输出目录。');
raw=fileread(fullfile(folder,'engine_process.json'));receipt=jsondecode(raw);
assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut && ~receipt.cancelled,'原引擎未成功结束。');
assert(strcmpi(folder,char(java.io.File(receipt.request.workDirectory).getCanonicalPath())),'必须从原运行目录查看。');
assert(strcmp(receipt.runtimeFlavor,'official') && isfield(receipt.physicsModel,'numericalStudy'),'缺少实际数值对照身份。');
for name={'case.cfg','mesh.su2','case_manifest.json','sst_initialization.json','sst_contract.json'}
    file=fullfile(folder,name{1});token=regexp(raw,[regexptranslate('escape',jsonencode(file)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
    assert(numel(token)==1 && strcmp(sha(file),token{1}{1}),'原输入身份缺失或改变：%s',name{1});
end
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));f=manifest.flight;init=manifest.initialization;
assert(strcmp(manifest.thermalBoundary.mode,'isothermal') && strcmp(init.kind,'same-mesh-SST-CFL-cap-comparison'));
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
result.thermalStateAudit=readSU2SSTStateAudit(folder,receipt);
result.convergenceAudit=auditSU2ResultConvergence(folder,result);
assert(result.convergenceAudit.checkpointIdentityConfirmed,result.convergenceAudit.reason);
assert(isfield(result.convergenceAudit,'numericalStudy'),'未审核的数值方案不能进入查看器。');
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
