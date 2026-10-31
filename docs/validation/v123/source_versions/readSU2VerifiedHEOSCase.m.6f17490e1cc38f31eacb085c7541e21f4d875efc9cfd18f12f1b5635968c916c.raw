function result=readSU2VerifiedHEOSCase(folder,state,volume)
%READSU2VERIFIEDHEOSCASE 在主工作台只读查看已有HEOS场，不伪造旧UI输入快照。
% 模块一：核对真实成功终态、原配置/网格/清单散列及当前物理边界。
% 模块二：原文件中的节点/墙面与当前原面投影一致；显示源场自己的实际网格。
% 模块三：审核具名迭代和独立EOS；预算可不同，绝不据此建立续算身份。
assert(~state.busy && state.boundaryConfirmed && ~isempty(volume),'请先完成体网格并确认边界。');
assert(isfield(state.cfg.su2,'thermodynamicModel') && strcmp(state.cfg.su2.thermodynamicModel,'HEOS-Air'), ...
    '请先明确选择HEOS Air及对应定比热源场。');
folder=char(java.io.File(folder).getCanonicalPath());root=fileparts(mfilename('fullpath'));
assert(startsWith(lower(folder),lower([fullfile(root,'outputs') filesep])),'请选择本项目输出目录。');
raw=fileread(fullfile(folder,'engine_process.json'));receipt=jsondecode(raw);
assert(strcmp(receipt.status,'terminal') && receipt.exitCode==0 && ~receipt.timedOut && ~receipt.cancelled,'原引擎未成功结束。');
assert(strcmpi(folder,char(java.io.File(receipt.request.workDirectory).getCanonicalPath())),'必须从原运行目录查看。');
assert(contains(receipt.runtimeFlavor,'rah-coolprop'),'所选场没有实际HEOS运行身份。');
for name={'case.cfg','mesh.su2','case_manifest.json','coolprop_initialization.json'}
    file=fullfile(folder,name{1});token=regexp(raw,[regexptranslate('escape',jsonencode(file)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
    assert(numel(token)==1 && strcmp(sha(file),token{1}{1}),'原输入身份缺失或改变：%s',name{1});
end
manifest=jsondecode(fileread(fullfile(folder,'case_manifest.json')));f=manifest.flight;init=manifest.initialization;
assert(strcmp(manifest.thermodynamicModel,'CoolProp-7.2.0-HEOS-Air') && strcmp(manifest.thermalBoundary.mode,'isothermal'));
assert(~isfield(state.cfg.su2,'wallThermalMode') || strcmp(state.cfg.su2.wallThermalMode,'isothermal'),'当前热边界不同。');
actual=[f.referenceMach,f.alpha,f.beta,f.temperature,f.pressure,f.wallTemperature];
expected=[state.cfg.flight.Mach,state.cfg.flight.alpha,state.cfg.flight.beta,state.freestream.temperature,state.freestream.pressure,state.cfg.flight.Twall];
assert(all(abs(actual-expected)<=1e-11*max(abs(expected),1)),'当前来流/攻角/给定Tw与原场不一致。');
assert(strcmpi(char(java.io.File(init.sourceDirectory).getCanonicalPath()), ...
    char(java.io.File(state.cfg.su2.heosSource).getCanonicalPath())),'当前选定的物性转换源与原场不同。');
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
pipeline=fullfile(folder,'workflow','pipeline_process.json');
if isfile(pipeline)
    p=jsondecode(fileread(pipeline));assert(strcmp(p.status,'terminal') && p.exitCode==0 && ~p.cancelled && ~p.timedOut);
    assert(strcmp(p.nativeReceiptSHA256,sha(fullfile(folder,'engine_process.json'))) && ...
        strcmp(p.auditSHA256,sha(fullfile(folder,'thermodynamics','audit.json'))),'后台流程与EOS记录不匹配。');
end
projection=buildSU2WallProjection(volume,state.mesh);
result=mapSU2WallResult(readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk')),projection);
result.V=state.mesh.V;result.F=state.mesh.F;result.receipt=receipt;result.flight=f;result.resultDirectory=folder;
result=attachSU2ThermalBoundary(result,manifest.thermalBoundary);
result.convergenceAudit=auditSU2ResultConvergence(folder,result);
assert(result.convergenceAudit.checkpointIdentityConfirmed,result.convergenceAudit.reason);
result.verifiedCaseReview=struct('sourceDirectory',folder,'manifestSHA256',sha(fullfile(folder,'case_manifest.json')), ...
    'meshSHA256',sha(fullfile(folder,'mesh.su2')),'sourceNodeCount',manifest.nodes,'sourceCellCount',manifest.cells, ...
    'outputSettingsIgnoredForReadOnlyView',true,'currentVolumeUsedOnlyForSurfaceProjection',true, ...
    'historicalUIInputSnapshotCreated',false,'solverRerun',false,'physicsQualified',false);
if isfile(pipeline),result.workflowReceipt=p;end
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
