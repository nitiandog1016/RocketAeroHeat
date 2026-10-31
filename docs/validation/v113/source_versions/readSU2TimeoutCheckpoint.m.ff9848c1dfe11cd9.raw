function result=readSU2TimeoutCheckpoint(folder,state,volume)
%READSU2TIMEOUTCHECKPOINT 显式重开本程序同工况超时检查点，仅供研究查看。
% 模块一：必须有原工作台输入快照和真实超时终态；不由缺文件猜测成功。
% 模块二：校验原快照/网格/配置SHA、当前输入及原投影，保护面元身份。
% 模块三：绑定保存迭代并核对HF/EOS，保留未收敛和未资格，不启动求解。
assert(~state.busy && state.boundaryConfirmed && ~isempty(volume),'请先生成网格并确认当前工况。');
folder=char(java.io.File(folder).getCanonicalPath());
raw=fileread(fullfile(folder,'engine_process.json'));receipt=jsondecode(raw);
assert(strcmp(receipt.status,'terminal') && receipt.timedOut && ~receipt.cancelled, ...
    '只接受真实超时任务；其他失败、取消或运行中的任务不能按本入口重开。');
assert(strcmpi(folder,char(java.io.File(receipt.request.workDirectory).getCanonicalPath())), ...
    '请从原运行目录重开，不能为复制文件伪造原进程身份。');
snapshot=fullfile(folder,'input_snapshot.mat');
for file={snapshot,fullfile(folder,receipt.request.meshFile),fullfile(folder,receipt.request.configFile)}
    token=regexp(raw,[regexptranslate('escape',jsonencode(file{1})) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
    assert(numel(token)==1 && strcmp(sha(file{1}),token{1}{1}),'原输入快照／网格／配置身份不完整或已改变。');
end
saved=load(snapshot,'frozen','projection');assert(isfield(saved,'projection'),'旧快照缺少原投影，不能按此入口推定映射。');
f=saved.frozen;
assert(isequaln(f.mesh,state.mesh) && isequaln(f.cfg,state.cfg) && isequaln(f.freestream,state.freestream), ...
    '原模型、边界或输出设置与当前已确认工况不一致。');
assert(isequal(saved.projection.sourceVertices,volume.vertices) && ...
    isequal(saved.projection.sourceCells,volume.bodyWall),'当前体网格与原投影不同。');
data=readSU2SurfaceResult(fullfile(folder,'surface_flow.vtk'));
result=mapSU2WallResult(data,saved.projection);result.V=f.mesh.V;result.F=f.mesh.F;
result.receipt=receipt;result.flight=f.cfg.flight;result.wallTemperatureInput=f.cfg.flight.Twall;
result.resultDirectory=folder;result.sourceRevision=state.revision;result.solveId=state.solveId;
result.convergenceAudit=auditSU2ResultConvergence(folder,result);
assert(result.convergenceAudit.checkpointIdentityConfirmed, ...
    ['无法确认完整检查点：' result.convergenceAudit.reason]);
result.explicitTimeoutReview=struct('inputSnapshotSHA256',sha(snapshot), ...
    'checkpointIteration',result.convergenceAudit.checkpointIteration,'solverRerun',false);
result.resultStatus='explicit-timeout-checkpoint; research-only';
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0,'缺少原输入文件。');guard=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
