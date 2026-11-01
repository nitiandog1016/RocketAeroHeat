function summary=readSU2UpdateTraceAudit(folder,receipt)
%READSU2UPDATETRACEAUDIT 绑定实际增量重放与全局节点审核，供科研工作台只读查看。
% 模块一：核对固定研究身份和一次更新，不能将观测工具作为普通求解器。
% 模块二：绑定真实回执与四个原生记录的散列，防止编号或阶段错配。
% 模块三：接口成功不代表残差收敛或区域热流取得物理资格。
assert(strcmp(receipt.runtimeFlavor,'rah-update-trace-v2'));
spec=receipt.physicsModel.updateTraceStudy;
assert(strcmp(spec.studyID,'v147-native-implicit-update-observation') && strcmp(spec.variant,'actual') && spec.enabled);
path=fullfile(fileparts(folder),'native_update_analysis.json');text=fileread(path);summary=jsondecode(text);
assert(summary.passed && summary.sourceLocalIteration==99 && summary.newLocalIteration==0);
assert(summary.nodes==335726 && summary.originalFaces==19752 && ~summary.physicsQualified && summary.newQualifiedArea==0);
for name={'flow','sst'}
    record=summary.replay.(name{1});
    assert(record.relaxationReplayAbsoluteError<2e-12 && record.updateReplayScaledError<2e-12);
end
assert(summary.replay.sst.savedStateScaledError<2e-12 && summary.maximumLogRMSDifference<2e-8);
files={'engine_process.json','sst_state_audit.json','restart_flow.csv','history.csv', ...
    'flow_update_before.csv','flow_update_after.csv','sst_update_before.csv','sst_update_after.csv'};
for i=1:numel(files)
    file=fullfile(folder,files{i});
    token=regexp(text,[regexptranslate('escape',jsonencode(file)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
    assert(numel(token)==1 && strcmp(sha(file),token{1}{1}),'实际更新审核缺少准确来源或文件已改变：%s',files{i});
end
summary.auditSHA256=sha(path);
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
