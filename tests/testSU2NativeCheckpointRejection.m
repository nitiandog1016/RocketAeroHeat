function testSU2NativeCheckpointRejection(out)
%TESTSU2NATIVECHECKPOINTREJECTION 检查具名证据缺失及散列错误时的拒绝行为。
% 模块一：复用真实完整输出，只改内存中的回执，不触碰已有文件或启动求解。
% 模块二：基础文件齐全但具名SHA缺失应保持未知；具名SHA错误应明确拒绝。
% 模块三：旧覆盖写配置不启用新绑定，保持历史数据原审核路径。
assert(~isfolder(out));mkdir(out);
folder='outputs/v111_validation/control01';
receipt=jsondecode(fileread(fullfile(folder,'engine_process.json')));
config=fileread(fullfile(folder,'case.cfg'));
missing=receipt;missing.outputSHA256=rmfield(missing.outputSHA256,'restart_flow_000599_csv');
a=identifySU2NativeCheckpoint(folder,missing,config);
assert(a.enabled && ~a.confirmed && isnan(a.localIteration));
changed=receipt;changed.outputSHA256.restart_flow_000599_csv=repmat('0',1,64);
rejected=false;reason='';
try,identifySU2NativeCheckpoint(folder,changed,config);
catch ex,rejected=contains(ex.message,'SHA');reason=ex.message;end
assert(rejected,'具名SHA错误没有按预期拒绝。');
old=strrep(config,'WRT_RESTART_OVERWRITE= NO','WRT_RESTART_OVERWRITE= YES');
b=identifySU2NativeCheckpoint(folder,receipt,old);assert(~b.enabled && ~b.confirmed);
report=struct('missingIndexedReceiptRemainsUnknown',true,'changedIndexedHashRejected',rejected, ...
    'rejectionReason',reason,'legacyConfigurationKeepsOldPath',true, ...
    'syntheticReceiptOnly',true,'originalFilesUnmodified',true,'solverInvoked',false);
fid=fopen(fullfile(out,'report.json'),'w','n','UTF-8');assert(fid>=0);
guard=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
fprintf('Native indexed evidence rejection checks passed without solver runs.\n');
end
