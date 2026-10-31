function summary=readSU2SSTStateAudit(folder,receipt)
%READSU2SSTSTATEAUDIT 绑定独立SST热状态审核，防止缺证据的研究快照绕过后验。
% 模块一：确认审核属于实际原生回执及其具名输出，不读取当前默认模型替代。
% 模块二：逐点热梯度/EOS/y+维持原门槛，单步阶段诊断保持独立身份。
% 模块三：接口通过只允许研究查看，不授权收敛、完整区域精度或壁温解。
path=fullfile(folder,'sst_state_audit.json');assert(isfile(path),'缺少SST同阶段热状态审核，不能接收研究结果。');
text=fileread(path);summary=jsondecode(text);
assert(isfield(receipt,'physicsModel') && isequaln(summary.model,receipt.physicsModel),'SST审核模型与原任务不符。');
assert(summary.interfacePassed && ~summary.physicsQualified && summary.newQualifiedArea==0);
assert(summary.heatGradientRelativeError<2e-6 && summary.EOSRelativeError<2e-6 && summary.yPlusRelativeError<2e-6);
file=fullfile(folder,'engine_process.json');original=fullfile(receipt.request.workDirectory,'engine_process.json');
token=regexp(text,[regexptranslate('escape',jsonencode(original)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
assert(numel(token)==1 && strcmp(sha(file),token{1}{1}),'SST审核未绑定本次原生终态。');
for name={'restart_flow.csv','surface_flow.vtk','history.csv'}
    assert(strcmp(sha(fullfile(folder,name{1})),receipt.outputSHA256.(matlab.lang.makeValidName(name{1}))), ...
        'SST审核后的原场已改变。');
end
if isfield(receipt.physicsModel,'diagnosticOnly') && receipt.physicsModel.diagnosticOnly
    assert(summary.diagnosticOnly && summary.checkpointLocalIteration==0 && ~summary.originalResidualThresholdPassed);
    assert(strcmp(summary.thermalStatePhase,'updated-flow-with-input-k; single diagnostic update only'));
end
summary.auditSHA256=sha(path);
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
