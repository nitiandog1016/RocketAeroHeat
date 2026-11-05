function review=readSU2BenchmarkReview(file)
%READSU2BENCHMARKREVIEW 核对公开算例的真实终态、历史步和审核图来源。
% 模块一：原进程与历史为状态依据，失败记录和源初场分别命名。
% 模块二：检查新保存场是否确实存在，禁止把最后历史步当成结果帧。
% 模块三：本函数只返回公开验证记录，当前箭体输入与求解状态由工作台保持。
root=fileparts(mfilename('fullpath'));file=char(java.io.File(file).getCanonicalPath());
assert(startsWith(lower(file),lower([fullfile(root,'outputs') filesep])),'请选择本项目公开算例审核记录。');
raw=fileread(file);review=jsondecode(raw);
assert(strcmp(review.schema,'RocketAeroHeat.axisymmetric-study-review.v1') && review.auditCompleted);
folder=char(java.io.File(review.directory).getCanonicalPath());
assert(startsWith(lower(folder),lower([fullfile(root,'outputs') filesep])));
processFile=fullfile(folder,'engine_process.json');manifestFile=fullfile(folder,'study_manifest.json');
bind(raw,processFile);bind(raw,manifestFile);
receipt=jsondecode(fileread(processFile));manifest=jsondecode(fileread(manifestFile));
assert(strcmp(receipt.status,'terminal') && review.exitCode==receipt.exitCode);
assert(all(structfun(@(x)logical(x),receipt.inputsStable)) && review.inputsStable);
assert(strcmp(manifest.kind,'NASA-ACC-SBLI-M7-axisymmetric-convergence-v1'));
assert(strcmp(review.model,manifest.benchmark.model) && isequaln(review.conditions,manifest.benchmark.conditions));
assert(review.sourceLocalIteration==manifest.sourceLocalIteration);
for name={'history.csv','engine.log'}
    assert(strcmp(sha(fullfile(folder,name{1})),receipt.outputSHA256.(matlab.lang.makeValidName(name{1}))));
end
history=readmatrix(fullfile(folder,'history.csv'),'NumHeaderLines',1);
assert(all(isfinite(history(:))) && isequal(history(:,3),(0:size(history,1)-1)'));
assert(review.lastRecordedIteration==history(end,3));
saved=dir(fullfile(folder,'restart_flow_*.csv'));
assert(review.hasSavedNewField==~isempty(saved),'新检查点存在性与审核记录不同。');
if isempty(saved)
    assert(isempty(review.savedLocalIteration) && isempty(review.currentFieldAssessment));
else
    name=sprintf('restart_flow_%06d.csv',review.savedLocalIteration);
    assert(isfile(fullfile(folder,name)) && ...
        strcmp(sha(fullfile(folder,name)),receipt.outputSHA256.(matlab.lang.makeValidName(name))));
end
assert(~review.physicsQualified && review.newQualifiedArea==0);
if receipt.exitCode~=0,assert(~review.validationPassed && ~review.runSucceeded);end
imageFile=char(java.io.File(review.historyPlot.path).getCanonicalPath());
assert(startsWith(lower(imageFile),lower([fileparts(file) filesep])) && ...
    strcmp(sha(imageFile),review.historyPlot.sha256),'审核图来源或字节已改变。');
review.reviewFile=file;review.reviewSHA256=sha(file);review.imageFile=imageFile;
end

function bind(raw,file)
token=regexp(raw,[regexptranslate('escape',jsonencode(file)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
assert(numel(token)==1 && strcmp(sha(file),token{1}{1}),'原运行身份缺失或已改变。');
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
