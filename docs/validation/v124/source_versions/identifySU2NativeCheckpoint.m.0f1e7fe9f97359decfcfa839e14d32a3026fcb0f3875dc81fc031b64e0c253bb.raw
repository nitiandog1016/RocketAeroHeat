function identity=identifySU2NativeCheckpoint(folder,receipt,config)
%IDENTIFYSU2NATIVECHECKPOINT 用SU2原生具名文件绑定迭代，避免热负荷平台歧义。
% 模块一：仅接受单区稳态、原生六位以上迭代后缀及同时完整的壁面/重启动。
% 模块二：具名文件须与回执SHA及当前基础文件逐字相同；不由mtime或退出0猜时刻。
% 模块三：多个时刻逐字同态时选择最新的完整具名对，明确这是本次运行局部迭代。
identity=struct('enabled',false,'confirmed',false,'localIteration',NaN,'method','native-indexed-pair-SHA', ...
    'reason','旧配置未启用原生具名检查点。','files',struct,'identicalPairIterations',[]);
enabled=@(key)~isempty(regexp(config,['(?m)^\s*' key '\s*=\s*NO\s*$'],'once'));
if ~enabled('WRT_RESTART_OVERWRITE') || ~enabled('WRT_SURFACE_OVERWRITE'),return;end
identity.enabled=true;identity.reason='缺少完整同迭代具名文件对，不能推断显示迭代。';
assert(isempty(regexp(config,'(?m)^\s*(TIME_DOMAIN|MULTIZONE_MESH)\s*=\s*YES','once')),'原生时刻绑定仅支持单区稳态。');
hashes=receipt.outputSHA256;
assert(isfield(hashes,'restart_flow_csv') && isfield(hashes,'surface_flow_vtk'),'缺少基础文件散列。');
assert(strcmp(sha(fullfile(folder,'restart_flow.csv')),hashes.restart_flow_csv),'基础重启动改变。');
assert(strcmp(sha(fullfile(folder,'surface_flow.vtk')),hashes.surface_flow_vtk),'基础壁面改变。');
list=dir(fullfile(folder,'restart_flow_*.csv'));iterations=[];pairs=cell(0,1);
for j=1:numel(list)
    token=regexp(list(j).name,'^restart_flow_([0-9]{6,})\.csv$','tokens','once');if isempty(token),continue;end
    restart=list(j).name;wall=['surface_flow_' token{1} '.vtk'];
    if ~isfile(fullfile(folder,wall)),continue;end
    rk=matlab.lang.makeValidName(restart);wk=matlab.lang.makeValidName(wall);
    if ~isfield(hashes,rk) || ~isfield(hashes,wk),continue;end
    rh=sha(fullfile(folder,restart));wh=sha(fullfile(folder,wall));
    assert(strcmp(rh,hashes.(rk)) && strcmp(wh,hashes.(wk)),'原生具名文件SHA与真实回执不一致。');
    if strcmp(rh,hashes.restart_flow_csv) && strcmp(wh,hashes.surface_flow_vtk)
        iterations(end+1)=str2double(token{1}); %#ok<AGROW>
        pairs{end+1}=struct('restart',restart,'surface',wall,'restartSHA256',rh,'surfaceSHA256',wh); %#ok<AGROW>
    end
end
if isempty(iterations),return;end
[latest,k]=max(iterations);identity.confirmed=true;identity.localIteration=latest;identity.files=pairs{k};
identity.identicalPairIterations=iterations;identity.reason='最新完整原生具名文件对与当前壁面/体场逐字相同；本次运行局部迭代。';
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
