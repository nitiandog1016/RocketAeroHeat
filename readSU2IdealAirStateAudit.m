function [pressure,temperature,summary]=readSU2IdealAirStateAudit(folder,receipt,rest)
%READSU2IDEALAIRSTATEAUDIT 绑定真实回执、独立物性和同阶段热梯度的只读入口。
% 模块一。原文件SHA及明确模型合同防止旧参考能量混入新模型。
% 模块二。读取原生流动阶段p/T；后更新k温度、给定Tw和资格单列。
assert(strcmp(receipt.runtimeFlavor,'rah-ideal-air-v1'));
index=jsondecode(fileread(fullfile(folder,'ideal_air_review.json')));
assert(strcmp(sha(fullfile(folder,'engine_process.json')),index.receiptSHA256));
base=char(java.io.File(folder).getCanonicalPath());
names={'auditPath','wallPath','thermoPath'};hashes={'auditSHA256','wallSHA256','thermoSHA256'};
for j=1:numel(names)
    p=char(java.io.File(fullfile(folder,index.(names{j}))).getCanonicalPath());
    assert(startsWith(lower(p),lower([base filesep])) && strcmp(sha(p),index.(hashes{j})));
end
text=fileread(fullfile(folder,index.auditPath));summary=jsondecode(text);
assert(summary.interfacePassed && ~summary.physicsQualified && summary.newQualifiedArea==0);
assert(isequaln(summary.thermodynamicModel,receipt.thermodynamicModel));
assert(summary.errors.phaseEnergy<2e-10 && summary.errors.EOS<2e-11 && summary.errors.heatGradient<2e-10);
assert(strcmp(summary.wallCSV_SHA256,index.wallSHA256));
for name={'restart_flow.csv','surface_flow.vtk','history.csv'}
    assert(strcmp(sha(fullfile(folder,name{1})),receipt.outputSHA256.(matlab.lang.makeValidName(name{1}))));
end
for name={'engine_process.json','restart_flow.csv'}
    original=fullfile(receipt.request.workDirectory,name{1});
    token=regexp(text,[regexptranslate('escape',jsonencode(original)) '\s*:\s*"([0-9a-f]{64})"'],'tokens');
    assert(numel(token)==1 && strcmp(sha(fullfile(folder,name{1})),token{1}{1}));
end
fid=fopen(fullfile(folder,'restart_flow.csv'),'r');assert(fid>=0);guard=onCleanup(@()fclose(fid));
header=fgetl(fid);clear guard;columns=regexp(header,'"([^"]+)"','tokens');columns=cellfun(@(x)x{1},columns,'UniformOutput',false);
if nargin<3,rest=readmatrix(fullfile(folder,'restart_flow.csv'),'NumHeaderLines',1);end
assert(size(rest,2)==numel(columns) && size(rest,1)==summary.nodes);
pressure=rest(:,strcmp(columns,'Pressure'));temperature=rest(:,strcmp(columns,'Temperature'));
rho=rest(:,strcmp(columns,'Density'));
assert(all(isfinite(rest),'all') && all(rho>0 & pressure>0 & temperature>=60 & temperature<=2000));
oracle=readtable(fullfile(folder,index.thermoPath),'VariableNamingRule','preserve');
assert(all(oracle.accepted==1) && isequal(oracle.row,rest(:,1)));
assert(max(abs(oracle.p-pressure)./max(1,abs(pressure)))<2e-11);
assert(max(abs(oracle.T-temperature))<2e-10);
assert(max(abs(pressure./(rho*receipt.thermodynamicModel.gasConstant.*temperature)-1))<2e-11);
summary.auditSHA256=index.auditSHA256;
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
