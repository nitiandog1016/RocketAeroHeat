function [pressure,temperature,summary]=readSU2CoolPropStateAudit(folder,receipt,rest)
%READSU2COOLPROPSTATEAUDIT 只读同库独立EOS审核，不在浏览结果时启动计算。
% 模块一：绑定原生守恒场、真实查询终态和精确物性输出SHA。
% 模块二：每行必须对应原节点；比内能零点由EOS定义，不套定比热正值假设。
assert(strcmp(receipt.runtimeFlavor,'rah-coolprop-v1'),'结果未声明受控物性运行包。');
path=fullfile(folder,'thermodynamics');summary=jsondecode(fileread(fullfile(path,'audit.json')));
assert(strcmp(summary.thermodynamicModel,'CoolProp-7.2.0-HEOS-Air') && summary.sourceStable);
assert(strcmp(sha(fullfile(folder,'restart_flow.csv')),summary.sourceStateSHA256));
assert(strcmp(sha(fullfile(path,'state_process.json')),summary.processReceiptSHA256));
process=jsondecode(fileread(fullfile(path,'state_process.json')));
assert(strcmp(process.status,'terminal') && process.exitCode==0 && ~process.timedOut);
file=fullfile(path,'properties.csv');assert(strcmp(sha(file),summary.propertyOutputSHA256));
data=readmatrix(file,'NumHeaderLines',1);
assert(size(data,1)==size(rest,1) && size(data,2)==17 && isequal(data(:,1),rest(:,1)));
nd=summary.spatialDimensions;assert(ismember(nd,[2 3]) && size(rest,2)==2*nd+3);
rho=rest(:,nd+2);u=rest(:,nd+3:2*nd+2)./rho;e=rest(:,end)./rho-sum(u.^2,2)/2;
assert(max(abs(data(:,4)./rho-1))<1e-8 && max(abs(data(:,5)-e)./max(abs(e),100000))<1e-8);
pressure=data(:,2);temperature=data(:,3);
assert(all(isfinite(data),'all') && all(pressure>0 & temperature>0 & data(:,7)>0 & data(:,8)>0));
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);g=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
