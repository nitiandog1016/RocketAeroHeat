function report=confirmRocketAxisymmetricSensitivity(kind)
%CONFIRMROCKETAXISYMMETRICSENSITIVITY 新网格已求根后的真实三步BE独立确认。
% 模块一：按新根实际SHA构造新确认链，插值和Newton不计作BE步。
% 模块二：沿用独立逐端点持久化，所有原标准不变，不覆盖旧工况。
root=fileparts(fileparts(mfilename('fullpath')));folder=fullfile(root,'outputs','v76_validation',kind);
path=fullfile(folder,'spatial_root_evidence.mat');
if ~exist(path,'file'),folder=fullfile(root,'validation','v76',kind);path=fullfile(folder,'spatial_root_evidence.mat');end
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
report=studyRocketAxisymmetricSensitivityRoot(fullfile(folder,'confirmed'),3,struct('kind',kind,'rootEvidenceSHA256',hash));
assert(report.acceptedSteps==3&&report.steadyConverged,'RocketSensitivity:Confirmation','新根尚未满足原三步BE定常确认。');
end
