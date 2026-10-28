function report=confirmRocketContactStudy(kind)
%CONFIRMROCKETCONTACTSTUDY 从新HLLC根真实SHA建立独立原三步BE确认链。
% 原HLLE初猜历史、新HLLC Newton准备和新BE计数分开，不改原接受标准。
root=fileparts(fileparts(mfilename('fullpath')));folder=fullfile(root,'outputs','v78_validation',kind);
if ~exist(fullfile(folder,'spatial_root_evidence.mat'),'file'),folder=fullfile(root,'validation','v78',kind);end
path=fullfile(folder,'spatial_root_evidence.mat');d=load(path,'c');
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
options=d.c.options;options.rootEvidenceSHA256=hash;
report=studyRocketContactRoot(fullfile(folder,'confirmed'),3,options);
assert(report.acceptedSteps==3&&report.steadyConverged,'RocketContact:Confirmation','原三步BE确认尚未通过。');
end
