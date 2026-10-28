function report=confirmRocketWallResolution(kind)
%CONFIRMROCKETWALLRESOLUTION 按新根真实SHA构造独立原三步BE确认链。
% 原Newton准备和插值不计BE，原门槛、不可变端点及完整五量账保持。
root=fileparts(fileparts(mfilename('fullpath')));folder=fullfile(root,'outputs','v77_validation',kind);
if ~exist(fullfile(folder,'spatial_root_evidence.mat'),'file'),folder=fullfile(root,'validation','v77',kind);end
path=fullfile(folder,'spatial_root_evidence.mat');d=load(path,'c');
fid=fopen(path,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
o=d.c.options;o.rootEvidenceSHA256=hash;
report=studyRocketWallResolutionRoot(fullfile(folder,'confirmed'),3,o);
assert(report.acceptedSteps==3&&report.steadyConverged,'RocketWallResolution:Confirmation','尚未满足原三步BE确认。');
end
