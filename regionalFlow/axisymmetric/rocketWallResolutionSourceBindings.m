function bindings=rocketWallResolutionSourceBindings
%ROCKETWALLRESOLUTIONSOURCEBINDINGS 在原完整方程依赖之外锁定新网格/初猜模块。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));bindings=rocketSensitivitySourceBindings;
names={'refineRocketAxisymmetricRadial','createRocketAxisymmetricWallResolutionCase','rocketWallResolutionSourceBindings'};
for k=1:numel(names)
    path=['regionalFlow/axisymmetric/' names{k} '.m'];fid=fopen(fullfile(root,path),'rb');assert(fid>=0);
    guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    bindings(end+1)=struct('path',path,'sha256',lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]))); %#ok<AGROW>
end
path='tests/solveRocketWallResolution.m';fid=fopen(fullfile(root,path),'rb');assert(fid>=0);
guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
bindings(end+1)=struct('path',path,'sha256',lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[])));
end
