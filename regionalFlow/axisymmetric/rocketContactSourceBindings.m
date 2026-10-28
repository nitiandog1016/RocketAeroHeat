function bindings=rocketContactSourceBindings
%ROCKETCONTACTSOURCEBINDINGS 锁定完整原方程/父根工厂和新受控HLLC工况及实际驱动。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));bindings=rocketWallResolutionSourceBindings;
paths={'regionalFlow/axisymmetric/createRocketAxisymmetricContactCase.m', ...
    'regionalFlow/axisymmetric/rocketContactSourceBindings.m','tests/solveRocketContactStudy.m'};
for k=1:numel(paths)
    fid=fopen(fullfile(root,paths{k}),'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    bindings(end+1)=struct('path',paths{k},'sha256',lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]))); %#ok<AGROW>
end
end
