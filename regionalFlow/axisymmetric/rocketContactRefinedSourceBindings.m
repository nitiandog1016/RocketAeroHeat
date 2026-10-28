function bindings=rocketContactRefinedSourceBindings
%ROCKETCONTACTREFINEDSOURCEBINDINGS 锁定原HLLC/新粗根插值及实际精网格驱动。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));bindings=rocketContactSourceBindings;
paths={'regionalFlow/axisymmetric/createRocketAxisymmetricContactRefinedCase.m', ...
    'regionalFlow/axisymmetric/rocketContactRefinedSourceBindings.m','tests/solveRocketContactRefinedStudy.m'};
for k=1:numel(paths)
    fid=fopen(fullfile(root,paths{k}),'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    bindings(end+1)=struct('path',paths{k},'sha256',lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]))); %#ok<AGROW>
end
end
