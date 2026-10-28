function bindings=rocketSensitivitySourceBindings
%ROCKETSENSITIVITYSOURCEBINDINGS 锁定新网格初猜、完整空间方程和实际求根依赖。
% 显示器不参与求解身份；旧v75根和旧检查点保持原字节，不重新绑定。
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));folder=fileparts(mfilename('fullpath'));
files=dir(fullfile(folder,'axisymmetric*.m'));names={files.name};
names=[names,{'buildAxisymmetricFVGeometry.m','buildRocketAxisymmetricDomain.m', ...
    'buildRocketAxisymmetricCollarDomain.m','createRocketAxisymmetricCase.m','createRocketAxisymmetricWallCase.m', ...
    'createRocketAxisymmetricCollarCase.m','createRocketAxisymmetricPressureCase.m','createRocketAxisymmetricFluxCase.m', ...
    'createRocketAxisymmetricLocalCase.m','createRocketAxisymmetricRootCase.m','createRocketAxisymmetricSensitivityCase.m', ...
    'mapRocketAxisymmetricInitialGuess.m','rocketSensitivitySourceBindings.m','solveAxisymmetricSteadyNewton.m', ...
    'advanceAxisymmetricSteadyLocal.m'}];
paths=cellfun(@(n)['regionalFlow/axisymmetric/' n],unique(names),'UniformOutput',false);
paths=[paths,{'buildRocketTestBody.m','standardAtmosphere.m','calcFreestream.m','createCaseConfig.m'}];
bindings=struct('path',{},'sha256',{});
for k=1:numel(paths)
    fid=fopen(fullfile(root,paths{k}),'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear guard;
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    bindings(k)=struct('path',paths{k},'sha256',lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[])));
end
end
