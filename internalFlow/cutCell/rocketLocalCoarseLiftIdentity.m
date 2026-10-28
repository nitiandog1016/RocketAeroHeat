function identity=rocketLocalCoarseLiftIdentity
%ROCKETLOCALCOARSELIFTIDENTITY 绑定局部保守方向及继承混合提升的完整公式。
% 模块一：新增局部提升单独绑定，不修改或重新解释旧473列缓存的提升来源。
identity=rocketMixedCoarseLiftIdentity;root=fileparts(mfilename('fullpath'));name='liftRocketLocalCoarse';
fid=fopen(fullfile(root,[name '.m']),'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));
identity.(name)=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
