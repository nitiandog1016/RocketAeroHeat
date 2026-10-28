function identity=rocketMixedCoarseLiftIdentity
%ROCKETMIXEDCOARSELIFTIDENTITY 绑定混合及原始变量提升公式，不跨公式复用缓存。
root=fileparts(mfilename('fullpath'));identity=struct;
for entry={'liftRocketMixedCoarse','liftRocketPrimitiveCoarse'}
    name=entry{1};fid=fopen(fullfile(root,[name '.m']),'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));identity.(name)=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
end
