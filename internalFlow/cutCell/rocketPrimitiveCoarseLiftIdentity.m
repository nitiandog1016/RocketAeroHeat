function digest=rocketPrimitiveCoarseLiftIdentity
%ROCKETPRIMITIVECOARSELIFTIDENTITY 绑定粗原始变量提升源码，防止缓存跨提升公式误用。
path=fullfile(fileparts(mfilename('fullpath')),'liftRocketPrimitiveCoarse.m');fid=fopen(path,'rb');assert(fid>=0);bytes=fread(fid,Inf,'*uint8');fclose(fid);
md=java.security.MessageDigest.getInstance('SHA-256');md.update(typecast(bytes,'int8'));digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
