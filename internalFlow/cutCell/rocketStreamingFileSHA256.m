function digest=rocketStreamingFileSHA256(filename)
%ROCKETSTREAMINGFILESHA256 分块绑定大MAT来源，避免把完整文件一次读入内存。
% 模块一：SHA仅用于来源身份，不从文件内容执行代码或获取外部场数据。
fid=fopen(filename,'rb');assert(fid>=0,'SteadyGuess:File','来源文件不可读。');guard=onCleanup(@()fclose(fid));
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),bytes=fread(fid,8*1024^2,'*uint8');if ~isempty(bytes),md.update(typecast(bytes,'int8'));end,end
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
