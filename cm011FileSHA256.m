function digest=cm011FileSHA256(path)
%CM011FILESHA256 分块计算 STL 或网格回执的文件指纹，绑定读取时的版本。
[fid,msg]=fopen(path,'rb');
assert(fid>=0,'RAH:CM011:File','无法读取 %s：%s',path,msg);
closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid)
    bytes=fread(fid,1048576,'*uint8');md.update(typecast(bytes,'int8'));
end
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
