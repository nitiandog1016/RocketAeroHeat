function digest=phaseAwareMeshIdentity(mesh)
%PHASEAWAREMESHIDENTITY 绑定原坐标、面元连接及排序，防止旧外缘包串用。
% 模块用途：以维度及double原始字节生成SHA；本接口只用于本机同一MATLAB数据链。
V=double(mesh.V); F=double(mesh.F);
validateattributes(V,{'double'},{'real','finite','ncols',3,'nonempty'});
validateattributes(F,{'double'},{'real','finite','integer','positive','ncols',3,'nonempty','<=',size(V,1)});
md=java.security.MessageDigest.getInstance('SHA-256');
for a={double([size(V) size(F)]),V,F}
    md.update(typecast(a{1}(:),'int8'));
end
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
