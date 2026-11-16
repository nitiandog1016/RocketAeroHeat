function digest=phaseAwarePayloadSHA(metadata,varargin)
%PHASEAWAREPAYLOADSHA 绑定原面场值/掩码与来源文字；固定double列优先字节序。
% 模块用途：生产器与工作台出口使用相同规则识别数据篡改或旧字段拼接。
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(uint8(unicode2native(metadata,'UTF-8')),'int8'));
for k=1:numel(varargin)
    value=double(varargin{k});dims=double(size(value));
    md.update(typecast(dims(:),'int8'));md.update(typecast(value(:),'int8'));
end
digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
