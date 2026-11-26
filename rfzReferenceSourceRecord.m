function source=rfzReferenceSourceRecord(samples)
% 模块一：读取参考样点专用的来源记录。各来源的工况编号与 Mach 单独保存。
folder=fullfile(fileparts(mfilename('fullpath')),'cases','rfz_rec');
record=fullfile(folder,'deg190_source_record.json');
assert(strcmp(cm011FileSHA256(record),'0c5e36109a158bb04425223b1eaedeeb01b6d484a05f37ddc61c87f2d9dfeb7d'), ...
    'RAH:RFZ:ReferenceSource','参考资料记录与核对版本不符。');
source=jsondecode(fileread(record));
assert(strcmp(source.caseFileSHA256,cm011FileSHA256(fullfile(folder,'deg190_reference.json'))), ...
    'RAH:RFZ:ReferenceSource','参考资料记录与当前发布工况不符。');

% 模块二：只把来源附给完全匹配的原始样点。摘要按六列顺序及双精度小端字节计算。
fields={'sampleID','authorX_m','authorY_m','authorZ_m','authorPressure_Pa','authorHeatFlux_Wm2','authorWallTemperature_K'};
assert(isequal(samples.Properties.VariableNames,fields) && height(samples)==source.samples.count && ...
    isequal(samples.sampleID,(1:height(samples))'),'RAH:RFZ:ReferenceSource','样点行数、字段或顺序有误。');
md=java.security.MessageDigest.getInstance('SHA-256');[~,~,endian]=computer;
for i=2:numel(fields)
    v=samples.(fields{i});assert(isa(v,'double') && iscolumn(v) && isreal(v),'RAH:RFZ:ReferenceSource','字段类型不符。');
    if endian=='B',v=swapbytes(v);end
    md.update(typecast(v,'int8'));
end
actual=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
assert(strcmp(actual,source.samples.orderedValuesLEFloat64SHA256),'RAH:RFZ:ReferenceSource','样点原值与资料记录不符。');
end
