function saveRFZReferenceCSV(file,samples)
% 模块用途：按作者行序导出，保留可往返双精度的有效数字。完整写入后替换目标。
names={'sampleID','authorX_m','authorY_m','authorZ_m','authorPressure_Pa','authorHeatFlux_Wm2','authorWallTemperature_K'};
assert(isequal(samples.Properties.VariableNames,names) && isequal(samples.sampleID,(1:height(samples))'), ...
    'RAH:RFZ:Reference','参考表的字段或样点编号不完整。');
folder=fileparts(file);if isempty(folder),folder=pwd;end
partial=[tempname(folder) '.csv'];cleanup=onCleanup(@()removePartial(partial)); %#ok<NASGU>
[fid,msg]=fopen(partial,'w','n','UTF-8');assert(fid>=0,'%s',msg);
try
    fprintf(fid,'%s\n',strjoin(names,','));
    fprintf(fid,'%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n',samples{:,:}');
    [msg,code]=ferror(fid);assert(code==0,'%s',msg);
catch ex
    fclose(fid);rethrow(ex);
end
assert(fclose(fid)==0,'参考样点文件未完整关闭。');
[ok,msg]=movefile(partial,file,'f');assert(ok,'%s',msg);
end

function removePartial(path)
if isfile(path),delete(path);end
end
