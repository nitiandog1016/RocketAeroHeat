function saveRFZReferenceBundle(file,samples)
% 模块一：校验真实样点并复用原有 CSV 导出，来源摘要绑定实际写入的文件。
source=rfzReferenceSourceRecord(samples);
folder=fileparts(file);if isempty(folder),folder=pwd;end
tmp=tempname(folder);[ok,msg]=mkdir(tmp);assert(ok,'%s',msg);
cleanup=onCleanup(@()clearStaging(tmp)); %#ok<NASGU>
csv=fullfile(tmp,'samples.csv');saveRFZReferenceCSV(csv,samples);
source.exportedCSV=struct('filename','samples.csv','sha256',cm011FileSHA256(csv),'rows',height(samples));
json=fullfile(tmp,'source.json');[fid,msg]=fopen(json,'w','n','UTF-8');assert(fid>=0,'%s',msg);
try
    fprintf(fid,'%s\n',jsonencode(source,'PrettyPrint',true));
    [msg,err]=ferror(fid);assert(err==0,'%s',msg);
catch ex
    fclose(fid);rethrow(ex);
end
assert(fclose(fid)==0,'来源文件未完整关闭。');

% 模块二：表与说明在同一压缩包内发布，避免部分写入或覆盖目标后才发现错误。
archive=fullfile(tmp,'reference.zip');
zip(archive,{'samples.csv','source.json'},tmp);
[ok,msg]=movefile(archive,file,'f');assert(ok,'%s',msg);
end

function clearStaging(folder)
for name={'samples.csv','source.json','reference.zip'}
    path=fullfile(folder,name{1});if isfile(path),delete(path);end
end
if isfolder(folder),rmdir(folder);end
end
