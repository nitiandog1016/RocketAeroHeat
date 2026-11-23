function exportShortSupportedSphere(g,path)
% 模块一：把 SI 几何、原球面编号、两个接头和七站的定位一起导出。
arguments
    g (1,1) struct
    path (1,:) char
end
[folder,name,ext]=fileparts(path);
if isempty(folder),folder=pwd;end
assert(isfolder(folder) && strcmpi(ext,'.zip') && ~isempty(name) && ~isfile(path), ...
    'RAH:Short:NewArchive','请选择未使用的 ZIP 文件名。');
tmp=tempname(folder);mkdir(tmp);
names={'geometry.json','summary.json','parents.csv','stations.csv','supported_m.stl'};
cleanup=onCleanup(@()discard(tmp,names)); %#ok<NASGU>
writeJSON(fullfile(tmp,names{1}),rmfield(g,{'summary','parents','stations'}));
writeJSON(fullfile(tmp,names{2}),g.summary);
writetable(g.parents,fullfile(tmp,names{3}));writetable(g.stations,fullfile(tmp,names{4}));
stlwrite(triangulation(g.F,g.V),fullfile(tmp,names{5}));
zip(fullfile(tmp,'complete.zip'),names,tmp);
% 同卷发布，只移动已经完整写好的压缩包。
java.nio.file.Files.move(java.io.File(fullfile(tmp,'complete.zip')).toPath(), ...
    java.io.File(fullfile(folder,[name ext])).toPath(),javaArray('java.nio.file.CopyOption',0));
end

function writeJSON(path,s)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(s));
end

function discard(folder,names)
% 模块二：只清理本次导出的临时目录中已知的文件。
for item=[names {'complete.zip'}]
    p=fullfile(folder,item{1});if isfile(p),delete(p);end
end
rmdir(folder);
end
