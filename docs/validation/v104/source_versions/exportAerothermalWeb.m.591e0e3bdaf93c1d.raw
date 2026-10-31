function path=exportAerothermalWeb(mesh,folder)
%EXPORTAEROTHERMALWEB 输出无需服务器和外部依赖的完整结果查看器目录。
% 模块一：复制本地网页；模块二：写入 JSON 数据脚本，NaN 转为 null。
if ~exist(folder,'dir'),mkdir(folder);end
root=fileparts(mfilename('fullpath'));
for name={'index.html','viewer.css','section.js','export.js','viewer.js'}
    copyfile(fullfile(root,'web',name{1}),fullfile(folder,name{1}));
end
payload=buildAerothermalViewerData(mesh);
fid=fopen(fullfile(folder,'case-data.js'),'w','n','UTF-8');
assert(fid>=0,'RocketAeroHeat:Viewer:Export','不能写入导出目录。');
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'window.ROCKET_CASE=%s;\n',jsonencode(payload,'ConvertInfAndNaN',true));
path=fullfile(folder,'index.html');
end
