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
if isfield(mesh,'shortSupportStationTable')
    writeShortStations(path,folder,mesh.shortSupportStationTable,mesh.evidenceIdentity.sourceFrame.sha256);
end
if isfield(mesh,'shortSupportAzimuth')
    writeAzimuth(path,folder,root,mesh.shortSupportAzimuth);
end
end

function writeAzimuth(path,folder,root,az)
% 模块用途：离线页与桌面共用同一个可交互图，完整 CSV 带实际和参考面号以及源帧。
base=['Short65-azimuth-' az.summary.frameSHA256];page=[base '.html'];sheet=[base '.csv'];
t=az.points;t.frameSHA256=repmat(string(az.summary.frameSHA256),height(t),1);
t.stlSHA256=repmat(string(az.stlSHA256),height(t),1);writetable(t,fullfile(folder,sheet),'Encoding','UTF-8');
template=fileread(fullfile(root,'web','short65-azimuth.html'));marker='/*SHORT65_OFFLINE_DATA*/';
assert(count(template,marker)==1);data=jsonencode(shortAzimuthPayload(az),'ConvertInfAndNaN',true);
writeHtml(fullfile(folder,page),strrep(template,marker,['window.SHORT65_DATA=' data ';']));
html=fileread(path);assert(count(html,'</aside>')==1);
block=['<details id="short65-azimuth"><summary>七站周向换热曲线</summary>' ...
    '<iframe loading="lazy" title="七站周向换热" style="width:100%;height:580px;border:0" src="' page '"></iframe>' ...
    '<p><a href="' page '">展开曲线</a> · <a download href="' sheet '">下载全周读数与取样状态 CSV</a></p></details>'];
writeHtml(path,strrep(html,'</aside>',[block '</aside>']));
end
function writeHtml(path,html)
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'无法写入周向换热页。');
close=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',html);
end

function writeShortStations(path,folder,stations,frame)
% 模块三：把七站对照放入本次离线三维页的侧栏，下载表保留原面号、算法、间距和源帧。
assert(height(stations)==7 && all(stations.gradientFrameSHA256==string(frame)));
name=['Short65-stations-' frame '.csv'];writetable(stations,fullfile(folder,name),'Encoding','UTF-8');
rows=strings(height(stations),1);
for j=1:height(stations)
    value='缺少法向取样';
    if isfinite(stations.gradientNu(j)),value=sprintf('%.4f',stations.gradientNu(j));end
    rows(j)=sprintf('<tr><td>%g</td><td>%g</td><td>%.4f</td><td>%s</td></tr>', ...
        stations.angleDeg(j),stations.reportedNu(j),stations.computedNu(j),value);
end
block=['<details id="short65-stations"><summary>七站换热对照</summary>' ...
    '<table><caption>局部 Nu · 以当地壁温取导热系数</caption><thead><tr>' ...
    '<th>极角 (°)</th><th>实验</th><th>壁反力</th><th>温度梯度</th></tr></thead><tbody>' ...
    char(join(rows,'')) '</tbody></table><p>梯度取自壁点与前两层流体温度。' ...
    '表 III 的实验 Nu 经整体量热标定。</p><a download href="' name '">下载七站热流、温度取样信息与算法 CSV</a></details>'];
html=fileread(path);assert(count(html,'</aside>')==1);
fid=fopen(path,'w','n','UTF-8');assert(fid>=0,'无法写入七站离线页。');
close=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',strrep(html,'</aside>',[block '</aside>']));
end
