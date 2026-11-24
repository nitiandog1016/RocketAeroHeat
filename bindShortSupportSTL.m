function [binding,g]=bindShortSupportSTL(mesh,source,geometry)
% 模块一：核对读取的文件、米制比例、刚体姿态与连续面记录，保留文件原号。
[ok,file]=fileattrib(geometry.stlFile);
assert(ok && strcmpi(file.Name,source.file) && strcmp(source.unit,'m') && source.scale>0 && ...
    source.scale==geometry.scale && strcmp(cm011FileSHA256(file.Name),source.sha256), ...
    'RAH:Short:Source','输入文件或尺度已改变，请重新导入。');
r=geometry.rotationMatrix;c=reshape(geometry.rotationCenter,1,3);
assert(geometry.isAligned && all(isfinite(c)) && norm(r*r'-eye(3),'fro')<1e-10 && abs(det(r)-1)<1e-10, ...
    'RAH:Short:Pose','请先确认实际模型朝向。');
folder=fullfile(fileparts(mfilename('fullpath')),'cases','short_sphere');
pair=readShortSupportTimePair(folder);archive=fullfile(folder,'supported65.zip');
assert(strcmp(cm011FileSHA256(archive),pair.bodySHA256) && ...
    strcmp(cm011FileSHA256(fullfile(folder,'supportedAngularTransfer65.npz')),pair.transferSHA256));
temp=tempname;mkdir(temp);cleanup=onCleanup(@()removeExtracted(temp)); %#ok<NASGU>
unzip(archive,temp);g=jsondecode(fileread(fullfile(temp,'geometry.json')));
g.stations=readtable(fullfile(temp,'stations.csv'));
n=size(g.F,1);v=mesh.V;f=mesh.F;
assert(size(f,1)==n && size(v,1)==size(g.V,1),'RAH:Short:Count','输入须有全部 %d 个原面和 %d 个共享顶点。',n,size(g.V,1));
xyz=fileCorners(source.file,n)*source.scale;xyz=(xyz-c)*r'+c;
actual=v(reshape(f.',[],1),:);difference=inf(n,1);
for shift=0:2
    idx=reshape(circshift(reshape(1:3*n,3,[]),shift,1),[],1);
    difference=min(difference,max(reshape(vecnorm(actual-xyz(idx,:),2,2),3,[]),[],1)');
end
assert(all(difference<=16*eps(max([abs(xyz(:));abs(c(:));.001]))),'RAH:Short:FileOrder','工作台三角面与同号 STL 记录不同。');
% 模块二：容许银球尺寸上两次 STL 单精度量化；逐点唯一且所有标准顶点均被使用。
offset=(max(v)+min(v)-max(g.V)-min(g.V))/2;aligned=v-offset;
tol=2*double(eps(single(max(abs(g.V(:))))));
[found,ids]=ismembertol(aligned,g.V,tol,'ByRows',true,'DataScale',1,'OutputAllIndices',true);
assert(all(found & cellfun(@numel,ids)==1),'RAH:Short:Vertex','角点不能唯一对应银球，请核对单位、朝向和形状。');
id=cell2mat(ids(:));d=vecnorm(aligned-g.V(id,:),2,2);
assert(numel(unique(id))==numel(id) && all(d<=tol),'RAH:Short:Vertex','顶点映射或坐标误差超界。');
[found,order]=ismember(sort(id(f),2),sort(g.F,2),'rows');
assert(all(found) && isequal(sort(order),(1:n)'),'RAH:Short:Topology','缺面、重复面或三角形连接不符。');
% 模块三：以循环绕序验方向，以实际边矢量验面积。文件编号只由记录顺序决定。
ref=g.F(order,:);mapped=id(f);
assert(all(all(mapped==ref,2) | all(mapped==ref(:,[2 3 1]),2) | all(mapped==ref(:,[3 1 2]),2)), ...
    'RAH:Short:Winding','导入三角面绕序与外法向不同。');
area=vecnorm(cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2),2,2)/2;
referenceArea=g.areaM2(order);relative=abs(area./referenceArea-1);
assert(all(relative<=4e-5),'RAH:Short:Area','面片面积超过该试件的单精度量化界限。');
assert(strcmp(cm011FileSHA256(file.Name),source.sha256),'RAH:Short:Source','核对期间 STL 已改变。');
binding=struct('source',source,'rotationMatrix',r,'rotationCenterM',c,'translationM',offset, ...
    'referenceFaceID',order,'areaM2',area,'referenceAreaM2',referenceArea,'regionID',g.regionID(order), ...
    'maxVertexErrorM',max(d),'maxRelativeAreaError',max(relative),'meshSHA256',pair.meshSHA256, ...
    'transferSHA256',pair.transferSHA256,'bodySHA256',pair.bodySHA256);
end

function xyz=fileCorners(file,n)
% 对二进制逐条跳过法向和属性字节；文本文件按原始 facet 读取，核查真实面数。
info=dir(file);fid=fopen(file,'rb','ieee-le');assert(fid>=0);close=onCleanup(@()fclose(fid)); %#ok<NASGU>
fseek(fid,80,'bof');count=fread(fid,1,'uint32');
if ~isempty(count) && info.bytes==84+50*count
    assert(count==n,'RAH:Short:FileCount','二进制 STL 面数不符。');
    fseek(fid,96,'bof');xyz=reshape(fread(fid,[9 n],'9*single=>double',14),3,[])';
else
    [xyz,count]=readSTLTextCorners(file);
    assert(count==n && size(xyz,1)==3*n,'RAH:Short:FileCount','文本 STL 面数不符。');
end
end

function removeExtracted(folder)
% 仅清理本次解包的已知输入。
for name={'geometry.json','summary.json','parents.csv','stations.csv','supported_m.stl'}
    file=fullfile(folder,name{1});if isfile(file),delete(file);end
end
rmdir(folder);
end
