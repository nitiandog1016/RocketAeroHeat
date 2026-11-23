function binding=bindThermistorSTL(mesh,reference,source,geometry)
% 模块一：从冻结 STL 的实际面记录复核工作台朝向、单位和面序。
assert(strcmp(source.unit,'m') && source.scale==geometry.scale && strcmp(cm011FileSHA256(source.file),source.sha256), ...
    'RAH:Bead:Source','STL 文件或单位已改变，请重新加载模型。');
r=geometry.rotationMatrix;c=reshape(geometry.rotationCenter,1,3);
assert(all(isfinite(c)) && norm(r*r'-eye(3),'fro')<1e-10 && abs(det(r)-1)<1e-10, ...
    'RAH:Bead:Pose','朝向必须是已确认的刚体旋转。');
xyz=facetCorners(source.file,size(reference.F,1))*source.scale;
xyz=(xyz-c)*r'+c;
if size(xyz,1)~=numel(mesh.F)
    % 个别导入器会略掉退化三角形。按文件本身的记录定位缺失、重复或损坏的面。
    mapThermistorSTL(struct('V',xyz,'F',reshape(1:size(xyz,1),3,[])'),reference);
    error('RAH:Bead:SourceFaces','文件有 %d 面，工作台有 %d 面。',size(xyz,1)/3,size(mesh.F,1));
end
actual=mesh.V(reshape(mesh.F.',[],1),:);
% STL 读取器应保留文件面序；三角形自身的循环换位保留同一法向。
d=inf(size(mesh.F,1),1);
for shift=0:2
    idx=reshape(circshift(reshape(1:size(xyz,1),3,[]),shift,1),[],1);
    d=min(d,max(reshape(vecnorm(actual-xyz(idx,:),2,2),3,[]),[],1)');
end
bad=find(~(d<=16*eps(max([abs(xyz(:));abs(c(:));1e-3]))),1);
assert(isempty(bad),'RAH:Bead:FileFace','工作台面 %d 与文件中同号三角形的坐标或绕序不同。',bad);
binding=mapThermistorSTL(mesh,reference);
assert(strcmp(cm011FileSHA256(source.file),source.sha256),'RAH:Bead:Source','几何核对期间 STL 被修改。');
binding.source=source;binding.rotationMatrix=r;binding.rotationCenterM=c;
end

function xyz=facetCorners(path,expected)
% 模块二：按文件顺序读出角点。先用记录长度辨别二进制文件，ASCII 按 facet 收束。
file=dir(path);fid=fopen(path,'rb','ieee-le');assert(fid>=0);guard=onCleanup(@()fclose(fid)); %#ok<NASGU>
fseek(fid,80,'bof');n=fread(fid,1,'uint32');
if ~isempty(n) && file.bytes==84+50*n
    assert(n>0 && n<=2*expected,'RAH:Bead:FileCount','STL 有 %d 面，此试件需要 %d 面。',n,expected);
    fseek(fid,96,'bof');raw=fread(fid,[9 n],'9*single=>double',14);
    xyz=reshape(raw,3,[])';return;
end
assert(file.bytes<32*1024^2,'RAH:Bead:FileSize','B 珠文本 STL 大小超出支持范围。');
[xyz,facets]=readSTLTextCorners(path);
assert(facets>0 && facets<=2*expected && size(xyz,1)==3*facets,'RAH:Bead:FileCount','STL 有 %d 个 facet 和 %d 个角点，无法对应登记试件。',facets,size(xyz,1));
assert(all(isfinite(xyz),'all'),'RAH:Bead:Vertex','STL 含非有限角点。');
end
