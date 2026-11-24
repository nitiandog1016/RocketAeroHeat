function region = bindRFZGeometry(mesh,filename,scale,fileHash)
% 模块一：识别 RFZ-REC 的冻结米制 STL，普通文件继续使用原工作流。
region=[];
if ~strcmp(fileHash,'f68c1b1416f58be69a94c062b8f2c8fabaccf9966bddeb6b6e8991b5d73ad083') || scale~=1
    return
end
folder=fullfile(fileparts(mfilename('fullpath')),'cases','rfz_rec');
info=jsondecode(fileread(fullfile(folder,'geometry.json')));
tableFile=fullfile(folder,'cad_faces.csv');
assert(strcmp(info.stlSHA256,fileHash) && strcmp(cm011FileSHA256(tableFile),info.regionsSHA256), ...
    'RAH:RFZ:Source','原面清单与已登记的 RFZ-REC 文件不符。');
region=readtable(tableFile,'TextType','string');
assert(height(region)==info.cadFaces && isequal(region.cadFace,(1:height(region))') && ...
    isequal(region.firstFace,[1;region.lastFace(1:end-1)+1]) && all(region.lastFace>=region.firstFace) && ...
    region.lastFace(end)==info.faces && size(mesh.F,1)==info.faces && size(mesh.V,1)==info.vertices);

% 模块二：逐面回读二进制 STL 的三个有序顶点，使分区编号和导入器的实际面号对应。
[fid,msg]=fopen(filename,'rb','ieee-le');assert(fid>=0,'%s',msg);
closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
assert(fseek(fid,80,'bof')==0 && fread(fid,1,'uint32')==info.faces);
data=fread(fid,[12 info.faces],'12*single=>double',2);
assert(size(data,2)==info.faces,'RAH:RFZ:Read','RFZ 面表读取不完整。');
for corner=1:3
    rows=3*corner+(1:3);
    assert(isequal(data(rows,:)',mesh.V(mesh.F(:,corner),:)), ...
        'RAH:RFZ:FaceOrder','导入后三角面的坐标或编号发生变化。');
end

% 模块三：从 MATLAB 实际面元积分，逐区核对 204 个 CAD 面的三角化面积。
regionID=repelem(region.cadFace,region.lastFace-region.firstFace+1);
areas=accumarray(regionID,mesh.area(:));
assert(max(abs(areas-region.stlAreaM2))<1e-8 && abs(sum(areas)-info.stlAreaM2)<1e-8, ...
    'RAH:RFZ:Area','RFZ 分区面积核对失败。');
end
