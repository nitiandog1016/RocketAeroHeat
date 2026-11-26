function [region,boundary] = bindRFZGeometry(mesh,filename,scale,fileHash)
% 模块一：识别 RFZ-REC 的冻结米制 STL，普通文件继续使用原工作流。
region=[];boundary=[];
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

% 模块四：按原 CAD 面号加入作者封盖的归属。人工多边形只存于独立清单。
source=fullfile(folder,'deg190_boundary.json');condition=fullfile(folder,'deg190_reference.json');
assert(strcmp(cm011FileSHA256(source),'b6647eb6dd7cb68e8b0db272240bdae0a9f4f686c660f63cf1eb62e1758de073'), ...
    'RAH:RFZ:Boundary','RFZ 边界清单与发布版本不符。');
b=jsondecode(fileread(source));
assert(strcmp(b.caseSHA256,cm011FileSHA256(condition)) && strcmp(b.stlSHA256,fileHash) && ...
    strcmp(b.cadTableSHA256,info.regionsSHA256) && strcmp(b.stepSHA256,info.stepSHA256));
boundary=struct('condition',jsondecode(fileread(condition)));
assert(strcmp(b.dataset,boundary.condition.dataset) && b.faces==size(mesh.F,1) && b.cadFaces==height(region));
caps=b.caps;rim=[caps.rimCADFace]';inner=cellfun(@(v)v(:),{caps.innerCADFaces},'UniformOutput',false);
allAssigned=[rim;vertcat(inner{:})];
assert(numel(unique(allAssigned))==numel(allAssigned) && all(ismember(allAssigned,region.cadFace)));
boundary.regions=region;
boundary.regions.referenceRole=repmat("其余 CAD 壁",height(region),1);
boundary.regions.referenceRole(rim)="实体唇环";
boundary.regions.referenceRole(vertcat(inner{:}))="封盖内侧 CAD 壁";
boundary.regions.rimCADFace=nan(height(region),1);
for j=1:numel(rim)
    boundary.regions.rimCADFace([rim(j);inner{j}])=rim(j);
end
keys={'other','rim','cavity'};labels=["其余 CAD 壁","实体唇环","封盖内侧 CAD 壁"];
for j=1:numel(keys)
    selected=boundary.regions.referenceRole==labels(j);g=b.groups.(keys{j});
    assert(nnz(selected)==g.cadFaces && sum(region.lastFace(selected)-region.firstFace(selected)+1)==g.faces && ...
        abs(sum(areas(selected))-g.areaM2)<1e-8);
end
boundary.caps=table(rim,[caps.polygonAreaM2]',[caps.referenceInteriorPoints]',[caps.farPoints]', ...
    'VariableNames',{'rimCADFace','polygonAreaM2','referenceInteriorPoints','farPoints'});
boundary.groups=b.groups;
end
