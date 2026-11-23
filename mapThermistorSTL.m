function map=mapThermistorSTL(mesh,body)
% 模块一：在用户已确认的坐标系中消除平移，逐角点核对带线 B 珠的米制几何。
% 不调整滚转或尺度。用户模型的顶点、面序及跨圆角的区域边界均保留。
tol=2e-10;
v=mesh.V;f=mesh.F;r=body.V;rf=body.F;
offset=(max(v)+min(v)-max(r)-min(r))/2;
aligned=v-offset;
[found,indices]=ismembertol(aligned,r,tol,'ByRows',true,'DataScale',1,'OutputAllIndices',true);
uniqueVertex=found & cellfun(@numel,indices)==1;
bad=find(any(~uniqueVertex(f),2),1);
assert(isempty(bad),'RAH:Bead:Vertex','导入面 %d 的角点无唯一对应。请核对单位、朝向和几何。',bad);
id=cellfun(@(x)x(1),indices(f));id=reshape(id,size(f));
distance=vecnorm(aligned(f(:),:)-r(id(:),:),2,2);
bad=find(any(reshape(distance,size(f))>tol,2),1);
assert(isempty(bad),'RAH:Bead:Vertex','导入面 %d 超出角点距离容差。',bad);
% 模块二：以三个顶点的完整集合定位三角形，再单独检查绕序和一一覆盖。
[found,order]=ismember(sort(id,2),sort(rf,2),'rows');
bad=find(~found,1);
assert(isempty(bad),'RAH:Bead:Triangle','导入面 %d 的三个角点未构成登记三角形。',bad);
[sorted,which]=sort(order);bad=find(diff(sorted)==0,1);
if ~isempty(bad)
    error('RAH:Bead:Duplicate','导入面 %d 和 %d 重复对应标准面 %d。',which(bad),which(bad+1),sorted(bad));
end
missing=setdiff(1:size(rf,1),order);
if ~isempty(missing),error('RAH:Bead:Missing','缺少标准面 %d 的导入三角形。',missing(1));end
n=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
standard=rf(order,:);rn=cross(r(standard(:,2),:)-r(standard(:,1),:),r(standard(:,3),:)-r(standard(:,1),:),2);
area=vecnorm(n,2,2)/2;refArea=vecnorm(rn,2,2)/2;
cosine=sum((n./(2*area)).*(rn./(2*refArea)),2);
bad=find(~(cosine>.999999 & abs(area./refArea-1)<=1e-5),1);
assert(isempty(bad),'RAH:Bead:NormalArea','导入面 %d 法向或实际面积与标准面不符。',bad);
map=struct('referenceFaceID',order,'translationM',offset,'areaM2',area,'referenceAreaM2',refArea, ...
    'maxVertexErrorM',max(distance),'minNormalCosine',min(cosine),'maxRelativeAreaError',max(abs(area./refArea-1)));
end
