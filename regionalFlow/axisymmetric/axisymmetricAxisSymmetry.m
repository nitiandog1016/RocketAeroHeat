function symmetry=axisymmetricAxisSymmetry(mesh,boundary)
%AXISYMMETRICAXISSYMMETRY 为显式偏斜轴边建立镜像虚采样梯度模板。
% 模块一：仅响应axis.symmetryGradient='mirror-parity'，缺省不改变旧算子。
% 模块二：WLS轴样本移到(xOwner,-rOwner)，真实r=0面及积分度量始终保留。
% 模块三：偶量镜像同值、径向速度镜像反号由调用者提供；不移动真实流场单元。
% 模块四：重建模板只使用原一阶邻接，验证既有Jacobian保守支持足够。
% 这仅为无旋轴对称的数值边界模板，不是鼻尖/尾迹物理验证或生产面积资格。
names={'iMin','iMax','jMin','jMax'};
assert(isstruct(boundary)&&isscalar(boundary)&&all(isfield(boundary,names)), ...
    'AxisFV:AxisSymmetryBoundary','镜像轴模板需要完整的四条显式边界。');
selected=false(4,1);
for k=1:4
    bc=boundary.(names{k});
    if ~isstruct(bc)||~isscalar(bc)||~isfield(bc,'symmetryGradient'),continue;end
    policy=bc.symmetryGradient;
    assert(isfield(bc,'type')&&strcmp(bc.type,'axis')&& ...
        ((ischar(policy)&&isrow(policy))||(isstring(policy)&&isscalar(policy)&&~ismissing(policy)))&& ...
        strcmp(policy,'mirror-parity'), ...
        'AxisFV:AxisSymmetryBoundary','symmetryGradient仅允许真实axis边显式选择mirror-parity。');
    selected(k)=true;
end
symmetry=struct('active',false);
if ~any(selected),return;end
assert(isstruct(mesh)&&isscalar(mesh)&&isfield(mesh,'kind')&& ...
    strcmp(mesh.kind,'axisymmetric-body-fitted-fv-v1'), ...
    'AxisFV:AxisSymmetryGeometry','镜像模板仅接受完整贴体轴对称FV几何。');
f=mesh.faces;nc=mesh.cellCount;nf=numel(f.owner);bf=mesh.boundaryFaces(:);
faces=find(ismember(f.boundarySide,find(selected)));
assert(~isempty(faces)&&all(ismember(find(selected),unique(f.boundarySide(faces)))), ...
    'AxisFV:AxisSymmetryGeometry','显式轴边没有对应的真实外边界面。');
[present,boundaryIndex]=ismember(faces,bf);owner=f.owner(faces);
assert(all(present)&&all(f.neighbour(faces)==0)&&numel(unique(owner))==numel(owner)&& ...
    all(f.areaPerRadian(faces)==0)&&all(f.centroid(faces,2)==0)&& ...
    all(abs(f.normal(faces,1))<1e-12)&&all(abs(abs(f.normal(faces,2))-1)<1e-12)&& ...
    all(isfinite(mesh.volumeCentroid(owner,:)),'all')&&all(mesh.volumeCentroid(owner,2)>0), ...
    'AxisFV:AxisSymmetryGeometry','镜像模板须为真实r=0、零旋转面积、唯一owner且体质心r>0的轴面。');

%% 一、仅替换梯度虚样本位移，不改真实面坐标、真实距离或通量面积
delta=f.delta;ghost=[mesh.volumeCentroid(owner,1),-mesh.volumeCentroid(owner,2)];
delta(faces,:)=ghost-mesh.volumeCentroid(owner,:);
inside=f.neighbour>0;nb=numel(bf);boundaryColumn=zeros(nf,1);boundaryColumn(bf)=nc+(1:nb);
sample=f.neighbour;sample(bf)=boundaryColumn(bf);
cellID=[f.owner;f.neighbour(inside)];sampleID=[sample;f.owner(inside)];
d=[delta;-delta(inside,:)];weights=1./sum(d.^2,2);
assert(all(isfinite(weights)&weights>0),'AxisFV:AxisSymmetryGeometry','镜像梯度位移必须有限且非零。');
ax=accumarray(cellID,weights.*d(:,1).^2,[nc 1]);
ar=accumarray(cellID,weights.*d(:,2).^2,[nc 1]);
cross=accumarray(cellID,weights.*d(:,1).*d(:,2),[nc 1]);detA=ax.*ar-cross.^2;
assert(all(isfinite(detA)&detA>1e-13.*ax.*ar), ...
    'AxisFV:AxisSymmetryGeometry','镜像轴WLS法方程欠秩或病态，未放宽原几何秩门槛。');
cx=weights.*(ar(cellID).*d(:,1)-cross(cellID).*d(:,2))./detA(cellID);
cr=weights.*(ax(cellID).*d(:,2)-cross(cellID).*d(:,1))./detA(cellID);
Gx=sparse([cellID;cellID],[sampleID;cellID],[cx;-cx],nc,nc+nb);
Gr=sparse([cellID;cellID],[sampleID;cellID],[cr;-cr],nc,nc+nb);

%% 二、虚样本仍只由同一owner原始量提供，不引入全轴求解或远端依赖
T=sparse((1:nb)',f.owner(bf),ones(nb,1),nb,nc);Z=[speye(nc);T];
oldSupport=spones(speye(nc)+(abs(mesh.gradientX)+abs(mesh.gradientR))*Z);
newSupport=spones((abs(Gx)+abs(Gr))*Z);
assert(nnz(newSupport-newSupport.*oldSupport)==0, ...
    'AxisFV:AxisSymmetryGeometry','镜像模板扩大了原一阶邻接，必须重新证明Jacobian支持。');
faceMask=false(nf,1);faceMask(faces)=true;
symmetry=struct('active',true,'kind','axisymmetric-mirror-parity-gradient-v1', ...
    'policy','mirror-parity','faceID',faces,'boundaryIndex',boundaryIndex,'ownerCell',owner, ...
    'faceMask',faceMask,'physicalPointXR',f.centroid(faces,:),'ghostPointXR',ghost, ...
    'ghostDisplacementXR',delta(faces,:),'physicalDisplacementXR',f.delta(faces,:), ...
    'gradientX',Gx,'gradientR',Gr,'cellStencilUnchanged',true, ...
    'physicalFaceGeometryChanged',false,'zeroRadiusRegularizedWithEpsilon',false);
end
