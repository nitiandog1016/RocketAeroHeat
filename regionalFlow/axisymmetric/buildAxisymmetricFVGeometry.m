function mesh=buildAxisymmetricFVGeometry(X,R,options)
%BUILDAXISYMMETRICFVGEOMETRY 建立贴体四边形轴对称有限体积几何（每弧度）。
% 模块一：保留输入节点，逐控制体计算平面A、V=integral(r dA)及两种质心。
% 模块二：共享直边只存一次，积分面积向量S=integral(r*n ds)，不是二维边长。
% 模块三：可显式配对i方向平移周期；生成线性精确的最小二乘梯度矩阵。
% 仅支持无折叠凸四边形、r>=0和固定网格；不是表面重网格器。
if nargin<3,options=struct;end
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),{'periodicI'})), ...
    'AxisFV:GeometryOptions','未知几何选项。');
if ~isfield(options,'periodicI'),options.periodicI=false;end
assert(islogical(options.periodicI)&&isscalar(options.periodicI),'AxisFV:GeometryOptions','periodicI须为logical。');
assert(isnumeric(X)&&isnumeric(R)&&isreal(X)&&isreal(R)&&isequal(size(X),size(R))&& ...
    ismatrix(X)&&all(size(X)>=2)&&all(isfinite(X),'all')&&all(isfinite(R)&R>=0,'all'), ...
    'AxisFV:GeometryInput','X/R须为同尺寸有限实数节点矩阵，r非负。');
X=double(X);R=double(R);[ni,nj]=size(X);ids=reshape(1:numel(X),ni,nj);
F=[reshape(ids(1:end-1,1:end-1),[],1),reshape(ids(2:end,1:end-1),[],1), ...
   reshape(ids(2:end,2:end),[],1),reshape(ids(1:end-1,2:end),[],1)];
vertices=[X(:),R(:)];nc=size(F,1);a=zeros(nc,1);v=a;ca=zeros(nc,2);cv=ca;
for c=1:nc
    p=vertices(F(c,:),:);origin=p(1,:);p=p-origin;
    s=p(:,1).*p([2:4 1],2)-p([2:4 1],1).*p(:,2);
    if sum(s)<0,F(c,:)=F(c,[1 4 3 2]);p=p([1 4 3 2],:);end
    q=p([2:4 1],:);edge=q-p;turn=edge(:,1).*edge([2:4 1],2)-edge(:,2).*edge([2:4 1],1);
    assert(all(turn>0),'AxisFV:GeometryFold','只支持严格凸、无折叠、无退化四边形。');
    s=p(:,1).*q(:,2)-q(:,1).*p(:,2);a(c)=sum(s)/2;
    mx=sum((p(:,1)+q(:,1)).*s)/6;mr=sum((p(:,2)+q(:,2)).*s)/6;
    mxr=sum((2*p(:,1).*p(:,2)+p(:,1).*q(:,2)+q(:,1).*p(:,2)+2*q(:,1).*q(:,2)).*s)/24;
    mrr=sum((p(:,2).^2+p(:,2).*q(:,2)+q(:,2).^2).*s)/12;
    v(c)=mr+origin(2)*a(c);ca(c,:)=origin+[mx,mr]/a(c);
    cv(c,:)=origin+[(mxr+origin(2)*mx),(mrr+origin(2)*mr)]/v(c);
end
assert(all(isfinite(v)&v>0)&&all(isfinite(cv),'all'),'AxisFV:GeometryVolume','控制体体积必须有限正值。');

%% 二、从完整拓扑构造共享面，法向指向owner外侧
ea=reshape(F.',[],1);eb=reshape(F(:,[2 3 4 1]).',[],1);
[~,~,groups]=unique(sort([ea,eb],2),'rows');nf=max(groups);
first=accumarray(groups,(1:numel(ea))',[],@min);last=accumarray(groups,(1:numel(ea))',[],@max);
count=accumarray(groups,1);assert(all(count<=2),'AxisFV:GeometryTopology','存在非流形边。');
paired=count==2;
assert(all(ea(first(paired))==eb(last(paired))&eb(first(paired))==ea(last(paired))), ...
    'AxisFV:GeometryTopology','相邻控制体的共享边方向必须相反，拒绝重叠/翻折。');
owner=floor((first-1)/4)+1;neighbour=zeros(nf,1);inside=count==2;neighbour(inside)=floor((last(inside)-1)/4)+1;
start=vertices(ea(first),:);finish=vertices(eb(first),:);delta=finish-start;
lengths=vecnorm(delta,2,2);normal=[delta(:,2),-delta(:,1)]./lengths;
rf=(start(:,2)+finish(:,2))/2;area=rf.*lengths;S=normal.*area;
point=(start+finish)/2;positive=rf>0;
point(positive,:)=((2*start(positive,2)+finish(positive,2)).*start(positive,:)+ ...
    (start(positive,2)+2*finish(positive,2)).*finish(positive,:))./(6*rf(positive));
side=zeros(nf,1);[ia,ja]=ind2sub([ni,nj],ea(first));[ib,jb]=ind2sub([ni,nj],eb(first));
side(~inside&ia==1&ib==1)=1;side(~inside&ia==ni&ib==ni)=2;
side(~inside&ja==1&jb==1)=3;side(~inside&ja==nj&jb==nj)=4;
assert(all(side(~inside)>0),'AxisFV:GeometryTopology','外边界不属于四条逻辑边。');
shift=zeros(nf,2);periodic=false(nf,1);
if options.periodicI
    f1=find(side==1);f2=find(side==2);[~,o]=sort(point(f1,2));f1=f1(o);[~,o]=sort(point(f2,2));f2=f2(o);
    d=point(f1,:)-point(f2,:);scale=max(1,max(abs(vertices),[],'all'));
    assert(numel(f1)==numel(f2)&&max(abs(d-d(1,:)),[],'all')<1e-12*scale&& ...
        abs(d(1,2))<1e-12*scale&&max(abs(area(f1)-area(f2)))<1e-12*scale^2&& ...
        max(abs(normal(f1,:)+normal(f2,:)),[],'all')<1e-12, ...
        'AxisFV:PeriodicGeometry','i方向周期边必须由单一轴向平移精确对应。');
    neighbour(f1)=owner(f2);shift(f1,:)=d;periodic(f1)=true;side(f1)=0;
    keep=true(nf,1);keep(f2)=false;
    owner=owner(keep);neighbour=neighbour(keep);point=point(keep,:);normal=normal(keep,:);
    lengths=lengths(keep);area=area(keep);S=S(keep,:);side=side(keep);shift=shift(keep,:);periodic=periodic(keep);
end
nf=numel(owner);inside=neighbour>0;idx=(1:nf)';
C=sparse([owner;neighbour(inside)],[idx;idx(inside)],[ones(nf,1);-ones(nnz(inside),1)],nc,nf);
B=sparse([owner;neighbour(inside)],[idx;idx(inside)],ones(nf+nnz(inside),1),nc,nf);
closure=C*S-[zeros(nc,1),a];
assert(max(abs(closure)./max(a,realmin),[],'all')<1e-9,'AxisFV:MetricIdentity','轴对称几何度量恒等式失败。');
dn=zeros(nf,2);dn(inside,:)=cv(neighbour(inside),:)+shift(inside,:)-cv(owner(inside),:);
dn(~inside,:)=point(~inside,:)-cv(owner(~inside),:);
normalDistance=sum(dn.*normal,2);
assert(all(normalDistance>0),'AxisFV:GeometryDistance','共享面/边界位移须指向外侧。');

%% 三、以体积质心及边界r加权采样点构造线性精确WLS；周期邻居带平移
boundary=find(~inside);nb=numel(boundary);boundaryColumn=zeros(nf,1);boundaryColumn(boundary)=nc+(1:nb);
sample=neighbour;sample(boundary)=boundaryColumn(boundary);
cellID=[owner;neighbour(inside)];sampleID=[sample;owner(inside)];d=[dn;-dn(inside,:)];
w=1./sum(d.^2,2);ax=accumarray(cellID,w.*d(:,1).^2,[nc 1]);ar=accumarray(cellID,w.*d(:,2).^2,[nc 1]);
cross=accumarray(cellID,w.*d(:,1).*d(:,2),[nc 1]);detA=ax.*ar-cross.^2;
assert(all(detA>1e-13.*ax.*ar),'AxisFV:GradientRank','梯度约束病态或欠秩。');
cx=w.*(ar(cellID).*d(:,1)-cross(cellID).*d(:,2))./detA(cellID);
cr=w.*(ax(cellID).*d(:,2)-cross(cellID).*d(:,1))./detA(cellID);
Gx=sparse([cellID;cellID],[sampleID;cellID],[cx;-cx],nc,nc+nb);
Gr=sparse([cellID;cellID],[sampleID;cellID],[cr;-cr],nc,nc+nb);
mesh=struct('kind','axisymmetric-body-fitted-fv-v1','X',X,'R',R,'vertices',vertices,'cells',F, ...
    'shape',[ni-1 nj-1],'cellCount',nc,'planarArea',a,'planarCentroid',ca, ...
    'volumePerRadian',v,'volumeCentroid',cv,'faces',struct('owner',owner,'neighbour',neighbour, ...
    'centroid',point,'normal',normal,'length',lengths,'areaPerRadian',area,'areaVectorPerRadian',S, ...
    'boundarySide',side,'neighbourShift',shift,'periodic',periodic,'delta',dn,'normalDistance',normalDistance), ...
    'boundaryFaces',boundary,'incidence',C,'cellFaceMultiplicity',B,'gradientX',Gx,'gradientR',Gr, ...
    'periodicI',options.periodicI,'metricClosureRelative',max(abs(closure)./max(a,realmin),[],'all'), ...
    'units','SI; volumes and surface areas per radian; full revolution multiply2*pi', ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false);
end
