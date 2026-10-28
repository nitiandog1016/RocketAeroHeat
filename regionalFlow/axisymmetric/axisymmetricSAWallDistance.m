function result=axisymmetricSAWallDistance(mesh,boundary)
%AXISYMMETRICSAWALLDISTANCE 控制体真体积质心到全部真实黏性壁线段的最短距。
% 模块一：按四条逻辑边的wall标签收集壁面，不把入口/出口虚构为端盖。
% 模块二：从owner单元的真实顶点边恢复壁段，核对r加权面质心和长度。
% 模块三：逐段投影到闭区间[0,1]，同时比较端点；保存最近真实face与坐标。
% 旋转子午壁的三维最短距严格等于同方位角XR平面最短距，限r>=0无旋几何。
% 本函数不使用最近墙节点、墙cell中心、r-R(X)或网格线长度替代最短壁距。
% 周期网格仅支持平直轴向墙；变壁形需周期镜像的完整搜索，当前明确拒绝。
names={'iMin','iMax','jMin','jMax'};
assert(isstruct(mesh)&&isscalar(mesh)&&all(isfield(mesh, ...
    {'volumeCentroid','cellCount','vertices','cells','faces','boundaryFaces'}))&& ...
    isstruct(boundary)&&isscalar(boundary)&&all(isfield(boundary,names)), ...
    'AxisFV:SAWallDistance','壁距需要完整有限体积几何和四边边界。');
f=mesh.faces;
assert(isstruct(f)&&all(isfield(f,{'owner','neighbour','boundarySide','centroid','length'})), ...
    'AxisFV:SAWallDistance','缺少共享面几何。');
points=mesh.volumeCentroid;vertices=mesh.vertices;cells=mesh.cells;
assert(isnumeric(points)&&isreal(points)&&isequal(size(points),[mesh.cellCount 2])&& ...
    all(isfinite(points),'all')&&all(points(:,2)>=0)&& ...
    isnumeric(vertices)&&isreal(vertices)&&size(vertices,2)==2&&all(isfinite(vertices),'all')&& ...
    all(vertices(:,2)>=0)&&isnumeric(cells)&&isequal(size(cells),[mesh.cellCount 4])&& ...
    all(isfinite(cells)&cells>=1&cells<=size(vertices,1)&cells==fix(cells),'all'), ...
    'AxisFV:SAWallDistance','壁距只接受真实有限体积质心、顶点和四边形拓扑。');
wallSides=false(4,1);
for s=1:4
    bc=boundary.(names{s});
    assert(isstruct(bc)&&isscalar(bc)&&isfield(bc,'type')&& ...
        ((ischar(bc.type)&&isrow(bc.type))||(isstring(bc.type)&&isscalar(bc.type))), ...
        'AxisFV:SAWallDistance','边界须明确type。');
    wallSides(s)=strcmp(bc.type,'wall');
end
nf=numel(f.owner);bf=mesh.boundaryFaces;
assert(isnumeric(bf)&&isvector(bf)&&all(isfinite(bf)&bf>=1&bf<=nf&bf==fix(bf))&& ...
    numel(unique(bf))==numel(bf)&&all(f.neighbour(bf)==0)&& ...
    numel(f.neighbour)==nf&&numel(f.boundarySide)==nf&&numel(f.length)==nf&& ...
    isequal(size(f.centroid),[nf 2])&&all(ismember(f.boundarySide(bf),1:4)), ...
    'AxisFV:SAWallDistance','边界面索引或标签无效。');
wallFaceID=sort(bf(wallSides(f.boundarySide(bf))));
assert(~isempty(wallFaceID),'AxisFV:SAWallDistance','没有真实wall线段，不能补壁距。');
wallStart=zeros(numel(wallFaceID),2);wallEnd=wallStart;

%% 二、独立恢复每个owner的四条顶点边，保留真实墙段与源face绑定
for j=1:numel(wallFaceID)
    face=wallFaceID(j);owner=f.owner(face);
    assert(isfinite(owner)&&owner>=1&&owner<=mesh.cellCount&&owner==fix(owner)&& ...
        all(isfinite(f.centroid(face,:)))&&isfinite(f.length(face))&&f.length(face)>0, ...
        'AxisFV:SAWallDistance','wall面owner、质心或长度无效。');
    ids=cells(owner,:);a=vertices(ids,:);b=vertices(ids([2 3 4 1]),:);
    radius=(a(:,2)+b(:,2))/2;centroid=(a+b)/2;positive=radius>0;
    centroid(positive,:)=((2*a(positive,2)+b(positive,2)).*a(positive,:)+ ...
        (a(positive,2)+2*b(positive,2)).*b(positive,:))./(6*radius(positive));
    edgeLength=vecnorm(b-a,2,2);scale=max(abs([a;b;f.centroid(face,:)]),[],'all');
    tolerance=256*eps(max(scale,realmin));
    match=find(vecnorm(centroid-f.centroid(face,:),2,2)<=tolerance& ...
        abs(edgeLength-f.length(face))<=tolerance);
    assert(isscalar(match)&&edgeLength(match)>0,'AxisFV:SAWallDistance', ...
        'wall面不唯一匹配owner真实顶点边，拒绝过时或篡改几何。');
    wallStart(j,:)=a(match,:);wallEnd(j,:)=b(match,:);
end

%% 周期范围：只有平直轴向墙的有限段并集与无限周期墙给出相同最短距离
periodic=false;
if isfield(mesh,'periodicI')
    assert(islogical(mesh.periodicI)&&isscalar(mesh.periodicI),'AxisFV:SAWallDistance','periodicI须明确为logical。');
    periodic=mesh.periodicI;
end
if periodic
    for side=1:4
        ids=find(f.boundarySide(wallFaceID)==side);
        if isempty(ids),continue;end
        nodes=[wallStart(ids,:);wallEnd(ids,:)];scale=max(abs([nodes;points]),[],'all');
        tolerance=256*eps(max(scale,realmin));
        assert(side>=3&&max(nodes(:,2))-min(nodes(:,2))<=tolerance&& ...
            all(points(:,1)>=min(nodes(:,1))-tolerance&points(:,1)<=max(nodes(:,1))+tolerance), ...
            'AxisFV:SAWallDistance','周期壁距当前只支持覆盖所有质心X的平直轴向墙；变壁形需真实周期镜像，未替代。');
    end
end

%% 三、投影闭线段；并列最近段按最小faceID确定，不添加虚构入口/出口墙
n=mesh.cellCount;distance2=inf(n,1);faceID=zeros(n,1);point=zeros(n,2);
projectionFraction=zeros(n,1);rawProjectionFraction=zeros(n,1);
for j=1:numel(wallFaceID)
    a=wallStart(j,:);edge=wallEnd(j,:)-a;length2=sum(edge.^2);
    assert(isfinite(length2)&&length2>0,'AxisFV:SAWallDistance','wall线段长度平方无效。');
    raw=((points-a)*edge')/length2;fraction=min(max(raw,0),1);candidate=a+fraction.*edge;
    candidate2=sum((points-candidate).^2,2);better=candidate2<distance2;
    distance2(better)=candidate2(better);faceID(better)=wallFaceID(j);point(better,:)=candidate(better,:);
    projectionFraction(better)=fraction(better);rawProjectionFraction(better)=raw(better);
end
distance=sqrt(distance2);
assert(all(isfinite(distance)&distance>0)&&all(faceID>0)&&all(isfinite(point),'all'), ...
    'AxisFV:SAWallDistance','控制体壁距必须有限严格为正，不用epsilon代替墙上点。');
result=struct('distance',distance,'faceID',faceID,'point',point, ...
    'projectionFraction',projectionFraction,'rawProjectionFraction',rawProjectionFraction, ...
    'wallFaceID',wallFaceID,'wallStart',wallStart,'wallEnd',wallEnd, ...
    'queryPoint',points,'queryPointKind','true-r-weighted-volume-centroid', ...
    'method','Euclidean minimum over all actual closed meridional wall segments; no end caps', ...
    'periodicGeometry',periodic,'periodicWallPolicy','straight-axial-walls-only; varying periodic walls rejected', ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false);
end
