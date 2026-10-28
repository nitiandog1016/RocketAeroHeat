function patches = partitionOrientedTriangles(V,F,grid)
%PARTITIONORIENTEDTRIANGLES 将定向三角曲面按六平面裁切到笛卡尔单元。
% 中文模块说明：本文件只处理几何交集，不建立流体控制体，也不推断流体侧。
% 输入采用 SI 坐标，V 为 N×3 顶点，F 为 M×3 有序三角面索引。
% grid.origin/spacing/cellCounts 定义闭合外域。每个单元的面内部按
% [下界,上界) 归属；最末外边界包含在最后一格。整片三角恰与内部
% 网格面共面时，只归坐标正向单元；这不是流体所有者的判定。
% 普通跨面三角用闭半空间裁切，公共边/点的零面积交集不作为面片。
% 每个非零面片保留原面编号、单元编号、面积、面积质心和定向面积向量。
% 不能用本草稿声称曲面封闭控制体、通流孔口或热流已经求解。

%% 输入检查：拒绝非有限量、非法拓扑、退化三角及外域不足
if ~isnumeric(V) || ~isreal(V) || ~ismatrix(V) || size(V,2)~=3 || isempty(V) || any(~isfinite(V(:)))
    error('CurvedWall:InvalidVertices','V 必须为非空有限实数 N×3 数组。');
end
if ~isnumeric(F) || ~isreal(F) || ~ismatrix(F) || size(F,2)~=3 || isempty(F) || any(~isfinite(F(:))) || ...
        any(F(:)~=fix(F(:))) || any(F(:)<1) || any(F(:)>size(V,1))
    error('CurvedWall:InvalidFaces','F 必须为有效的非空整数 M×3 索引。');
end
if ~isstruct(grid) || ~isscalar(grid) || ~all(isfield(grid,{'origin','spacing','cellCounts'}))
    error('CurvedWall:InvalidGrid','grid 必须包含 origin、spacing、cellCounts。');
end
if ~isnumeric(grid.origin) || ~isvector(grid.origin) || ~isnumeric(grid.spacing) || ~isvector(grid.spacing) || ...
        ~isnumeric(grid.cellCounts) || ~isvector(grid.cellCounts)
    error('CurvedWall:InvalidGrid','网格字段必须为数值向量。');
end
origin = double(grid.origin(:)'); h = double(grid.spacing(:)'); counts = double(grid.cellCounts(:)');
if numel(h)==1, h=repmat(h,1,3); end
if numel(origin)~=3 || numel(h)~=3 || numel(counts)~=3 || ...
        ~isreal(origin) || ~isreal(h) || ~isreal(counts) || ...
        any(~isfinite([origin h counts])) || any(h<=0) || any(counts<1) || any(counts~=fix(counts)) || ...
        prod(counts)>flintmax || any(~isfinite(origin+h.*counts)) || any(origin+h.*counts<=origin)
    error('CurvedWall:InvalidGrid','网格坐标、间距或正整数单元数无效。');
end
V = double(V); F = double(F);
e1 = V(F(:,2),:)-V(F(:,1),:); e2 = V(F(:,3),:)-V(F(:,1),:);
crosses = cross(e1,e2,2); twiceArea = sqrt(sum(crosses.^2,2));
edgeProduct = sqrt(sum(e1.^2,2)).*sqrt(sum(e2.^2,2));
if any(~isfinite(twiceArea)) || any(twiceArea<=64*eps.*edgeProduct) || any(edgeProduct==0)
    error('CurvedWall:DegenerateTriangle','输入含零面积、近共线或不可解析的三角面。');
end
q = (V-origin)./h;
domainTol = 64*eps(max([1 counts]));
if any(any(q < -domainTol | q > counts+domainTol))
    error('CurvedWall:OutsideDomain','完整三角曲面必须位于指定外域内。');
end

%% 面片存储：编号采用 MATLAB 列主序，x 索引变化最快
capacity = max(32,3*size(F,1)); used=0;
area=zeros(capacity,1); centroid=zeros(capacity,3); areaVector=zeros(capacity,3);
originalFaceID=zeros(capacity,1); cellID=zeros(capacity,1); subs=zeros(capacity,3);
polygons=cell(capacity,1); candidateCount=0; zeroAreaIntersections=0;
coplanarTriangleAxisCount=zeros(1,3);

%% 逐原面裁切：仅遍历三角包围盒相交的单元，避免整域笛卡尔积
for faceID=1:size(F,1)
    tri=q(F(faceID,:),:); lower=zeros(1,3); upper=zeros(1,3); locked=false(1,3);
    for axis=1:3
        lo=min(tri(:,axis)); hi=max(tri(:,axis));
        % 共面识别要求原坐标严格相等，不能把一个薄斜三角当作平面。
        constantAxis=all(V(F(faceID,:),axis)==V(F(faceID,1),axis));
        nearest=round(tri(1,axis));
        if constantAxis && abs(tri(1,axis)-nearest)<=domainTol && nearest>=0 && nearest<=counts(axis)
            owner=min(counts(axis),max(1,nearest+1));
            lower(axis)=owner; upper(axis)=owner; locked(axis)=true;
            coplanarTriangleAxisCount(axis)=coplanarTriangleAxisCount(axis)+1;
        else
            lower(axis)=max(1,min(counts(axis),floor(lo)+1));
            upper(axis)=max(1,min(counts(axis),ceil(hi)));
            % 非网格面的常坐标也只归包含它的唯一单元。
            if constantAxis, upper(axis)=lower(axis); end
        end
    end
    originalNormal=crosses(faceID,:)/twiceArea(faceID);
    for iz=lower(3):upper(3)
        for iy=lower(2):upper(2)
            for ix=lower(1):upper(1)
                ijk=[ix iy iz]; poly=tri; candidateCount=candidateCount+1;
                for axis=1:3
                    if locked(axis), continue; end
                    poly=clipPlane(poly,axis,ijk(axis)-1,true);
                    if size(poly,1)<3, break; end
                    poly=clipPlane(poly,axis,ijk(axis),false);
                    if size(poly,1)<3, break; end
                end
                if size(poly,1)<3
                    zeroAreaIntersections=zeroAreaIntersections+1; continue;
                end
                % 在局部物理坐标中计算，避免大平移损伤叉积精度。
                local=poly.*h; anchor=local(1,:); a=0; moment=zeros(1,3); av=zeros(1,3);
                for k=2:size(local,1)-1
                    crossPart=cross(local(k,:)-anchor,local(k+1,:)-anchor)/2;
                    signedPart=dot(crossPart,originalNormal);
                    if signedPart<0 && abs(signedPart)>256*eps*twiceArea(faceID)
                        error('CurvedWall:OrientationLost','裁切产生与原面相反的非零面积。');
                    end
                    part=max(0,signedPart);
                    a=a+part; moment=moment+part*(anchor+local(k,:)+local(k+1,:))/3;
                    av=av+crossPart;
                end
                % 只丢弃严格零面积，不用全局面积门限静默删掉细小面片。
                if a==0
                    zeroAreaIntersections=zeroAreaIntersections+1; continue;
                end
                if ~isfinite(a) || any(~isfinite(moment)) || any(~isfinite(av))
                    error('CurvedWall:NonfinitePatch','裁切面片的几何矩非有限。');
                end
                used=used+1;
                if used>capacity
                    extra=capacity; capacity=capacity+extra;
                    area(capacity,1)=0; centroid(capacity,3)=0; areaVector(capacity,3)=0;
                    originalFaceID(capacity,1)=0; cellID(capacity,1)=0; subs(capacity,3)=0; polygons{capacity,1}=[];
                end
                area(used)=a; centroid(used,:)=origin+moment/a; areaVector(used,:)=av;
                originalFaceID(used)=faceID; subs(used,:)=ijk;
                cellID(used)=ix+(iy-1)*counts(1)+(iz-1)*counts(1)*counts(2);
                polygons{used}=origin+local;
            end
        end
    end
end

%% 输出合同：曲面几何分割与物理/流体控制体求解明确分开
patches=struct('area',area(1:used),'centroid',centroid(1:used,:), ...
    'areaVector',areaVector(1:used,:),'originalFaceID',originalFaceID(1:used), ...
    'cellID',cellID(1:used),'cellSubscripts',subs(1:used,:), ...
    'polygonVertices',{polygons(1:used)},'grid',struct('origin',origin,'spacing',h,'cellCounts',counts));
patches.algorithm='six-plane Sutherland-Hodgman clipping of original oriented triangles';
patches.orientationConvention='original ordered triangle orientation; no fluid-side normal inference or reversal';
patches.ownershipPolicy='[lower,upper); coordinate-positive cell for internal coplanar planes; final boundary belongs to last cell';
patches.fluidOwnerDetermined=false;
patches.physicalSolution=false; patches.eligibleForSurfaceCoverage=false;
patches.fluidVolumesImplemented=false; patches.openFaceAperturesImplemented=false; patches.perCellClosureImplemented=false;
patches.candidateCount=candidateCount; patches.zeroAreaIntersectionCount=zeroAreaIntersections;
patches.coplanarTriangleAxisCount=coplanarTriangleAxisCount;
patches.normalizedOwnershipTolerance=domainTol;
patches.maxAcceptedDomainRoundingExcursion=max([0; -q(:); reshape(q-counts,[],1)]);
end

function out=clipPlane(poly,axis,bound,isLower)
% 中文模块说明：Sutherland–Hodgman 单半空间裁切，保持原顶点环序。
% 交点的被裁轴直接设为平面坐标，防止同一交点的舍入漂移。
if isempty(poly), out=zeros(0,3); return; end
out=zeros(size(poly,1)+2,3); used=0;
prev=poly(end,:); prevD=prev(axis)-bound;
if isLower, prevIn=prevD>=0; else, prevIn=prevD<=0; end
for k=1:size(poly,1)
    curr=poly(k,:); currD=curr(axis)-bound;
    if isLower, currIn=currD>=0; else, currIn=currD<=0; end
    if currIn~=prevIn
        t=prevD/(prevD-currD); point=prev+t*(curr-prev); point(axis)=bound;
        if used==0 || any(point~=out(used,:)), used=used+1; out(used,:)=point; end
    end
    if currIn && (used==0 || any(curr~=out(used,:)))
        used=used+1; out(used,:)=curr;
    end
    prev=curr; prevD=currD; prevIn=currIn;
end
out=out(1:used,:);
if size(out,1)>1 && all(out(1,:)==out(end,:)), out(end,:)=[]; end
end
