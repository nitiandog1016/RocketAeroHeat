function projection=buildSU2WallProjection(volume,surface)
%BUILDSU2WALLPROJECTION 将三维壁节点双网格分片保守投影回全部原三角面。
% 模块一：每个壁多边形按顶点/边中点/面中心分成节点面积片，采用分片常值。
% 模块二：在原三角面平面中裁剪交集；非共面候选不投影，不最近面补漏。
% 模块三：逐原面面积与逐源节点面积双向核对，保证常数场及有符号积分守恒。
% 该表示为节点面积片常值，不宣称高阶点值精度；几何映射不等于物理解资格。
assert(isequal(surface.V,volume.sourceVertices)&&isequal(surface.F,volume.sourceFaces), ...
    'RocketAeroHeat:SU2:ProjectionSource','必须是建体网格时完整的同一原表面。');
V=volume.vertices;nf=size(surface.F,1);nv=size(V,1);scale=max(1,max(abs(V(:))));tol=1e-10*scale;
if isfield(volume,'wallSourceFaceID')
    % 原三角壁面保持时，三个中位对偶片各占面面积1/3，无需全环带裁剪搜索。
    F=double(surface.F);assert(isequal(volume.wallSourceFaceID,(1:nf)') && ...
        isequal(vertcat(volume.bodyWall{:}),F) && isequal(V(1:size(surface.V,1),:),surface.V));
    av=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2)/2;areas=vecnorm(av,2,2);
    areaMap=sparse(repmat((1:nf)',3,1),F(:),repmat(areas/3,3,1),nf,nv);
    nodeArea=full(sum(areaMap,1))';nodeAreaVector=zeros(nv,3);
    for j=1:3,nodeAreaVector(:,j)=accumarray(F(:),repmat(av(:,j)/3,3,1),[nv 1]);end
    projection=struct('kind','original-triangle-conservative-median-dual-constant-v1', ...
        'weights',spdiags(1./areas,0,nf,nf)*areaMap,'sourceNodeArea',nodeArea,'sourceNodeAreaVector',nodeAreaVector, ...
        'engineNodeArea',vecnorm(nodeAreaVector,2,2),'targetFaceArea',areas,'wallVertexIDs',find(nodeArea>0), ...
        'sourceVertices',V,'sourceCells',{volume.bodyWall},'targetVertices',surface.V,'targetFaces',F, ...
        'maximumTargetAreaRelativeError',max(abs(full(sum(areaMap,2))-areas)./areas), ...
        'maximumSourceAreaRelativeError',0,'representation','exact original triangle median-dual patches; three equal-area vertex weights', ...
        'engineAreaRule','SU2 integrated HeatFlux uses magnitude of summed nodal area vectors; differs from scalar patch area at curved edges', ...
        'productionQualified',false,'newQualifiedRocketArea',0);
    return;
end
areas=zeros(nf,1);nodeArea=zeros(nv,1);nodeAreaVector=zeros(nv,3);patches={};nodes=[];bands=[];
for k=1:numel(volume.bodyWall)
    f=volume.bodyWall{k};P=V(f,:);center=mean(P,1);n=numel(f);
    for j=1:n
        previous=mod(j-2,n)+1;next=mod(j,n)+1;
        poly=[P(j,:);(P(j,:)+P(next,:))/2;center;(P(previous,:)+P(j,:))/2];
        patches{end+1,1}=poly;nodes(end+1,1)=f(j);bands(end+1,1)=volume.wallSourceBand(k); %#ok<AGROW>
        nodeArea(f(j))=nodeArea(f(j))+norm(cross(poly(2,:)-poly(1,:),poly(3,:)-poly(1,:)))/2+ ...
            norm(cross(poly(3,:)-poly(1,:),poly(4,:)-poly(1,:)))/2;
        nodeAreaVector(f(j),:)=nodeAreaVector(f(j),:)+(cross(poly(2,:)-poly(1,:),poly(3,:)-poly(1,:))+ ...
            cross(poly(3,:)-poly(1,:),poly(4,:)-poly(1,:)))/2;
    end
end
rows=[];columns=[];overlap=[];
for i=1:nf
    tri=surface.V(surface.F(i,:),:);a=tri(1,:);e1=tri(2,:)-a;e1=e1/norm(e1);
    normal=cross(tri(2,:)-a,tri(3,:)-a);areas(i)=norm(normal)/2;normal=normal/norm(normal);e2=cross(normal,e1);
    basis=[e1;e2]';target=(tri-a)*basis;
    candidates=find(bands==surface.param.faceBandID(i));
    for j=candidates'
        poly=patches{j};
        if max(abs((poly-a)*normal'))>tol,continue;end
        plane=(poly-a)*basis;
        if any(max(plane,[],1)<min(target,[],1)-tol) || any(min(plane,[],1)>max(target,[],1)+tol),continue;end
        clipped=clipTriangle(plane,target);if size(clipped,1)<3,continue;end
        area=abs(sum(clipped(:,1).*clipped([2:end 1],2)-clipped(:,2).*clipped([2:end 1],1)))/2;
        if area<=eps(scale)^2,continue;end
        rows(end+1,1)=i;columns(end+1,1)=nodes(j);overlap(end+1,1)=area; %#ok<AGROW>
    end
end
areaMap=sparse(rows,columns,overlap,nf,nv);rowError=max(abs(full(sum(areaMap,2))-areas)./areas);
active=nodeArea>0;colError=max(abs(full(sum(areaMap(:,active),1))'-nodeArea(active))./nodeArea(active));
assert(rowError<2e-8 && colError<2e-8,'RocketAeroHeat:SU2:ProjectionCoverage', ...
    '保守投影面积不闭合：原面%.3g、源节点%.3g；禁止补最近面。',rowError,colError);
weights=spdiags(1./areas,0,nf,nf)*areaMap;
projection=struct('kind','original-triangle-conservative-median-dual-constant-v1', ...
    'weights',weights,'sourceNodeArea',nodeArea,'sourceNodeAreaVector',nodeAreaVector, ...
    'engineNodeArea',vecnorm(nodeAreaVector,2,2),'targetFaceArea',areas,'wallVertexIDs',find(active), ...
    'sourceVertices',V,'sourceCells',{volume.bodyWall},'targetVertices',surface.V,'targetFaces',surface.F, ...
    'maximumTargetAreaRelativeError',rowError,'maximumSourceAreaRelativeError',colError, ...
    'representation','piecewise constant on boundary median-dual node patches; area-mean target faces', ...
    'engineAreaRule','SU2 integrated HeatFlux uses magnitude of summed nodal area vectors; differs from scalar patch area at curved edges', ...
    'productionQualified',false,'newQualifiedRocketArea',0);
end

function polygon=clipTriangle(polygon,triangle)
% Sutherland-Hodgman平面凸裁剪，交点由真实线段计算，不移动原网格顶点。
for k=1:3
    A=triangle(k,:);B=triangle(mod(k,3)+1,:);e=B-A;
    if isempty(polygon),return;end
    result=zeros(0,2);P=polygon(end,:);dp=e(1)*(P(2)-A(2))-e(2)*(P(1)-A(1));
    for j=1:size(polygon,1)
        Q=polygon(j,:);dq=e(1)*(Q(2)-A(2))-e(2)*(Q(1)-A(1));
        if (dp>=0)~=(dq>=0),result(end+1,:)=P+dp/(dp-dq)*(Q-P);end %#ok<AGROW>
        if dq>=0,result(end+1,:)=Q;end %#ok<AGROW>
        P=Q;dp=dq;
    end
    polygon=result;
end
end
