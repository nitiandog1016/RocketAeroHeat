function mapping=intersectShockNormalsWithBody(V,F,positions,directions)
%INTERSECTSHOCKNORMALSWITHBODY 将共同激波的内法线精确关联到原三角面元。
% 模块一：验证闭合、一致朝外的原始三角拓扑；不重网格或平滑真实特征。
% 模块二：逐射线Moller-Trumbore求最近正向交点，保留面号、重心权和并列面。
% 模块三：内部起点、离体射线、掠切均保持NaN；几何命中不代表激波层匹配。
validateattributes(V,{'double'},{'real','finite','2d','ncols',3,'nonempty'});
validateattributes(F,{'double'},{'real','finite','integer','positive','2d','ncols',3,'nonempty'});
validateattributes(positions,{'double'},{'real','finite','2d','ncols',3,'nonempty'});
validateattributes(directions,{'double'},{'real','finite','size',size(positions)});
assert(max(F,[],'all')<=size(V,1),'ShockLayer:BodyMesh','原面元含越界节点。');
assert(max(abs(vecnorm(directions,2,2)-1))<1e-10,'ShockLayer:Ray','激波内法线须为单位向量，不自动归一。');
edges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];[~,~,group]=unique(sort(edges,2),'rows');
counts=accumarray(group,1);orientation=accumarray(group,2*(edges(:,1)<edges(:,2))-1);
assert(all(counts==2)&&all(orientation==0),'ShockLayer:BodyMesh','要求闭合、无非流形边且面元定向一致的原几何。');
v0=V(F(:,1),:);e1=V(F(:,2),:)-v0;e2=V(F(:,3),:)-v0;
normal=cross(e1,e2,2);twiceArea=vecnorm(normal,2,2);
assert(all(twiceArea>0),'ShockLayer:BodyMesh','原几何含退化三角形。');
volume=sum(dot(v0,cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
assert(volume>0,'ShockLayer:BodyMesh','原几何整体朝向不是朝外。');
normal=normal./twiceArea;span=norm(max(V,[],1)-min(V,[],1));
distanceToleranceM=1e-9*span;barycentricTolerance=1e-9;count=size(positions,1);
hit=false(count,1);faceID=nan(count,1);distance=faceID;barycentric=nan(count,3);point=barycentric;
normalAtBody=barycentric;coincidentFaces=cell(count,1);reason=repmat({''},count,1);

%% 二、只取最近的进入几何交点，不把出体点或另一个壳体当作壁面入口
for j=1:count
    D=directions(j,:);O=positions(j,:);
    P=cross(repmat(D,size(F,1),1),e2,2);det=dot(e1,P,2);
    usable=abs(det)>64*eps*twiceArea;
    relative=O-v0;u=nan(size(det));v=u;t=u;
    u(usable)=dot(relative(usable,:),P(usable,:),2)./det(usable);
    Q=cross(relative,e1,2);
    v(usable)=Q(usable,:)*D'./det(usable);
    t(usable)=dot(e2(usable,:),Q(usable,:),2)./det(usable);
    good=usable & u>=-barycentricTolerance & v>=-barycentricTolerance & ...
        u+v<=1+barycentricTolerance & t>128*eps*span;
    if ~any(good)
        reason{j}='内法线未与原闭合箭体相交；不能外推壁面或补编号';continue;
    end
    nearest=min(t(good));ids=find(good & abs(t-nearest)<=distanceToleranceM);
    chosen=min(ids);coincidentFaces{j}=ids;
    if dot(normal(chosen,:),D)>=-1e-10
        reason{j}='最近交点为出体或掠切：起点可能在几何内部，禁止当作进入壁面的映射';continue;
    end
    weights=[1-u(chosen)-v(chosen),u(chosen),v(chosen)];
    actual=weights*V(F(chosen,:),:);rayPoint=O+t(chosen)*D;
    assert(norm(actual-rayPoint)<distanceToleranceM,'ShockLayer:BodyHit','重心重建与真实射线交点不一致。');
    hit(j)=true;faceID(j)=chosen;distance(j)=t(chosen);barycentric(j,:)=weights;
    point(j,:)=actual;normalAtBody(j,:)=normal(chosen,:);
end
mapping=struct('method','原三角面Moller-Trumbore最近进入交点；保留原编号/原节点', ...
    'hit',hit,'originalFaceID',faceID,'distanceM',distance,'bodyPosition',point, ...
    'barycentric',barycentric,'bodyNormal',normalAtBody,'coincidentOriginalFaces',{coincidentFaces}, ...
    'missingReason',{reason},'rayPosition',positions,'rayDirection',directions, ...
    'distanceToleranceM',distanceToleranceM,'barycentricTolerance',barycentricTolerance, ...
    'inputFaceCount',size(F,1),'signedBodyVolumeM3',volume, ...
    'shockLayerThicknessSolved',false,'bodyMatchedByFlowSolution',false, ...
    'flowValuesAssignedToBody',false,'productionQualified',false,'newProductionSolvedAreaM2',0);
end
