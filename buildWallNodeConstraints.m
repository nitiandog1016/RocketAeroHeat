function thermal=buildWallNodeConstraints(volume,faceInput)
% 模块一：把每个壁三角形的给定温度收集到实际共享节点。温界交线保留两端值。
faces=volume.wall;
given=faceInput.givenWallK;
n=size(volume.vertices,1);
assert(height(faceInput)==size(faces,1) && isequal(faceInput.faceID,(1:height(faceInput))') && ...
    iscolumn(given) && all(isfinite(given) & given>0),'RAH:WallNodes:Input','原面温度或编号无效。');
corners=double(faces(:));cornerK=repmat(given,3,1);
incidence=accumarray(corners,1,[n,1]);
lower=accumarray(corners,cornerK,[n,1],@min);
upper=accumarray(corners,cornerK,[n,1],@max);
onWall=incidence>0;
lower(~onWall)=NaN;upper(~onWall)=NaN;
conflict=onWall & lower~=upper;
prescribed=lower;prescribed(conflict)=NaN;
wallIDs=find(onWall);
thermal.wallNodes=table(wallIDs,incidence(wallIDs),lower(wallIDs),upper(wallIDs),prescribed(wallIDs), ...
    'VariableNames',{'vertexID','incidentFaces','minimumGivenK','maximumGivenK','uniqueGivenK'});
% 数值相同的相邻面仍可能具有不同的给温来源，节点表保留全部来源。
sources=string(faceInput.temperatureSource);
assert(~any(ismissing(sources)) && all(strlength(sources)>0),'RAH:WallNodes:Input','原面缺少给温来源。');
origin=strings(n,1);
for source=reshape(unique(sources,'stable'),1,[])
    ids=unique(faces(sources==source,:));
    append=ids(origin(ids)~="");
    origin(append)=origin(append)+" / ";
    origin(ids)=origin(ids)+source;
end
thermal.wallNodes.temperatureSource=categorical(origin(wallIDs));

% 模块二：从相邻四面体中找出与原面相对的唯一节点，按它实际所在的边界分类。
cells=volume.tetrahedra(volume.wallOwner,:);
for corner=1:3
    assert(all(sum(cells==faces(:,corner),2)==1),'RAH:WallNodes:Owner','原面顶点与所列相邻单元不符。');
end
opposite=sum(double(cells),2)-sum(double(faces),2);
assert(all(opposite>=1 & opposite<=n) && all(sum(cells==opposite,2)==1) && ...
    all(faces~=opposite,'all'),'RAH:WallNodes:Owner','不能从相邻体单元恢复唯一对顶点。');
onBox=false(n,1);onBox(volume.outer(:))=true;
assert(~any(onWall & onBox),'RAH:WallNodes:Shell','物理壁和外盒包含相同节点。');
location=ones(size(opposite));
location(onBox(opposite))=2;location(onWall(opposite))=3;
locationNames={'空气内部','外盒','实体壁'};

% 模块三：一次连续节点温度在四个同温壁节点上只能是常数。交线上的多值条件单列。
mixedCount=sum(conflict(cells),2);
code=location;
code(onWall(opposite) & lower(opposite)~=given)=4;
code(mixedCount>0)=5;
names={'内部对点待求温','外盒对点','四壁点同温','四壁点异温','节点温界冲突'};
thermal.faces=faceInput(:,{'faceID','cadFace','areaM2','givenWallK','temperatureSource'});
thermal.faces.ownerCellID=volume.wallOwner;
thermal.faces.oppositeNodeID=opposite;
thermal.faces.oppositeLocation=categorical(location,1:3,locationNames);
thermal.faces.oppositeWallMinK=lower(opposite);
thermal.faces.oppositeWallMaxK=upper(opposite);
thermal.faces.mixedWallNodesInCell=mixedCount;
thermal.faces.nodalConstraint=categorical(code,1:5,names);
end
