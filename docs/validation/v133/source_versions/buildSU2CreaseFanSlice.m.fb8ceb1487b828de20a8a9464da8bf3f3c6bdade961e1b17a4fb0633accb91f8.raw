function [section,info]=buildSU2CreaseFanSlice(volume,parentID)
%BUILDSU2CREASEFANSLICE 从实际棱柱连接提取尾缘一列截面，不生成虚构压力或温度。
% 模块一：将原壁点的侧/底两分支编号对应到真实体单元。
% 模块二：提取该列的三角横截面、去除相邻环段的重复面，保留原节点来源。
r=volume.meshConstruction;
assert(isfield(r,'constructionOperation') && startsWith(r.constructionOperation,'closed-convex-rim-'));
if nargin<2
    parentID=r.rimNodeIDs(1);
    if any(r.rimNodeIDs==17),parentID=17;end
end
column=find(r.rimNodeIDs==parentID);assert(isscalar(column),'请选择当前棱边中的原壁点编号。');
eligible=unique([r.sideNodeIDs(:,column);r.baseNodeIDs(:,column)])+1;
W=vertcat(volume.cells{r.originalPrisms+1:r.prisms});
tri=[W(:,1:3);W(:,4:6)];tri=tri(all(ismember(tri,eligible),2),:);
[~,keep]=unique(sort(tri,2),'rows','stable');tri=tri(keep,:);
assert(~isempty(tri),'实际体单元中未找到该列的棱柱横截面。');
[ids,~,map]=unique(tri(:));F=reshape(map,size(tri));V=volume.vertices(ids,:);
area=vecnorm(cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2),2,2)/2;
assert(all(area>0) && numel(ids)==numel(eligible));
section=struct('V',V,'F',F);
info=struct('sourceVertexID',parentID+1,'sourceVertexIDZeroBased',parentID, ...
    'sourceVolumeNodeIDs',ids,'sourceVolumeTriangles',tri,'faces',size(F,1), ...
    'area',sum(area),'geometryOnly',true,'physicsQualified',false,'newQualifiedArea',0);
end
