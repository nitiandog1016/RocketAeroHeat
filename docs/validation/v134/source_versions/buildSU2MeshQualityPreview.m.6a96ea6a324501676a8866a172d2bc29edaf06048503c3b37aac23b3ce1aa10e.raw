function [mesh,info]=buildSU2MeshQualityPreview(volume)
%BUILDSU2MESHQUALITYPREVIEW 独立计算四面体SICN，显示低质量单元及邻近原壁面。
% 模块一：从实际体连接计算等边参考Jacobian，不依赖历史质量摘要。
% 模块二：取SICN<0.1单元的外露三角面，逐面保留体单元编号与原坐标。
% 模块三：灰色原壁面只作几何定位，不赋予体质量或物理场量。
cellIDs=find(volume.cellTypes==10);T=vertcat(volume.cells{cellIDs});P=volume.vertices;
a=P(T(:,2),:)-P(T(:,1),:);b=(P(T(:,3),:)-P(T(:,1),:)-a/2)*2/sqrt(3);
c=(P(T(:,4),:)-P(T(:,1),:)-a/2-b*sqrt(3)/6)/sqrt(2/3);
bc=cross(b,c,2);ca=cross(c,a,2);ab=cross(a,b,2);det=dot(a,bc,2);
quality=3*det./sqrt((sum(a.^2+b.^2+c.^2,2)).*(sum(bc.^2+ca.^2+ab.^2,2)));
assert(all(isfinite(quality)) && all(quality>0) && all(quality<=1+1e-12));
bad=find(quality<.1);assert(~isempty(bad),'没有SICN<0.1的四面体。');
selected=T(bad,:);tris=[selected(:,[1 3 2]);selected(:,[1 2 4]);selected(:,[2 3 4]);selected(:,[3 1 4])];
owners=repmat(bad,4,1);[~,~,group]=unique(sort(tris,2),'rows');count=accumarray(group,1);
assert(all(count<=2));keep=count(group)==1;tris=tris(keep,:);owners=owners(keep);
points=P(unique(selected(:)),:);extent=max(points)-min(points);padding=max(.01,.5*max(extent));
if isfield(volume.meshConstruction,'normalHeightsM')
    padding=max(padding,1.5*max(volume.meshConstruction.normalHeightsM));
end
low=min(points)-padding;high=max(points)+padding;
SF=volume.sourceFaces;SV=volume.sourceVertices;centers=(SV(SF(:,1),:)+SV(SF(:,2),:)+SV(SF(:,3),:))/3;
context=find(all(centers>=low & centers<=high,2));
allTriangles=[tris;SF(context,:)];[nodeIDs,~,mapped]=unique(allTriangles(:));faces=reshape(mapped,size(allTriangles));
mesh=struct('V',P(nodeIDs,:),'F',faces);
values=[quality(owners);nan(numel(context),1)];sourceCell=[cellIDs(owners);zeros(numel(context),1)];
info=struct('method','SICN = 3 det(J) / (norm(J,F) norm(cofactor(J),F))', ...
    'threshold',.1,'minimumSICN',min(quality),'lowQualityCells',numel(bad), ...
    'sourceCellIDs',cellIDs(bad),'sourceTetRows',bad,'cellSICN',quality(bad), ...
    'cellCenters',(P(selected(:,1),:)+P(selected(:,2),:)+P(selected(:,3),:)+P(selected(:,4),:))/4, ...
    'sourceVolumeNodeIDs',nodeIDs,'displaySourceCellIDs',sourceCell, ...
    'displayValues',values,'contextOriginalFaceIDs',context,'visibleCellFaces',size(tris,1), ...
    'geometryOnly',true,'physicsQualified',false,'newQualifiedArea',0);
end
