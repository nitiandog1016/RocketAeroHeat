function [input,wall] = buildRFZSurfaceInput(mesh,geometry,mode,cavityK)
% 模块一：绑定完整原面、已确认的刚体朝向和来源。所有温度都是输入量。
reference=mesh.rfzBoundary;
regions=reference.regions;
assert(geometry.isAligned && geometry.scale==1 && isequal(regions.cadFace,mesh.rfzRegions.cadFace), ...
    'RAH:RFZ:Input','需要已确认的米制 RFZ 原壁。');
ids=regions.cadFace;
counts=regions.lastFace-regions.firstFace+1;
assert(isequal(ids,(1:height(regions))') && isequal(regions.firstFace,[1;regions.lastFace(1:end-1)+1]) && ...
    all(counts>0) && regions.lastFace(end)==size(mesh.F,1),'RAH:RFZ:Input','原面编号区间不连续。');
cadFace=repelem(ids,counts);
area=mesh.area(:);
assert(all(isfinite(area) & area>0) && max(abs(accumarray(cadFace,area)-regions.stlAreaM2))<1e-8, ...
    'RAH:RFZ:Input','当前网格分区面积与原面表不符。');

% 模块二：按区间展开边界输入，保留封盖内侧的全部实体面。
exteriorK=reference.condition.publishedCondition.prescribedIsothermalWallK;
mode=validatestring(mode,{'matchExterior','specified'});
if strcmp(mode,'matchExterior')
    cavityK=exteriorK;
else
    validateattributes(cavityK,{'numeric'},{'scalar','real','finite','positive'});
end
isCavity=regions.referenceRole(cadFace)=="封盖内侧 CAD 壁";
temperature=repmat(exteriorK,numel(cadFace),1);
temperature(isCavity)=cavityK;
pose=[geometry.rotationMatrix(:);geometry.rotationCenter(:)];
input=struct('schema','rfz-rec-surface-input-v1','domain','fullCAD','dataset',reference.condition.dataset, ...
    'stlSHA256',reference.stlSHA256,'cadTableSHA256',reference.cadTableSHA256,'referenceCaseSHA256',reference.referenceCaseSHA256, ...
    'meshSHA256',phaseAwareMeshIdentity(mesh),'poseDecimal',{cellstr(compose('%.17g',pose))}, ...
    'faces',numel(cadFace),'cavityFaces',nnz(isCavity),'exteriorWallK',exteriorK,'cavityMode',mode,'cavityWallK',cavityK);
input.wallSHA256=phaseAwarePayloadSHA([input.domain '/' mode],cadFace,temperature,double(isCavity));
if nargout<2,return;end
source=repmat("论文给定外壁温度",numel(cadFace),1);
if strcmp(mode,'matchExterior'),source(isCavity)="项目假设：内壁随外壁同温";else,source(isCavity)="用户单独给定内壁温度";end
wall=table((1:numel(cadFace))',cadFace,area,temperature,source,'VariableNames', ...
    {'faceID','cadFace','areaM2','givenWallK','temperatureSource'});
end
