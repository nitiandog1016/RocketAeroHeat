function c=createTmrAxisymmetricOutletCase(parentEvidenceFile)
%CREATETMRAXISYMMETRICOUTLETCASE 冻结实际出口剖面的明确HLLE外储库对照。
% 模块一：保留v66网格/Q/气体/入口/外域/墙，仅替换iMax数值储库为冻结owner剖面。
% 模块二：四原始量和SA nuTilde均逐源face给定；不是自适应零梯度或外部CFD输入。
% 模块三：只接受原出口面采样点，不插值、不外推；剖面来源及原始值完整保存。
% 模块四：用于边界影响诊断，不作为外部真实压力或实验拟合，不增加生产面积。
if nargin<1,parentEvidenceFile='';end
c=createTmrAxisymmetricSensitivityCase('baseline',parentEvidenceFile);
mesh=c.mesh;info=axisymmetricFVState(c.Q(:,1:4),c.gas);
ids=find(mesh.faces.boundarySide==2);owner=mesh.faces.owner(ids);
points=mesh.faces.centroid(ids,:);values=info(owner,1:4);N=c.Q(owner,5)./c.Q(owner,1);
profile=struct('kind','frozen-v66-owner-outlet-reservoir-sensitivity-v1','faceID',ids,'ownerCell',owner, ...
    'centroidXR',points,'primitiveRhoUxUrT',values,'nuTildeM2s',N, ...
    'source',c.parentBinding,'valueKind','prescribed numerical sensitivity reservoir, not measured exterior truth', ...
    'dynamicZeroGradient',false,'interpolated',false,'extrapolated',false,'experimentalValidation',false);
referenceBoundary=c.boundary.iMax;
c.boundary.iMax.primitive=@(point)sample(point,points,values);
c.boundary.iMax.nuTilde=@(point)sample(point,points,N);
c.boundaryStudy=struct('kind','outlet-profile','changedBoundary','iMax', ...
    'referenceBoundary',referenceBoundary,'profile',profile, ...
    'purpose','diagnose incoming acoustic/reservoir influence; not choose pressure to fit experimental points', ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0);
c.gridStudy.kind='outlet-profile';
end

function value=sample(point,points,values)
% 模块用途：显式按原面坐标取冻结值，拒绝换网格或任意空间查询造成的隐式外推。
assert(isnumeric(point)&&isreal(point)&&isequal(size(point),[1 2])&&all(isfinite(point)), ...
    'AxisFV:OutletProfilePoint','出口剖面采样需要一个有限实坐标。');
distance=vecnorm(points-point,2,2);[minimum,index]=min(distance);
scale=max(abs([point(:);points(index,:)']));
assert(minimum<=64*eps(scale)&&nnz(distance<=64*eps(scale))==1, ...
    'AxisFV:OutletProfilePoint','只能在原出口面的明确采样点给定值，不插值或外推。');
value=values(index,:);
end
