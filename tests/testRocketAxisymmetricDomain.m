function report=testRocketAxisymmetricDomain(folder)
%TESTROCKETAXISYMMETRICDOMAIN 完整真实箭体C域的几何、积分、绑定及拒绝校核。
% 模块一：保留原鼻尖/接缝/底面和全部原面band，不把几何通过计成物理解。
% 模块二：独立三角Duffy/Gauss求r权积分，不复用构造器的多边形矩公式。
% 模块三：独立圆台账、轴/墙/外场拓扑、仿射梯度及实际首层距离检查。
% 模块四：层数/首间距/外域可独立改变，非星形核、突起和损坏源必须拒绝。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'regionalFlow','axisymmetric'));
if nargin<1,folder=fullfile(root,'outputs','v68_validation','domain');end
bindings=sourceBindings(root);
[surface,truth]=buildRocketTestBody(struct('protrusionEnabled',false));
[mesh,domain]=buildRocketAxisymmetricDomain(surface,truth);
assert(domain.options.radialLayers==12&&mesh.shape(2)==12&&~domain.endcapsAdded&& ...
    ~domain.physicalSolution&&~domain.experimentalValidation&&~domain.eligibleForSurfaceCoverage&& ...
    domain.newProductionSolvedAreaM2==0&&~domain.nearWallResolutionValidated&&~domain.farfieldAdequacyValidated);

%% 一、鼻/侧壁/收缩/后体/真实底面全部保留；候选面映射不冒充场量投影
original=domain.originalWallPolylineXR;wall=domain.wallPolylineXR;
previous=0;originalWallNodeIndex=zeros(size(original,1),1);
for i=1:size(original,1)
    id=find(wall(:,1)==original(i,1)&wall(:,2)==original(i,2));
    assert(isscalar(id)&&id>previous,'RocketDomainTest:NativeWall','漏掉、重复或移动原生母线节点。');
    originalWallNodeIndex(i)=id;previous=id;
end
assert(isequal(original(1,:),[surface.V(surface.param.noseCapNodeID,1) 0])&&isequal(original(end,:),[0 0])&& ...
    isequal(domain.sourceVertices,surface.V)&&isequal(domain.sourceFaces,surface.F)&& ...
    isequal(domain.sourceFaceBandID,surface.param.faceBandID(:)));
regions=unique(domain.wallRegion);
assert(all(ismember({'sphericalNose';'cone';'cylinder';'contraction';'afterbody';'base'},regions)));
assert(isequal(sort(unique(domain.wallSourceFaceBandID)),sort(domain.originalWallFaceBandID))&& ...
    all(domain.kernelDistanceToOriginalWallM>0)&&all(diff(domain.wallPolarAngleRad)>0));
allCandidates=[];
for i=1:numel(domain.wallBoundaryFaceID)
    face=domain.wallBoundaryFaceID(i);ids=domain.wallSourceTriangleIDs{i};
    assert(isequal(ids,find(surface.param.faceBandID(:)==domain.wallSourceFaceBandID(i)))&& ...
        mesh.faces.owner(face)==sub2ind(mesh.shape,i,1)&&mesh.faces.boundarySide(face)==3);
    endpoints=domain.wallEndpointsXR(i,:,:);
    assert(isequal(reshape(endpoints(1,:,1),1,2),wall(i,:))&& ...
        isequal(reshape(endpoints(1,:,2),1,2),wall(i+1,:)));
    allCandidates=[allCandidates;ids]; %#ok<AGROW>
end
assert(isequal(unique(allCandidates),(1:size(surface.F,1))'),'RocketDomainTest:SourceFaces','真实墙面源三角面有遗漏。');
base=strcmp(domain.wallRegion,'base');assert(nnz(base)>=domain.options.baseRadialSegments&& ...
    all(wall([find(base);find(base)+1],1)==0)&& ...
    abs(domain.fullBaseAreaM2-pi*truth.baseRadius^2)<512*eps*pi*truth.baseRadius^2);
% 各源band的完整旋转面积与原折线直接积分相同；不乘任何热流修正系数。
originalArea=pi*(original(1:end-1,2)+original(2:end,2)).*hypot(diff(original(:,1)),diff(original(:,2)));
grouped=accumarray(domain.wallSourceOriginalSegment,2*pi*domain.wallAreaPerRadianM2,[numel(originalArea) 1]);
assert(max(abs(grouped-originalArea)./originalArea)<3e-12);

%% 二、独立Duffy三角积分完整A、V、两个r权一阶矩
moments=independentMoments(mesh.vertices,mesh.cells);
relativeArea=max(abs(moments(:,1)-mesh.planarArea)./moments(:,1));
relativeVolume=max(abs(moments(:,2)-mesh.volumePerRadian)./moments(:,2));
centroid=[moments(:,3)./moments(:,2),moments(:,4)./moments(:,2)];
centroidError=max(abs(centroid-mesh.volumeCentroid),[],'all');
assert(relativeArea<3e-11&&relativeVolume<3e-11&&centroidError<3e-11, ...
    'RocketDomainTest:IndependentMoments','独立真实r权积分与原几何不一致。');
outer=domain.outerBoundaryPolylineXR;expectedOuter=independentFrusta(outer);expectedBody=independentFrusta(original);
expectedFluid=expectedOuter-expectedBody;actualFluid=2*pi*sum(mesh.volumePerRadian);
assert(abs(actualFluid-expectedFluid)<=3e-12*(abs(actualFluid)+abs(expectedFluid))&& ...
    abs(2*pi*domain.bodyVolumePerRadianM3-expectedBody)<=3e-12*expectedBody);
assert(size(mesh.vertices,1)-numel(mesh.faces.owner)+mesh.cellCount==1, ...
    'RocketDomainTest:Topology','完整C域平面拓扑不是单一矩形同胚。');
inside=mesh.faces.neighbour>0;bf=mesh.boundaryFaces;
assert(numel(bf)==2*sum(mesh.shape)&&nnz(inside)==numel(mesh.faces.owner)-numel(bf)&& ...
    all(full(sum(abs(mesh.incidence),1))'==1+inside));
axisIDs=[domain.axisBoundaryFaceID{1};domain.axisBoundaryFaceID{2}];
assert(isequal(sort(axisIDs),find(ismember(mesh.faces.boundarySide,[1 2])))&& ...
    all(mesh.faces.areaPerRadian(axisIDs)==0)&&all(mesh.faces.centroid(axisIDs,2)==0)&& ...
    all(mesh.faces.normalDistance>0)&&all(mesh.volumeCentroid(:,2)>0)&& ...
    any(abs(domain.axisFaceOwnerXOffsetM)>1e-10*mesh.volumeCentroid(mesh.faces.owner(axisIDs),2)), ...
    'RocketDomainTest:Axis','轴几何必须真实r=0，且不能伪称符合旧正交模板。');
samples=[mesh.volumeCentroid;mesh.faces.centroid(bf,:)];linear=1+.2*samples(:,1)+.3*samples(:,2);
gradientError=max([abs(mesh.gradientX*linear-.2);abs(mesh.gradientR*linear-.3)]);
assert(gradientError<3e-10,'RocketDomainTest:Gradient','真实面样本下几何WLS不再现线性场。');

%% 三、独立几何参数对照保留同一原墙与映射，不把更大域当远场验收
[refined,refinedDomain]=buildRocketAxisymmetricDomain(surface,truth,struct('radialLayers',24));
[expanded,expandedDomain]=buildRocketAxisymmetricDomain(surface,truth,struct('outerRadius',1.5*domain.options.outerRadius));
[nearWall,nearWallDomain]=buildRocketAxisymmetricDomain(surface,truth,struct('firstRadialSpacing',.5*domain.options.firstRadialSpacing));
for value={refinedDomain,expandedDomain,nearWallDomain}
    item=value{1};assert(isequal(item.wallPolylineXR,wall)&& ...
        isequal(item.wallSourceFaceBandID,domain.wallSourceFaceBandID)&& ...
        isequal(item.wallSourceTriangleIDs,domain.wallSourceTriangleIDs)&& ...
        ~item.physicalSolution&&~item.farfieldAdequacyValidated);
end
assert(refined.cellCount==2*mesh.cellCount&&expanded.cellCount==mesh.cellCount&&nearWall.cellCount==mesh.cellCount&& ...
    expandedDomain.fluidVolumePerRadianM3>domain.fluidVolumePerRadianM3&& ...
    max(abs(nearWallDomain.wallNormalDistanceM./domain.wallNormalDistanceM-.5))<.01);

%% 四、非法源及域参数不得静默删突起、漏封闭区或改用等面积半径
refusals=struct('label',{},'identifier',{});
reject('unknown-domain-option',@()buildRocketAxisymmetricDomain(surface,truth,struct('guessExterior',true)));
reject('one-radial-layer',@()buildRocketAxisymmetricDomain(surface,truth,struct('radialLayers',1)));
reject('outside-star-kernel',@()buildRocketAxisymmetricDomain(surface,truth,struct('kernelX',.01)));
reject('outer-circle-inside-body',@()buildRocketAxisymmetricDomain(surface,truth,struct('outerRadius',.5*max(hypot(original(:,1)-domain.kernel(1),original(:,2))))));
reject('first-spacing-does-not-fit',@()buildRocketAxisymmetricDomain(surface,truth,struct('firstRadialSpacing',domain.options.outerRadius)));
[featured,featuredTruth]=buildRocketTestBody(struct('resolution',24));
reject('nonaxisymmetric-protrusion',@()buildRocketAxisymmetricDomain(featured,featuredTruth));
bad=surface;bad.V(surface.param.sectionNodeIDs{1}(1),2)=1.01*bad.V(surface.param.sectionNodeIDs{1}(1),2);
reject('changed-actual-ring-radius',@()buildRocketAxisymmetricDomain(bad,truth));
bad=surface;bad.param.noseCapNodeID=surface.param.sectionNodeIDs{end}(1);
reject('missing-original-nose-tip',@()buildRocketAxisymmetricDomain(bad,truth));
bad=surface;bad.param.faceBandID(1)=0;
reject('unknown-source-face-band',@()buildRocketAxisymmetricDomain(bad,truth));
bad=surface;a=find(surface.param.faceBandID==1,1);b=find(surface.param.faceBandID==2,1);
bad.param.faceBandID([a b])=bad.param.faceBandID([b a]);
reject('face-band-label-swapped-between-real-rings',@()buildRocketAxisymmetricDomain(bad,truth));
bad=surface;bad.F(end,:)=[];bad.param.faceBandID(end)=[];badTruth=truth;badTruth.faceBaseRegion(end)=[];
reject('open-original-surface',@()buildRocketAxisymmetricDomain(bad,badTruth));
assert(numel(refusals)==11&&isequaln(bindings,sourceBindings(root)),'RocketDomainTest:Source','测试期间源码改变。');
report=struct('status','passed-complete-native-rocket-exterior-geometry-only', ...
    'matlabEngineProcessID',feature('getpid'),'cellCount',mesh.cellCount,'nodeShape',size(mesh.X), ...
    'wallSegments',numel(domain.wallBoundaryFaceID),'originalWallSegments',numel(domain.originalWallFaceBandID), ...
    'originalWallNodesRetainedExactly',true,'originalWallNodeIndex',originalWallNodeIndex, ...
    'allOriginalSourceFacesBoundAsCandidates',true,'fieldProjectionPerformed',false, ...
    'independentRelativePlanarAreaError',relativeArea,'independentRelativeRWeightedVolumeError',relativeVolume, ...
    'independentRWeightedCentroidErrorM',centroidError,'geometricAffineGradientError',gradientError, ...
    'independentFluidVolumeM3',expectedFluid,'actualFluidVolumeM3',actualFluid, ...
    'fullRevolvedWallAreaM2',domain.fullRevolvedWallAreaM2,'fullBaseAreaM2',domain.fullBaseAreaM2, ...
    'sourceTriangulatedWallAreaM2',domain.sourceTriangulatedWallAreaM2,'analyticBodyWallAreaM2',domain.analyticBodyWallAreaM2, ...
    'minimumWallNormalDistanceM',min(domain.wallNormalDistanceM),'maximumWallNormalDistanceM',max(domain.wallNormalDistanceM), ...
    'minimumStarKernelDistanceM',min(domain.kernelDistanceToOriginalWallM), ...
    'maximumAxisOwnerXOffsetM',max(abs(domain.axisFaceOwnerXOffsetM)), ...
    'refinedCellCount',refined.cellCount,'expandedCellCount',expanded.cellCount,'refusals',refusals, ...
    'axisMirrorTemplateRequired',true,'flowResidualEvaluated',false,'nearWallResolutionValidated',false, ...
    'farfieldAdequacyValidated',false,'physicalSolution',false,'experimentalValidation',false, ...
    'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0,'sourceBindings',bindings);
if ~exist(folder,'dir'),mkdir(folder);end
save(fullfile(folder,'rocket_domain_evidence.mat'),'report','surface','truth','mesh','domain', ...
    'refined','refinedDomain','expanded','expandedDomain','nearWall','nearWallDomain');
fid=fopen(fullfile(folder,'rocket_domain_summary.json'),'w','n','UTF-8');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(report));clear cleanup;
fprintf('Complete rocket C-domain geometry: %d cells/%d wall segments, exact original seams and real base, independent r-weighted moments, 11 refusals; no flow/physical claim.\n', ...
    mesh.cellCount,numel(domain.wallBoundaryFaceID));
    function reject(label,action)
        try,action();catch exception
            assert(~isempty(exception.identifier),'RocketDomainTest:Refusal','拒绝必须有明确错误标识。');
            refusals(end+1)=struct('label',label,'identifier',exception.identifier);return; %#ok<AGROW>
        end
        error('RocketDomainTest:NoRefusal','非法合同%s未拒绝。',label);
    end
end

function moments=independentMoments(vertices,cells)
% 模块用途：在两真实三角形上直接积分1/r/xr/r²，与多边形矩实现独立。
g=[(1-sqrt(3/5))/2,.5,(1+sqrt(3/5))/2];w=[5/18,4/9,5/18];
moments=zeros(size(cells,1),4);
for triangle={[1 2 3],[1 3 4]}
    ids=triangle{1};A=vertices(cells(:,ids(1)),:);B=vertices(cells(:,ids(2)),:);C=vertices(cells(:,ids(3)),:);
    E=B-A;D=C-A;jac=abs(E(:,1).*D(:,2)-E(:,2).*D(:,1));
    for i=1:3
        for j=1:3
            point=A+g(i)*E+(1-g(i))*g(j)*D;weight=jac*(1-g(i))*w(i)*w(j);
            moments=moments+weight.*[ones(size(weight)),point(:,2),point(:,1).*point(:,2),point(:,2).^2];
        end
    end
end
end
function volume=independentFrusta(points)
% 模块用途：沿真实X区间逐圆台求三维体积，不从FV单元或构造器账取数。
volume=0;
for i=1:size(points,1)-1
    height=points(i,1)-points(i+1,1);a=points(i,2);b=points(i+1,2);
    volume=volume+pi*height*(a*a+a*b+b*b)/3;
end
end
function out=sourceBindings(root)
paths={'buildRocketTestBody.m','regionalFlow/axisymmetric/buildAxisymmetricFVGeometry.m', ...
    'regionalFlow/axisymmetric/buildRocketAxisymmetricDomain.m','tests/testRocketAxisymmetricDomain.m'};
out=struct('path',{},'sha256',{});
for i=1:numel(paths)
    path=fullfile(root,paths{i});fid=fopen(path,'rb');assert(fid>=0);cleanup=onCleanup(@()fclose(fid));
    bytes=fread(fid,Inf,'*uint8');clear cleanup;md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    digest=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));out(i)=struct('path',paths{i},'sha256',digest);
end
end
