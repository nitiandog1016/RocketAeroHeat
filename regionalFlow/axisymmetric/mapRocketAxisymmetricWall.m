function projection=mapRocketAxisymmetricWall(item)
%MAPROCKETAXISYMMETRICWALL 将细分真实墙段保守映射回原生箭体三角面。
% 模块一：重建原生完整外域，禁止把面元候选集合当成实际投影。
% 模块二：在母线参数t/周向角theta平面裁切原三角形，精确积分r(t)L dt dtheta。
% 模块三：鼻尖和底心的单一扇面对应完整参数矩形，正确处理极点退化。
% 模块四：输出真实曲面面积矩阵/加权矩阵；显示弦面面积分列，不改热流原值。
% 曲面指折线母线旋转面；只保证标量积分，不保证由弦面法向积分的压力合力/力矩。
% 此模块仅建立几何投影；数值状态、物理适用性及生产资格由调用者另行审计。
required={'mesh','domain','surfaceMesh','geometryTruth','options'};
assert(isstruct(item)&&isscalar(item)&&all(isfield(item,required)), ...
    'RocketMap:Inputs','需要原完整工况中的几何、真实域及源表面。');
geometry=struct('protrusionEnabled',false);
if isfield(item.options,'geometry'),geometry=item.options.geometry;end
if ~isfield(geometry,'protrusionEnabled'),geometry.protrusionEnabled=false;end
[nativeSurface,nativeTruth]=buildRocketTestBody(geometry);
assert(isequaln(nativeSurface,item.surfaceMesh)&&isequaln(nativeTruth,item.geometryTruth), ...
    'RocketMap:NativeSource','只能映射同一输入重新生成的原生V/F，不接受任意环带重排/替换。');
if isfield(item.domain,'meshFactory')
    assert(isfield(item.domain,'collar')&&isfield(item.domain.collar,'options'), ...
        'RocketMap:Binding','未知网格工厂或缺少真实近壁层配置。');
    if strcmp(item.domain.meshFactory,'buildRocketAxisymmetricCollarDomain')
        [mesh,domain]=buildRocketAxisymmetricCollarDomain(item.surfaceMesh,item.geometryTruth,item.domain.options,item.domain.collar.options);
    elseif strcmp(item.domain.meshFactory,'refineRocketAxisymmetricRadial')
        % 逐值重建完整父几何再真实二分，不能只因源墙相同就跳过完整网格绑定。
        assert(isfield(item.domain.collar,'parentRequestedOptions')&& ...
            item.domain.collar.options.radialRefinementFactor==2&&mod(item.domain.options.radialLayers,2)==0, ...
            'RocketMap:Binding','二分域缺少明确父网格输入。');
        parentOptions=item.domain.options;parentOptions.radialLayers=parentOptions.radialLayers/2;
        parentOptions.firstRadialSpacing=parentOptions.firstRadialSpacing*2;
        [parentMesh,parentDomain]=buildRocketAxisymmetricCollarDomain(item.surfaceMesh,item.geometryTruth, ...
            parentOptions,item.domain.collar.parentRequestedOptions);
        [mesh,domain]=refineRocketAxisymmetricRadial(parentMesh,parentDomain);
    else,error('RocketMap:Binding','未知网格工厂不得省略完整几何重建。');end
else
    [mesh,domain]=buildRocketAxisymmetricDomain(item.surfaceMesh,item.geometryTruth,item.domain.options);
end
assert(isequaln(mesh,item.mesh)&&isequaln(domain,item.domain), ...
    'RocketMap:Binding','完整域/原生表面不能逐值重建。');
surface=item.surfaceMesh;V=surface.V;F=surface.F;meta=item.geometryTruth;
original=domain.originalWallPolylineXR;wall=domain.wallPolylineXR;
parent=domain.wallSourceOriginalSegment;bands=domain.originalWallFaceBandID;
nWall=numel(parent);nFace=size(F,1);
faceParent=zeros(nFace,1);parameterAreas=zeros(nFace,1);poleSector=false(nFace,1);
rowIDs=cell(nFace,1);columnIDs=rowIDs;overlapAreas=rowIDs;
scale=max(1,max(abs(V),[],'all'));tol=2048*eps(scale);

%% 一、逐原始三角形建立真实参数域；极点扇区不能少算一半面积
for f=1:nFace
    p=find(bands==surface.param.faceBandID(f));
    assert(isscalar(p),'RocketMap:SourceBand','源三角形band不能唯一对应原母线。');
    faceParent(f)=p;A=original(p,:);B=original(p+1,:);E=B-A;L=norm(E);
    ids=F(f,:);points=[V(ids,1),meta.nodeBaselineRadius(ids)];
    isA=vecnorm(points-A,2,2)<=tol;isB=vecnorm(points-B,2,2)<=tol;
    assert(all(xor(isA,isB))&&any(isA)&&any(isB), ...
        'RocketMap:SourceVertices','源三角形节点必须位于该原母线的两个端环/极点。');
    t=double(isB);nonpole=meta.nodeBaselineRadius(ids)>0;
    angles=atan2(V(ids(nonpole),3),V(ids(nonpole),2));
    relative=atan2(sin(angles-angles(1)),cos(angles-angles(1)));
    angularSpan=max(relative)-min(relative);localTheta=relative-min(relative);
    assert(angularSpan>0&&angularSpan<pi&& ...
        all(min(abs(localTheta),abs(localTheta-angularSpan))<=4096*eps), ...
        'RocketMap:SourceAngles','源三角形须横跨同一相邻周向角扇区。');
    if nnz(nonpole)==2
        assert(nnz(~nonpole)==1&&(A(2)==0||B(2)==0), ...
            'RocketMap:Pole','退化参数域只用于真实鼻尖或底心。');
        polygon=[0 0;1 0;1 angularSpan;0 angularSpan];poleSector(f)=true;
    else
        assert(nnz(nonpole)==3&&A(2)>0&&B(2)>0, ...
            'RocketMap:Pole','普通环带三角形不能包含未说明的极点。');
        polygon=[t,localTheta];
    end
    parameterAreas(f)=integrateRadius(polygon,A(2),B(2),L);
    candidates=find(parent==p);weights=zeros(numel(candidates),1);
    for k=1:numel(candidates)
        j=candidates(k);
        lo=dot(wall(j,:)-A,E)/L^2;hi=dot(wall(j+1,:)-A,E)/L^2;
        assert(lo>=-tol/L&&hi<=1+tol/L&&hi>lo, ...
            'RocketMap:WallInterval','墙段不能超出或反向跨越其真实原母线。');
        % 只将已证明属于原端点的几何舍入还原为参数0/1，不接触任何流场状态。
        if norm(wall(j,:)-A)<=tol,lo=0;end
        if norm(wall(j+1,:)-B)<=tol,hi=1;end
        clipped=clipAt(clipAt(polygon,lo,true),hi,false);
        weights(k)=integrateRadius(clipped,A(2),B(2),L);
    end
    keep=weights>0;rowIDs{f}=repmat(f,nnz(keep),1);
    columnIDs{f}=candidates(keep);overlapAreas{f}=weights(keep);
end
areaMatrix=sparse(vertcat(rowIDs{:}),vertcat(columnIDs{:}),vertcat(overlapAreas{:}),nFace,nWall);
sourceArea=full(sum(areaMatrix,2));wallArea=2*pi*domain.wallAreaPerRadianM2;
mappedWallArea=full(sum(areaMatrix,1))';
rowRelative=max(abs(sourceArea-parameterAreas)./parameterAreas);
columnRelative=max(abs(mappedWallArea-wallArea)./wallArea);
assert(all(sourceArea>0)&&all(parameterAreas>0)&&all(nonzeros(areaMatrix)>0)&& ...
    rowRelative<=2e-11&&columnRelative<=2e-11, ...
    'RocketMap:AreaLedger','真实源参数面或墙段曲面面积没有完整一一覆盖。');
weights=spdiags(1./sourceArea,0,nFace,nFace)*areaMatrix;
assert(max(abs(full(sum(weights,2))-1))<=4096*eps, ...
    'RocketMap:Weights','逐面几何平均权重没有保持常值。');
crosses=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
facetArea=.5*vecnorm(crosses,2,2);

%% 二、完整投影账；p/q等场量可按sourceValues=weights*wallValues映射
projection=struct('kind','native-rocket-source-face-curved-area-conservative-projection-v1', ...
    'sourceVertices',V,'sourceFaces',F,'sourceFaceID',(1:nFace)', ...
    'sourceFaceBandID',surface.param.faceBandID(:),'sourceOriginalSegment',faceParent, ...
    'sourceFaceIsPoleSector',poleSector,'sourceCurvedAreaM2',sourceArea, ...
    'sourceDisplayFacetAreaM2',facetArea,'wallBoundaryFaceID',domain.wallBoundaryFaceID, ...
    'wallAreaM2',wallArea,'overlapAreaM2',areaMatrix,'weights',weights, ...
    'rowAreaMaxRelativeError',rowRelative,'columnAreaMaxRelativeError',columnRelative, ...
    'completeCurvedAreaM2',sum(sourceArea),'completeFacetAreaM2',sum(facetArea), ...
    'projectionMethod','clip source (t,theta) parameter domain at actual wall intervals; integrate r(t)*L exactly', ...
    'poleMethod','one native pole triangle corresponds to full angular sector parameter rectangle', ...
    'valueMethod','sourceValue=weights*wallValue; true curved-area integral conserved; signed q unchanged before averaging', ...
    'displayFacetIntegralIsCurvedIntegral',false,'fieldValuesProvided',false, ...
    'surfaceDefinition','exact integration of piecewise-linear generating meridian; not analytic spherical-nose surface', ...
    'pressureForceOrMomentConservationEstablished',false, ...
    'originalSourceFacesPreserved',true,'sourceIDsAreProjectionCandidatesOnly',false, ...
    'physicalSolution',false,'experimentalValidation',false,'eligibleForSurfaceCoverage',false, ...
    'newProductionSolvedAreaM2',0);
end

function out=clipAt(polygon,cut,keepAbove)
% 模块用途：按母线参数直线裁切凸参数三角形/矩形，保留真实周向范围。
if isempty(polygon),out=zeros(0,2);return;end
out=zeros(0,2);previous=polygon(end,:);previousInside=isInside(previous(1));
for k=1:size(polygon,1)
    current=polygon(k,:);currentInside=isInside(current(1));
    if currentInside~=previousInside
        alpha=(cut-previous(1))/(current(1)-previous(1));
        point=previous+alpha*(current-previous);point(1)=cut;out(end+1,:)=point; %#ok<AGROW>
    end
    if currentInside,out(end+1,:)=current;end %#ok<AGROW>
    previous=current;previousInside=currentInside;
end
    function value=isInside(t)
        if keepAbove,value=t>=cut;else,value=t<=cut;end
    end
end

function value=integrateRadius(polygon,r0,r1,L)
% 模块用途：Green精确面积及一阶t矩，积分线性半径；不采用面心近似权重。
if size(polygon,1)<3,value=0;return;end
next=polygon([2:end 1],:);crosses=polygon(:,1).*next(:,2)-next(:,1).*polygon(:,2);
area=sum(crosses)/2;momentT=sum((polygon(:,1)+next(:,1)).*crosses)/6;
if area<0,area=-area;momentT=-momentT;end
value=L*(r0*area+(r1-r0)*momentT);
assert(isfinite(value)&&value>=0,'RocketMap:ParameterArea','参数曲面面积非有限或负值。');
end
