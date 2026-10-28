function [Q,audit]=mapRocketAxisymmetricInitialGuess(parent,parentQ,item)
%MAPROCKETAXISYMMETRICINITIALGUESS 同一完整母线上的内部数值初猜插值。
% 模块一：逐流向条带按体质心至墙心的实际距离作五量凸线性插值。
% 模块二：父域内端点采用最近父格；父外域之外保留新工况的自由流初猜。
% 模块三：保存两个供体、权重和域外标志；不宣称保守重映射或已解流场。
assert(isequal(parent.domain.wallPolylineXR,item.domain.wallPolylineXR)&& ...
    isequaln(parent.gas,item.gas)&&isequaln(parent.boundary,item.boundary)&& ...
    isequal(parent.residualScale,item.residualScale)&& ...
    isequal(size(parentQ),[parent.mesh.cellCount 5]), ...
    'RocketSensitivity:Mapping','初猜须来自同一完整母线、气体、边界和五量尺度。');
ni=parent.mesh.shape(1);assert(item.mesh.shape(1)==ni,'RocketSensitivity:Mapping','不能跨不同墙分段静默插值。');
polygon=[parent.domain.wallPolylineXR;flipud(parent.domain.outerBoundaryPolylineXR)];
p=item.mesh.volumeCentroid;inside=inpolygon(p(:,1),p(:,2),polygon(:,1),polygon(:,2));
Q=item.Q;donorIDs=zeros(item.mesh.cellCount,2);weights=zeros(item.mesh.cellCount,2);
endpoint=false(item.mesh.cellCount,1);wall=parent.mesh.faces.centroid(parent.domain.wallBoundaryFaceID,:);
for i=1:ni
    old=(i:ni:parent.mesh.cellCount)';new=(i:ni:item.mesh.cellCount)';
    distance=vecnorm(parent.mesh.volumeCentroid(old,:)-wall(i,:),2,2);
    query=vecnorm(p(new,:)-wall(i,:),2,2);
    assert(all(diff(distance)>0),'RocketSensitivity:Mapping','父条带质心距离须严格递增。');
    for k=1:numel(new)
        id=new(k);if ~inside(id),continue;end
        x=query(k);
        if x<=distance(1),a=1;b=1;t=0;endpoint(id)=true;
        elseif x>=distance(end),a=numel(old);b=a;t=0;endpoint(id)=true;
        else,a=find(distance<=x,1,'last');b=a+1;t=(x-distance(a))/(distance(b)-distance(a));end
        donorIDs(id,:)=[old(a),old(b)];weights(id,:)=[1-t,t];
        Q(id,:)=(1-t)*parentQ(old(a),:)+t*parentQ(old(b),:);
    end
end
axisymmetricFVState(Q(:,1:4),item.gas);assert(all(Q(:,5)>=0),'RocketSensitivity:Mapping','不能裁负SA来制造合法初猜。');
audit=struct('kind','same-wall-strip-physical-distance-convex-linear-five-quantity-initial-guess-v1', ...
    'donorCellIDs',donorIDs,'donorWeights',weights,'insideParentDomain',inside, ...
    'endpointConstantExtension',endpoint,'outsideParentDomainUsesUniformFreestream',~inside, ...
    'parentCellCount',parent.mesh.cellCount,'newCellCount',item.mesh.cellCount, ...
    'parentConservedIntegralPerRadian',sum(parentQ.*parent.mesh.volumePerRadian,1), ...
    'newInitialConservedIntegralPerRadian',sum(Q.*item.mesh.volumePerRadian,1), ...
    'conservativeRemappingClaimed',false,'initialGuessIsBoundaryCondition',false, ...
    'externalCFDImported',false,'stateClipped',false,'physicalSolution',false,'newProductionSolvedAreaM2',0);
end
