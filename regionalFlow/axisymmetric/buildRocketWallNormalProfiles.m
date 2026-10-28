function profiles=buildRocketWallNormalProfiles(item,result)
%BUILDROCKETWALLNORMALPROFILES 从实际五量控制体读取每个真实墙段的法向剖面。
% 模块一：取真实r权体质心、墙法向投影距离及同条带T/u/p/rho/SA黏度代理。
% 模块二：给定Tw/墙速单列，内部温度不是Te、Tr或求解壁温；不拟合边界层厚度。
% 模块三：附各径向内面实际Euler能量通量相对算术物理通量之差与实际黏性功/热总能通量。
% 该差仅作离散通量诊断，不冒称独立误差估计、纯数值热流或物理模型验证。
Q=result.Q;mesh=item.mesh;g=item.gas;P=axisymmetricFVState(Q(:,1:4),g);
assert(isequal(P,result.finalInfo.primitive(:,1:6)),'RocketWallProfile:State','必须来自同一完整五量端点。');
ni=mesh.shape(1);nj=mesh.shape(2);wall=item.domain.wallPolylineXR;
center=mesh.faces.centroid(item.domain.wallBoundaryFaceID,:);edge=diff(wall);tangent=edge./vecnorm(edge,2,2);
normal=[tangent(:,2),-tangent(:,1)];profiles=struct([]);f=mesh.faces;
radial=find(f.neighbour>0&f.neighbour-f.owner==ni);
assert(numel(radial)==ni*(nj-1),'RocketWallProfile:Topology','完整径向共享面数量不符。');
for i=1:ni
    ids=(i:ni:mesh.cellCount)';distance=sum((mesh.volumeCentroid(ids,:)-center(i,:)).*normal(i,:),2);
    assert(all(distance>0)&&all(diff(distance)>0),'RocketWallProfile:Distance','同条带法向质心距离须严格递增。');
    faces=radial(mod(f.owner(radial)-1,ni)+1==i);[~,order]=sort(f.owner(faces));faces=faces(order);
    left=f.owner(faces);right=f.neighbour(faces);n=f.normal(faces,:);
    ul=sum(P(left,2:3).*n,2);ur=sum(P(right,2:3).*n,2);
    centralEnergy=.5*((Q(left,4)+P(left,5)).*ul+(Q(right,4)+P(right,5)).*ur);
    actualEuler=result.finalInfo.faceConvectiveFluxPerRadian(faces,4)./f.areaPerRadian(faces);
    actualViscous=result.finalInfo.faceViscousFluxPerRadian(faces,4)./f.areaPerRadian(faces);
    row=struct('wallSegment',i,'region',item.domain.wallRegion{i},'wallCentroidXR',center(i,:), ...
        'cellIDs',ids,'cellCentroidXR',mesh.volumeCentroid(ids,:),'distanceM',distance, ...
        'temperatureK',P(ids,4),'densityKgM3',P(ids,1),'pressurePa',P(ids,5), ...
        'speedMps',vecnorm(P(ids,2:3),2,2),'tangentialVelocityMps',sum(P(ids,2:3).*tangent(i,:),2), ...
        'muTOverMu',result.finalInfo.sa.cellMuT(ids)./result.finalInfo.transport.cellMu(ids), ...
        'wallTemperatureInputK',item.inputContract.wallTemperatureInputK,'wallVelocityInputMps',[0 0], ...
        'wallBoundaryValueKind','prescribed','edgeTemperatureSolved',false,'recoveryTemperatureSolved',false,'wallTemperatureSolved',false, ...
        'radialFaceID',faces,'radialFaceNormalDistanceM',sum((f.centroid(faces,:)-center(i,:)).*normal(i,:),2), ...
        'EulerEnergyFluxMinusArithmeticPhysicalWm2',actualEuler-centralEnergy, ...
        'viscousTotalEnergyFluxWm2',actualViscous,'physicalSolution',false);
    if isempty(profiles),profiles=row;else,profiles(end+1)=row;end %#ok<AGROW>
end
end
