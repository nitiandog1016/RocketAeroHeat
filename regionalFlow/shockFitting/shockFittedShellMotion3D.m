function motion=shockFittedShellMotion3D(shell,nodes,faceNormalSpeed)
%SHOCKFITTEDSHELLMOTION3D 将逐面物理激波速度投影到共享节点，固定真实壁面。
% 模块一：构造真实面重心速度C*c，每个面同时耦合三个共享节点。
% 模块二：以sqrt(面积)加权的稀疏LSQR求全局径向最小二乘，不以逐节点平均冒称全局投影。
% 模块三：同一节点只有一份速度，内层线性随动，原壁面严格保持0。
% 模块四：记录原方程驻定性与逐面残差；未表示的高频面速度不能冒称已经相容。
mesh=shell.mesh;N=numel(shell.originalVertexIDs);outer=mesh.boundaryFaces(shell.outerBoundaryMask);
assert(isequal(size(faceNormalSpeed),[numel(outer) 1])&&all(isfinite(faceNormalSpeed)), ...
    'ShockFit:ShockSpeed','每个实际外壳面需要一个物理法向速度。');
g=shockFittedALEGeometry(mesh,nodes,nodes,1);normal=g.normal(outer,:);area=g.averageAreaM2(outer);
body=nodes(1:N,:);outside=nodes(shell.layers*N+(1:N),:);rays=outside-body;
assert(all(vecnorm(rays,2,2)>0),'ShockFit:ShockGap','外壳节点穿过或落在真实壁面上。');
directions=rays./vecnorm(rays,2,2);faces=mesh.faces(outer,:)-shell.layers*N;
coefficient=zeros(numel(outer),3);
for j=1:3
    dotn=sum(directions(faces(:,j),:).*normal,2);
    assert(all(dotn>0),'ShockFit:ShockRays','壁面射线与外壳外法向不一致。');
    coefficient(:,j)=dotn/3;
end
C=sparse(repmat((1:numel(outer))',3,1),faces(:),coefficient(:),numel(outer),N);
assert(all(sum(abs(C),1)>0),'ShockFit:ShockMotion','存在未被外壳面约束的节点。');
A=spdiags(sqrt(area),0,numel(area),numel(area))*C;rhs=sqrt(area).*faceNormalSpeed;
if norm(rhs)==0,radialSpeed=zeros(N,1);flag=0;iterations=0;
else,[radialSpeed,flag,~,iterations]=lsqr(A,rhs,1e-12,4*N);end
residual=A*radialSpeed-rhs;
stationarity=norm(A'*residual)/max(realmin,norm(A,'fro')*norm(rhs));
assert(stationarity<2e-10&&all(isfinite(radialSpeed)),'ShockFit:ShockMotionSolve', ...
    '共享节点全局最小二乘未达到原法方程驻定条件，不提交运动。');
outerVelocity=radialSpeed.*directions;
velocity=zeros(size(nodes));
for layer=1:shell.layers,velocity(layer*N+(1:N),:)=layer/shell.layers*outerVelocity;end
actual=sum((outerVelocity(faces(:,1),:)+outerVelocity(faces(:,2),:)+outerVelocity(faces(:,3),:))/3.*normal,2);
motion=struct('nodeVelocityMS',velocity,'faceNormalSpeedMS',actual,'requestedFaceNormalSpeedMS',faceNormalSpeed, ...
    'projectionResidualMS',actual-faceNormalSpeed,'radialNodeSpeedMS',radialSpeed, ...
    'normalEquationStationarityRelative',stationarity,'leastSquaresFlag',flag,'leastSquaresIterations',iterations, ...
    'weightedProjectionRmsMS',norm(residual)/sqrt(sum(area)), ...
    'method','真实三节点面重心法向速度的全局面积加权稀疏LSQR/固定壁面/射线内层随动', ...
    'perFaceIndependentShapeFit',false,'bodyNodesFixed',all(velocity(1:N,:)==0,'all'), ...
    'physicalQualification',false);
end
