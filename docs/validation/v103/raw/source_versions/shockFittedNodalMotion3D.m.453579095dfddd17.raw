function motion=shockFittedNodalMotion3D(shell,nodes,Q,upstream,gas,maximumArrayBytes)
%SHOCKFITTEDNODALMOTION3D 在共享激波节点上求单一相容速度，内层随动且壁面固定。
% 模块一：真实外壳三角面面积向量给出节点法向，不用解析球面法向替代当前几何。
% 模块二：每节点重构下游状态及J+相容速度，再沿固定壁面射线移动唯一节点。
% 模块三：仍记录对原单元常值逐面Riemann目标的误差，不能把节点恒等式当逐面合格。
if nargin<6,maximumArrayBytes=128*1024^2;end
mesh=shell.mesh;N=numel(shell.originalVertexIDs);outer=mesh.boundaryFaces(shell.outerBoundaryMask);
nodeIDs=shell.layers*N+(1:N)';faces=mesh.faces(outer,:)-shell.layers*N;
g=shockFittedALEGeometry(mesh,nodes,nodes,1);vector=g.areaVector0(outer,:);nodalVector=zeros(N,3);
for j=1:3
    for k=1:3,nodalVector(:,k)=nodalVector(:,k)+accumarray(faces(:,j),vector(:,k),[N 1]);end
end
assert(all(vecnorm(nodalVector,2,2)>0),'ShockFit:NodalNormal','外壳节点法向退化。');
normal=nodalVector./vecnorm(nodalVector,2,2);
trace=reconstructShockFittedNodalTrace3D(Q,mesh,nodes,gas,nodeIDs,maximumArrayBytes);
characteristic=shockFittedCharacteristicSpeed3D(trace.Q,normal,upstream,gas);
ray=nodes(nodeIDs,:)-nodes(1:N,:);ray=ray./vecnorm(ray,2,2);cosine=sum(ray.*normal,2);
assert(all(cosine>0),'ShockFit:NodalRay','节点射线与实际外壳法向反向。');
radial=characteristic.normalSpeedMS./cosine;outsideVelocity=radial.*ray;
velocity=zeros(size(nodes));
for layer=1:shell.layers,velocity(layer*N+(1:N),:)=layer/shell.layers*outsideVelocity;end
actual=sum((outsideVelocity(faces(:,1),:)+outsideVelocity(faces(:,2),:)+outsideVelocity(faces(:,3),:))/3.*g.normal(outer,:),2);
reference=shockFittedBoundarySpeed3D(Q(mesh.owner(outer),:),g.normal(outer,:),upstream,gas);
residual=actual-reference.normalSpeedMS;area=g.averageAreaM2(outer);
facePrimitive=(trace.primitiveRP(faces(:,1),:)+trace.primitiveRP(faces(:,2),:)+trace.primitiveRP(faces(:,3),:))/3;
motion=struct('nodeVelocityMS',velocity,'faceNormalSpeedMS',actual, ...
    'requestedFaceNormalSpeedMS',reference.normalSpeedMS,'projectionResidualMS',residual, ...
    'radialNodeSpeedMS',radial,'nodalNormal',normal,'nodalNormalSpeedMS',characteristic.normalSpeedMS, ...
    'normalEquationStationarityRelative',NaN,'leastSquaresFlag',NaN,'leastSquaresIterations',NaN, ...
    'weightedProjectionRmsMS',sqrt(sum(area.*residual.^2)/sum(area)), ...
    'nodalTrace',trace,'nodalCharacteristic',characteristic,'referenceCellFaceShock',reference, ...
    'reconstructedFacePrimitiveRP',facePrimitive, ...
    'method','共享节点流体侧仿射重构/J+与RH相容/实际面积法向/固定壁面射线及内层随动', ...
    'referenceFaceResidualMeaning','对原单元常值逐面精确Riemann目标的误差；本方法不求该逐面最小二乘', ...
    'perFaceIndependentShapeFit',false,'bodyNodesFixed',all(velocity(1:N,:)==0,'all'), ...
    'physicalQualification',false);
end
