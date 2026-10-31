function cap=buildShockFittedSphericalCap(radius,thetaMaxDeg,rings,azimuths,gap,layers)
%BUILDSHOCKFITTEDSPHERICALCAP 构造完整周向球鼻局部域；不使用复制的轴对称流场。
% 模块一：鼻尖沿+X、固定壁面为真实三角面，共享极点及周向接缝。
% 模块二：初始外壳径向偏移；后续由下游相容速度推进，不把偏移壳称作已求激波。
% 模块三：保留解析球半径和面元弦高误差，粗网格不能冒称几何精度已收敛。
validateattributes(radius,{'double'},{'real','finite','positive','scalar'});
validateattributes(thetaMaxDeg,{'double'},{'real','finite','positive','scalar','<',90});
validateattributes(rings,{'double'},{'real','finite','integer','scalar','>=',2});
validateattributes(azimuths,{'double'},{'real','finite','integer','scalar','>=',6});
V=zeros(1+rings*azimuths,3);V(1,:)=[radius 0 0];theta=linspace(0,deg2rad(thetaMaxDeg),rings+1);
for ring=1:rings
    phi=(0:azimuths-1)'*2*pi/azimuths;t=theta(ring+1);
    V(1+(ring-1)*azimuths+(1:azimuths),:)=radius*[repmat(cos(t),azimuths,1),sin(t)*cos(phi),sin(t)*sin(phi)];
end
F=zeros(azimuths+(rings-1)*2*azimuths,3);q=0;
for j=1:azimuths,q=q+1;F(q,:)=[1 1+j 1+mod(j,azimuths)+1];end
for ring=1:rings-1,for j=1:azimuths
    next=mod(j,azimuths)+1;a=1+(ring-1)*azimuths+j;b=1+(ring-1)*azimuths+next;
    c=1+ring*azimuths+j;d=1+ring*azimuths+next;
    F(q+(1:2),:)=[a c d;a d b];q=q+2;
end,end
normal=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
centers=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;negative=sum(normal.*centers,2)<0;
F(negative,[2 3])=F(negative,[3 2]);shell=buildShockFittedSurfaceShell(V,F,(1:size(F,1))',gap,layers);
N=numel(shell.originalVertexIDs);body=V(shell.originalVertexIDs,:);shell.vertexNormals=body/radius;
for layer=0:layers,shell.nodes(layer*N+(1:N),:)=body*(1+layer*gap/(layers*radius));end
% 解析径向射线保证侧边界是穿过球心的多边形锥面，避免外壳法向平均改变出口几何。
g=shockFittedALEGeometry(shell.mesh,shell.nodes,shell.nodes,1); %#ok<NASGU>
cap=struct('V',V,'F',F,'shell',shell,'radiusM',radius,'thetaMaxDeg',thetaMaxDeg,'rings',rings, ...
    'azimuths',azimuths,'maximumFacetCentroidSagM',max(radius-vecnorm(centers,2,2)), ...
    'originalBodyVerticesExact',true,'analyticSphereGeometryExact',false,'productionQualified',false);
end
