function geometry=shockFittedALEGeometry(mesh,nodes0,nodes1,dt)
%SHOCKFITTEDALEGEOMETRY 线性节点运动下共享三角面的精确扫掠体积及离散GCL。
% 模块一：面积向量是时间二次式，用Simpson精确积分；不把端点面积当平均面积。
% 模块二：同一共享面只算一次扫掠量，owner和neighbour使用相反符号。
% 模块三：回代每单元V1-V0=Σ扫掠体积；几何恒等式通过不代表物理流场已解。
assert(isa(nodes0,'double')&&isa(nodes1,'double')&&isreal(nodes0)&&isreal(nodes1)&& ...
    isequal(size(nodes0),[mesh.nodeCount 3])&&isequal(size(nodes1),size(nodes0))&& ...
    all(isfinite(nodes0),'all')&&all(isfinite(nodes1),'all'),'ShockFit:MeshMotion','节点尺寸或数值不符。');
validateattributes(dt,{'double'},{'real','finite','positive','scalar'});
v0=volumes(nodes0);v1=volumes(nodes1);vm=volumes((nodes0+nodes1)/2);
assert(all(v0>0&v1>0&vm>0),'ShockFit:CellVolume','运动路径的端点或中点单元体积非正。');
% 体积是三次多项式；额外检查区间内导数根，拒绝只在中间翻转又恢复的单元。
va=volumes(nodes0+(nodes1-nodes0)/3);vb=volumes(nodes0+2*(nodes1-nodes0)/3);
coeff=[ones(4,1),[0;1/3;2/3;1],[0;1/3;2/3;1].^2,[0;1/3;2/3;1].^3]\[v0';va';vb';v1'];
for j=1:mesh.cellCount
    roots_=roots([3*coeff(4,j),2*coeff(3,j),coeff(2,j)]);
    for t=roots_'
        if isreal(t)&&t>0&&t<1
            assert(coeff(:,j)'*[1;t;t^2;t^3]>0,'ShockFit:CellVolume','单元在运动中途翻转。');
        end
    end
end
A0=area(nodes0);Am=area((nodes0+nodes1)/2);A1=area(nodes1);Abar=(A0+4*Am+A1)/6;
centroidVelocity=((nodes1(mesh.faces(:,1),:)-nodes0(mesh.faces(:,1),:))+ ...
    (nodes1(mesh.faces(:,2),:)-nodes0(mesh.faces(:,2),:))+ ...
    (nodes1(mesh.faces(:,3),:)-nodes0(mesh.faces(:,3),:)))/(3*dt);
swept=dt*sum(Abar.*centroidVelocity,2);areaMagnitude=vecnorm(Abar,2,2);
assert(all(areaMagnitude>0),'ShockFit:FaceArea','平均面积向量退化。');
normal=Abar./areaMagnitude;normalSpeed=swept./(dt*areaMagnitude);
closure=mesh.incidence*Abar;areaSum=abs(mesh.incidence)*areaMagnitude;
gcl=(v1-v0)-mesh.incidence*swept;
relative=abs(gcl)./max(v0,v1);areaRelative=vecnorm(closure,2,2)./areaSum;
assert(max(relative)<2e-11&&max(areaRelative)<2e-12,'ShockFit:GCL','离散几何守恒或面积闭合没有通过。');
geometry=struct('volume0',v0,'volume1',v1,'areaVector0',A0,'areaVector1',A1, ...
    'averageAreaVector',Abar,'averageAreaM2',areaMagnitude,'normal',normal, ...
    'sweptVolumeM3',swept,'gridNormalSpeedMS',normalSpeed,'gclResidualM3',gcl, ...
    'gclRelative',relative,'areaClosureRelative',areaRelative, ...
    'faceCentroidMidM',((nodes0(mesh.faces(:,1),:)+nodes1(mesh.faces(:,1),:))+ ...
    (nodes0(mesh.faces(:,2),:)+nodes1(mesh.faces(:,2),:))+ ...
    (nodes0(mesh.faces(:,3),:)+nodes1(mesh.faces(:,3),:)))/6, ...
    'method','共享三角面线性节点运动/二次面积Simpson积分/单元扫掠体积离散GCL');

    function v=volumes(X)
        t=mesh.tetrahedra;a=X(t(:,1),:);b=X(t(:,2),:);c=X(t(:,3),:);d=X(t(:,4),:);
        v=sum((b-a).*cross(c-a,d-a,2),2)/6;
    end
    function A=area(X)
        f=mesh.faces;A=.5*cross(X(f(:,2),:)-X(f(:,1),:),X(f(:,3),:)-X(f(:,1),:),2);
    end
end
