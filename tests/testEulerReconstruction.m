function testEulerReconstruction()
%TESTEULERRECONSTRUCTION 独立解析积分、仿射场及间断正状态检验。
[X,R]=ndgrid([-2,-1.2,0],[0,.1,.4]);g=struct('X',X,'R',R,'nx',2,'nr',2);
geo=axisymmetricCellGeometry(X,R);dx=diff(X(:,1));lo=R(1,1:end-1);hi=R(1,2:end);
exactVolume=dx*.5.*(hi.^2-lo.^2);
exactR=repmat(2/3*(hi.^3-lo.^3)./(hi.^2-lo.^2),2,1);
assert(max(abs(geo.volume-exactVolume),[],'all')<1e-14);
assert(max(abs(geo.centerR-exactR),[],'all')<1e-14);
assert(max(abs(geo.centerX-.5*(X(1:end-1,1:end-1)+X(2:end,1:end-1))),[],'all')<1e-14);
assert(max(abs(geo.areaCenterR-repmat(.5*(lo+hi),2,1)),[],'all')<1e-14);
shifted=axisymmetricCellGeometry(X+1e6,R);
assert(max(abs(shifted.volume-geo.volume),[],'all')<1e-10);
assert(max(abs(shifted.centerX-1e6-geo.centerX),[],'all')<1e-9);
rec=prepareEulerReconstruction(g);q=affine(geo.centerX,geo.centerR);
for limited=false
    [faces,source]=reconstructEulerPrimitive(q,rec,limited);
    xx={geo.faceAX(1:end-1,:),geo.faceAX(2:end,:),geo.faceBX(:,1:end-1),geo.faceBX(:,2:end)};
    rr={geo.faceAR(1:end-1,:),geo.faceAR(2:end,:),geo.faceBR(:,1:end-1),geo.faceBR(:,2:end)};
    for k=1:4,assert(max(abs(faces{k}-affine(xx{k},rr{k})),[],'all')<1e-12);end
    expected=affine(geo.areaCenterX,geo.areaCenterR);
    assert(max(abs(source-expected(:,:,4)),[],'all')<1e-12);
end
% 非正交网格上不限幅线性重构也应准确；不把面中心当作节点中点。
[X,R]=ndgrid(linspace(0,1,8),linspace(.1,1,7));X=X+.2*R+.1*R.^2;
g=struct('X',X,'R',R,'nx',7,'nr',6);rec=prepareEulerReconstruction(g);geo=rec.geometry;
[faces,source]=reconstructEulerPrimitive(affine(geo.centerX,geo.centerR),rec,false);
expected=affine(geo.faceBX(:,1:end-1),geo.faceBR(:,1:end-1));
assert(max(abs(faces{3}-expected),[],'all')<1e-10);
expected=affine(geo.areaCenterX,geo.areaCenterR);assert(max(abs(source-expected(:,:,4)),[],'all')<1e-10);
% 极强跳变：缩小斜率使所有面的密度和压力保持正值，原单元值不能被改写。
q=affine(geo.centerX,geo.centerR);q(1:3,:,1)=1e-8;q(1:3,:,4)=1e-10;before=q;
[faces,source,d]=reconstructEulerPrimitive(q,rec,true);
assert(isequal(q,before) && d.limitedCellCount>0 && d.minimumLimiter>=0);
for k=1:4,assert(all(faces{k}(:,:,[1 4])>0,'all') && all(isfinite(faces{k}),'all'));end
assert(all(source>0,'all'));
% 平滑限幅保持常量；弱超调应小于独立有理函数上界 sqrt(epsilon²/8)。
constant=repmat(reshape([2,20,4,10],1,1,4),g.nx,g.nr,1);
[faces,source]=reconstructEulerPrimitive(constant,rec,true);
for k=1:4,assert(isequal(faces{k},constant));end
assert(all(source==10,'all'));
q=affine(geo.centerX,geo.centerR);[faces,~]=reconstructEulerPrimitive(q,rec,true);
lo=min(reshape(q,[],4),[],1);hi=max(reshape(q,[],4),[],1);
allowance=sqrt(max(rec.smoothingLengthCubed)/8).*max(abs(reshape(q,[],4)),[],1);
for k=1:4
 values=reshape(faces{k},[],4);
 assert(all(values>=lo-allowance-1e-12 & values<=hi+allowance+1e-12,'all'));
end
fprintf('Euler reconstruction passed: analytic weighted moments, translation, affine fields, source point and positive faces.\n');
end
function q=affine(x,r)
q=cat(3,2+.2*x+.1*r,20+3*x-2*r,4+2*x+3*r,10+2*x+4*r);
end
