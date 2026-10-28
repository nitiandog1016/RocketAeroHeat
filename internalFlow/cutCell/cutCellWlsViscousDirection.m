function [value,detail]=cutCellWlsViscousDirection(Q,directionU,baseInfo,cut,gas,boundary,scale,cache,baseBoundary,boundaryDirection)
%CUTCELLWLSVISCOUSDIRECTION 原完整WLS、位移修正和静止墙黏性通量的方向导数。
% 模块一：保守扰动经EOS转为速度/温度扰动，开放边界样本提升可显式计入。
% 模块二：保留完整位移梯度、Newton应力及Fourier项，墙采用当前零切向约束。
% 模块三：黏性功分开保存U*d牵引和dU*原牵引；后者不能按冻结功速度省略。
% 本函数只返回黏性部分；没有原NS非线性更新或面积资格。
assert(isequaln(cache.geometry,cut)&&isequal(cache.incoming,baseBoundary.incoming)&&isequal(cache.wall,baseBoundary.wall), ...
    'CutDirection:Cache','几何或入流分支改变，不能沿用旧WLS缓存。');
assert(strcmp(cutCellViscousWorkPolicy(boundary),'arithmetic')&& ...
    strcmp(cutCellWallVelocityPolicy(boundary.solidWall),'dirichlet-tangent')&& ...
    strcmp(cutCellWallThermalPolicy(boundary.solidWall),'dirichlet-tangent')&&all(boundary.solidWall.velocity==0), ...
    'CutDirection:Policy','当前方向模块要求算术内面功速度及显式静止Dirichlet墙。');
scale=reshape(scale,1,5);dq=directionU.*scale;n=size(Q,1);cv=gas.R/(gas.gamma-1);k=gas.mu*gas.gamma*cv/gas.Pr;
rho=Q(:,1);u=Q(:,2:4)./rho;T=(Q(:,5)-.5*rho.*sum(u.^2,2))./(rho*cv);
assert(isequal(size(Q),[n 5])&&isequal(size(directionU),size(Q))&&all(isfinite(dq),'all')&&all(rho>0));
assert(max(abs([rho u T]-baseInfo.primitive),[],'all')<=1e-13*max(1,max(abs(baseInfo.primitive),[],'all')), ...
    'CutDirection:Base','原保守态与保存的基态primitive不一致。');
du=(dq(:,2:4)-u.*dq(:,1))./rho;
dT=(dq(:,5)-sum(u.*dq(:,2:4),2)+(.5*sum(u.^2,2)-cv*T).*dq(:,1))./(rho*cv);dphi=[du dT];
f=cut.openFaces;w=cut.walls;nf=numel(f.area);nw=numel(w.area);v=cut.cells.volume(:);
given=zeros(nf,5);sensitivity=false;
if ~isempty(boundaryDirection)
    assert(isequal(boundaryDirection.incoming,baseBoundary.incoming)&&isequal(boundaryDirection.wall,baseBoundary.wall));
    given=boundaryDirection.primitiveDerivative;sensitivity=boundaryDirection.boundarySensitivityIncluded;
end
assert(isequal(size(given),[nf 5])&&all(isfinite(given),'all'));

%% 模块四：WLS内部稀疏梯度加上非冻结的开放边界样本方向。
gradient=zeros(4,3,n);ids=cache.incomingFaceIDs;owner=cache.incomingOwner;
for d=1:3
    selected=cache.gradientOperator{d}*dphi;
    for a=1:4
        selected(:,a)=selected(:,a)+accumarray(owner,cache.boundaryCoefficientXYZ(:,d).*given(ids,a+1),[n 1]);
    end
    gradient(:,d,:)=reshape(selected',4,1,n);
end
stressFace=zeros(nf,3);heatFace=zeros(nf,1);workGradientFace=heatFace;workVelocityFace=heatFace;
active=find(f.area>0);
for first=1:2048:numel(active)
    ids=active(first:min(first+2047,numel(active)));mi=f.cellMinus(ids);pi=f.cellPlus(ids);inside=mi>0&pi>0;
    G=zeros(4,3,numel(ids));velocity=zeros(numel(ids),3);dvelocity=velocity;
    if any(inside)
        ii=find(inside);m=mi(ii);p=pi(ii);delta=cut.cells.centroid(p,:)-cut.cells.centroid(m,:);
        G(:,:,ii)=correct(.5*(gradient(:,:,m)+gradient(:,:,p)),delta,dphi(p,:)-dphi(m,:));
        velocity(ii,:)=.5*(u(m,:)+u(p,:));dvelocity(ii,:)=.5*(du(m,:)+du(p,:));
    end
    if any(~inside)
        ii=find(~inside);faces=ids(ii);cells=max(mi(ii),pi(ii));G(:,:,ii)=gradient(:,:,cells);
        incoming=baseBoundary.incoming(faces);pick=ii(incoming);selected=faces(incoming);cellsIn=cells(incoming);
        if ~isempty(pick)
            delta=f.centroid(selected,:)-cut.cells.centroid(cellsIn,:);
            G(:,:,pick)=correct(G(:,:,pick),delta,given(selected,2:5)-dphi(cellsIn,:));
        end
        pick=ii(~incoming);selected=faces(~incoming);
        if ~isempty(pick)
            normal=permute(f.normal(selected,:),[2 3 1]);normalGradient=pagemtimes(G(:,:,pick),normal);
            G(:,:,pick)=G(:,:,pick)-pagemtimes(normalGradient,permute(f.normal(selected,:),[3 2 1]));
        end
        velocity(ii,:)=baseBoundary.primitive(faces,2:4);dvelocity(ii,:)=given(faces,2:4);
    end
    normal=f.normal(ids,:);traction=force(G,normal);stressFace(ids,:)=traction;
    heatFace(ids)=k*sum(reshape(G(4,:,:),3,[])'.*normal,2);
    workGradientFace(ids)=sum(velocity.*traction,2);
    workVelocityFace(ids)=sum(dvelocity.*baseInfo.viscousFaceFlux(ids,2:4),2);
end

%% 模块五：真实正法向距及固定墙值的导数，静止墙黏性功严格为零。
stressWall=zeros(nw,3);heatWall=zeros(nw,1);wallIDs=find(w.ownerCell>0);
for first=1:2048:numel(wallIDs)
    ids=wallIDs(first:min(first+2047,numel(wallIDs)));owner=w.ownerCell(ids);normal=w.normal(ids,:);
    dn=sum((w.centroid(ids,:)-cut.cells.centroid(owner,:)).*normal,2);assert(all(dn>0),'CutDirection:Wall','真实墙距非正。');
    G=zeros(4,3,numel(ids));
    for a=1:3,G(a,:,:)=reshape((-du(owner,a)./dn.*normal)',1,3,[]);end
    G(4,:,:)=reshape((-dT(owner)./dn.*normal)',1,3,[]);
    stressWall(ids,:)=force(G,normal);heatWall(ids)=k*sum(reshape(G(4,:,:),3,[])'.*normal,2);
end
momentum=zeros(n,5);momentum(:,2:4)=divergence(stressFace,stressWall);heat=zeros(n,5);heat(:,5)=divergence(heatFace,heatWall);
workGradient=zeros(n,5);workGradient(:,5)=divergence(workGradientFace,zeros(nw,1));
workVelocity=zeros(n,5);workVelocity(:,5)=divergence(workVelocityFace,zeros(nw,1));
momentum=momentum./scale;heat=heat./scale;workGradient=workGradient./scale;workVelocity=workVelocity./scale;
value=momentum+heat+workGradient+workVelocity;assert(all(isfinite(value),'all'),'CutDirection:Nonfinite','WLS黏性方向非有限。');
detail=struct('kind','selected-branch-full-WLS-static-wall-and-arithmetic-work-viscous-direction', ...
    'momentumU',momentum,'heatU',heat,'workGradientU',workGradient,'workVelocityU',workVelocity,'gradientDerivative',gradient, ...
    'faceViscousFluxDerivative',[zeros(nf,1) stressFace heatFace+workGradientFace+workVelocityFace], ...
    'wallViscousFluxDerivative',[zeros(nw,1) stressWall heatWall], ...
    'boundarySensitivityIncluded',sensitivity,'solidWallVelocityOrTemperatureDifferentiated',false, ...
    'oldStressMatrixReused',false,'fullNSJacobianClaimed',false,'newOriginalNSEvaluations',0,'newProductionSolvedAreaM2',0);
    function corrected=correct(G,delta,difference)
        d2=sum(delta.^2,2);assert(all(d2>0),'CutDirection:Distance','位移修正距离为零。');
        error=difference'-reshape(pagemtimes(G,permute(delta,[2 3 1])),4,[]);
        corrected=G+pagemtimes(reshape(error./d2',4,1,[]),permute(delta,[3 2 1]));
    end
    function traction=force(G,normal)
        V=G(1:3,:,:);S=gas.mu*(V+pagetranspose(V));trace=V(1,1,:)+V(2,2,:)+V(3,3,:);
        for a=1:3,S(a,a,:)=S(a,a,:)-(2/3)*gas.mu*trace;end
        traction=reshape(pagemtimes(S,permute(normal,[2 3 1])),3,[])';
    end
    function out=divergence(faceFlux,wallFlux)
        out=zeros(n,size(faceFlux,2));m=f.cellMinus;p=f.cellPlus;left=f.area>0&m>0;right=f.area>0&p>0;wall=w.ownerCell>0;
        for a=1:size(faceFlux,2)
            out(:,a)=accumarray(m(left),f.area(left)./v(m(left)).*faceFlux(left,a),[n 1])- ...
                accumarray(p(right),f.area(right)./v(p(right)).*faceFlux(right,a),[n 1])+ ...
                accumarray(w.ownerCell(wall),w.area(wall)./v(w.ownerCell(wall)).*wallFlux(wall,a),[n 1]);
        end
    end
end
