function [rate,info]=rocketFiniteStepPatchResidual(Q,patch)
%ROCKETFINITESTEPPATCHRESIDUAL 真实内部控制体及两圈耦合的缩小复现算子。
% 公式与原NS的内部精确Godunov/WLS/Newton/Fourier/算术黏性功一致。
% 核心外状态为诊断夹具；不新增真实边界条件，不授予整箭求解资格。
g=patch.gas;rho=Q(:,1);u=Q(:,2:4)./rho;internal=Q(:,5)-.5*rho.*sum(u.^2,2);
assert(all(isfinite(Q),'all')&&all(rho>0),'CutNS:State','密度或状态无效，不裁剪。');
assert(all(internal>0),'CutNS:InternalEnergy','内能非正，不裁剪。');
T=internal./(rho*g.R/(g.gamma-1));
assert(all(T>=g.temperatureRange(1)&T<=g.temperatureRange(2)),'CutNS:Temperature','温域保持80至2000K。');
phi=[u T];n=size(Q,1);A=zeros(3,3,n);B=zeros(3,4,n);centroid=patch.centroid;
for f=1:numel(patch.gradientArea)
    m=patch.gradientMinus(f);p=patch.gradientPlus(f);d=centroid(p,:)-centroid(m,:);v=phi(p,:)-phi(m,:);
    weight=patch.gradientArea(f)/dot(d,d);aa=weight*(d'*d);bb=weight*(d'*v);
    A(:,:,m)=A(:,:,m)+aa;A(:,:,p)=A(:,:,p)+aa;B(:,:,m)=B(:,:,m)+bb;B(:,:,p)=B(:,:,p)+bb;
end
gradient=nan(4,3,n);
for id=patch.gradientCells(:)'
    assert(rcond(A(:,:,id))>1e-13,'CutNS:GradientRank','原梯度约束退化。');
    gradient(:,:,id)=(A(:,:,id)\B(:,:,id))';
end
m=patch.fluxMinus;p=patch.fluxPlus;normal=patch.fluxNormal;
[convection,audit]=cutCellExactEulerFlux(Q(m,:),Q(p,:),normal,g);
viscous=zeros(size(convection));k=g.mu*g.gamma*g.R/(g.gamma-1)/g.Pr;
for f=1:numel(m)
    d=centroid(p(f),:)-centroid(m(f),:);v=phi(p(f),:)-phi(m(f),:);
    G=.5*(gradient(:,:,m(f))+gradient(:,:,p(f)));G=G+((v(:)-G*d')/dot(d,d))*d;
    D=G(1:3,:);stress=g.mu*(D+D'-(2/3)*trace(D)*eye(3));traction=(stress*normal(f,:)')';
    velocity=.5*(u(m(f),:)+u(p(f),:));viscous(f,:)=[0 traction dot(velocity,traction)+k*dot(G(4,:),normal(f,:))];
end
flux=convection-viscous;integrated=zeros(n,5);
% 与生产NS保持原面顺序，以同一双精度加法次序复现核心行。
for f=1:numel(m)
    value=patch.fluxArea(f)*flux(f,:);integrated(m(f),:)=integrated(m(f),:)-value;integrated(p(f),:)=integrated(p(f),:)+value;
end
core=patch.core;rate=integrated(core,:)./patch.volume(core);
outSign=double(ismember(m,core))-double(ismember(p,core));net=sum(outSign.*patch.fluxArea.*flux,1);
balance=sum(integrated(core,:),1)+net;denom=max(sum(abs(integrated(core,:)),1),realmin);
info=struct('gradient',gradient(:,:,patch.gradientCells),'faceFlux',flux,'convectiveFlux',convection, ...
    'viscousFlux',viscous,'sampleRegion',audit.sampleRegion,'minimumDensity',min(rho),'temperatureRange',[min(T) max(T)], ...
    'spatialConservationRelative',max(abs(balance)./denom),'boundaryNetFlux',net, ...
    'patchSteadyFluxRatio',max(abs(net)./max(sum(abs(outSign.*patch.fluxArea.*flux),1),realmin)), ...
    'physicalSolution',false,'newProductionSolvedAreaM2',0);
end
