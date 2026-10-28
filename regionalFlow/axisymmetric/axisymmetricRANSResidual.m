function [rate,info]=axisymmetricRANSResidual(Q,mesh,gas,boundary)
%AXISYMMETRICRANSRESIDUAL 真正五变量轴对称SA-noft2/NS耦合空间方程。
% 模块一：前四量使用原EOS；第五量W=rho*nuTilde，不虚构湍动能或裁负状态。
% 模块二：真实最短壁距、SA完整输运与NS有效应力/导热在同一次残差计算。
% 模块三：共享一次质量通量及完整五量边界/源账，不作分裂事后更新。
gas=axisymmetricFVGas(gas);
assert(isfield(gas,'transport'),'AxisSA:Gas','当前RANS切片要求显式Sutherland输运输入。');
assert(isnumeric(Q)&&isreal(Q)&&isequal(size(Q),[mesh.cellCount 5])&&all(isfinite(Q),'all')&&all(Q(:,5)>=0), ...
    'AxisSA:State','五量Q须有限且W非负，未裁值。');
P=axisymmetricFVState(Q(:,1:4),gas);N=Q(:,5)./P(:,1);
assert(all(isfinite(N)&N>=0),'AxisSA:State','nuTilde非有限或负值。');
distance=axisymmetricSAWallDistance(mesh,boundary);
[flow,info]=axisymmetricNavierStokesResidual(Q(:,1:4),mesh,gas,boundary,struct('nuTilde',N,'wallDistance',distance.distance));
sa=info.sa;flux=sa.faceFluxIntegralPerRadian;source=sa.cellSourceIntegral;
scalar=(-mesh.incidence*flux+source)./mesh.volumePerRadian;rate=[flow scalar];
info.faceFluxIntegralPerRadian=[info.faceFluxIntegralPerRadian flux];
info.cellSourceIntegral=[info.cellSourceIntegral source];
info.netBoundaryFlux=sum(info.faceFluxIntegralPerRadian(mesh.boundaryFaces,:),1);
info.sourceIntegral=sum(info.cellSourceIntegral,1);
info.conservationDefect=sum(rate.*mesh.volumePerRadian,1)+info.netBoundaryFlux-info.sourceIntegral;
info.sa.wallDistanceBinding=distance;
info.method=[info.method '; coupled fifth conservative variable rho*nuTilde; complete density-gradient correction'];
info.quantityOrder={'mass','axialMomentum','radialMomentum','totalEnergy','rhoNuTilde'};
info.quantityUnits={'kg/(s rad)','N/rad','N/rad','W/rad','kg m^2/(s^2 rad)'};
info.turbulenceImplemented=true;info.turbulenceValidated=false;
end
