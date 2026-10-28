function [left,right,audit]=cutCellLineFaceLinearization(QL,QR,normal,distance,gas)
%CUTCELLLINEFACELINEARIZATION 冻结波速/黏性功速度的两点面通量近似导数。
% 模块一：Rusanov导数只作预条件；正式残差仍用原完整HLLC/WLS/Fourier。
% 模块二：真实正法向质心距；Newton应力含三个速度及法向压缩项，导热含dT/dQ。
% 模块三：记录冻结lambda与黏性功速度，不能将本近似当作完整NS Jacobian。
N=size(QL,1);distance=distance(:);assert(isequal(size(QR),size(QL))&&numel(distance)==N&&all(distance>0&isfinite(distance)), ...
    'ShellLine:Distance','真实相邻质心法向距须为正。');
[left,HL,PL]=cutCellConservativeDerivatives(QL,normal,gas);
[right,HR,PR]=cutCellConservativeDerivatives(QR,normal,gas);
lambda=max(abs(sum(PL(:,2:4).*normal,2))+sqrt(gas.gamma*gas.R*PL(:,5)), ...
    abs(sum(PR(:,2:4).*normal,2))+sqrt(gas.gamma*gas.R*PR(:,5)));
left=.5*left;right=.5*right;
for a=1:5,left(a,a,:)=left(a,a,:)+reshape(.5*lambda,1,1,N);right(a,a,:)=right(a,a,:)-reshape(.5*lambda,1,1,N);end
velocity=.5*(PL(:,2:4)+PR(:,2:4));k=gas.mu*gas.gamma*gas.R/(gas.gamma-1)/gas.Pr;
for b=1:5
    for a=1:3
        dl=zeros(N,1);dr=dl;
        for c=1:3
            factor=gas.mu./distance.*(double(a==c)+normal(:,a).*normal(:,c)/3);
            dl=dl+factor.*reshape(HL(c,b,:),N,1);dr=dr+factor.*reshape(HR(c,b,:),N,1);
        end
        left(a+1,b,:)=left(a+1,b,:)+reshape(dl,1,1,N);
        right(a+1,b,:)=right(a+1,b,:)-reshape(dr,1,1,N);
        left(5,b,:)=left(5,b,:)+reshape(velocity(:,a).*dl,1,1,N);
        right(5,b,:)=right(5,b,:)-reshape(velocity(:,a).*dr,1,1,N);
    end
    left(5,b,:)=left(5,b,:)+reshape(k./distance.*reshape(HL(4,b,:),N,1),1,1,N);
    right(5,b,:)=right(5,b,:)-reshape(k./distance.*reshape(HR(4,b,:),N,1),1,1,N);
end
audit=struct('waveSpeedFrozen',lambda,'viscousWorkVelocityFrozen',velocity, ...
    'normalDistanceM',distance,'fullOriginalNSJacobian',false,'physicalDiscretizationReplaced',false);
end
