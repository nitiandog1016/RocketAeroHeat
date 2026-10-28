function [A,H,P]=cutCellConservativeDerivatives(Q,normal,gas)
%CUTCELLCONSERVATIVEDERIVATIVES 定比热五量Euler法向通量及速度/温度的保守态导数。
% 模块一：完整三个动量分量；P=[rho,u,v,w,T]，H=d[u,v,w,T]/dQ。
% 模块二：A=dF_Euler(Q,n)/dQ，法向来自真实三维共享面，不投影成轴对称。
% 模块三：只供显式标明的近似隐式预条件；不改正式HLLC/热梯度空间方程。
N=size(Q,1);assert(isequal(size(Q),[N 5])&&isequal(size(normal),[N 3])&& ...
    all(isfinite(Q),'all')&&all(isfinite(normal),'all')&&all(abs(sum(normal.^2,2)-1)<1e-10), ...
    'ShellLine:State','五量/单位法向合同不符。');
rho=Q(:,1);assert(all(rho>0),'CutNS:State','密度非正，不裁剪。');
u=Q(:,2:4)./rho;u2=sum(u.^2,2);cv=gas.R/(gas.gamma-1);e=Q(:,5)-.5*rho.*u2;
assert(all(e>0),'CutNS:InternalEnergy','内能非正，不裁剪。');T=e./(rho*cv);
assert(all(T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)),'CutNS:Temperature','温度越过原物性域。');
P=[rho,u,T];H=zeros(4,5,N);A=zeros(5,5,N);gm=gas.gamma-1;
dp=[.5*gm*u2,-gm*u,repmat(gm,N,1)];un=sum(u.*normal,2);enthalpy=(Q(:,5)+gm*e)./rho;
for a=1:3
    H(a,1,:)=reshape(-u(:,a)./rho,1,1,N);H(a,a+1,:)=reshape(1./rho,1,1,N);
    A(1,a+1,:)=reshape(normal(:,a),1,1,N);
    A(a+1,1,:)=reshape(-u(:,a).*un+normal(:,a).*dp(:,1),1,1,N);
    for b=1:3
        A(a+1,b+1,:)=reshape(double(a==b)*un+u(:,a).*normal(:,b)+normal(:,a).*dp(:,b+1),1,1,N);
    end
    A(a+1,5,:)=reshape(gm*normal(:,a),1,1,N);
end
dT=[(.5*u2-cv*T)./(rho*cv),-u./(rho*cv),1./(rho*cv)];
H(4,:,:)=permute(dT,[3 2 1]);
A(5,1,:)=reshape(un.*(dp(:,1)-enthalpy),1,1,N);
for b=1:3,A(5,b+1,:)=reshape(un.*dp(:,b+1)+enthalpy.*normal(:,b),1,1,N);end
A(5,5,:)=reshape(gas.gamma*un,1,1,N);
assert(all(isfinite(A),'all')&&all(isfinite(H),'all'),'ShellLine:State','保守态导数非有限。');
end
