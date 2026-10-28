function P=axisymmetricFVState(Q,gas)
%AXISYMMETRICFVSTATE 从真实四守恒量得到[rho ux ur T p a]，无旋/不裁值。
gas=axisymmetricFVGas(gas);
assert(isnumeric(Q)&&isreal(Q)&&ismatrix(Q)&&size(Q,2)==4&&all(isfinite(Q),'all')&& ...
    all(Q(:,1)>0),'AxisFV:State','Q须为有限实数nc×4，密度为正。');
Q=double(Q);
u=Q(:,2:3)./Q(:,1);e=Q(:,4)./Q(:,1)-sum(u.^2,2)/2;T=e/gas.cv;
assert(all(isfinite(T)&T>=gas.temperatureRange(1)&T<=gas.temperatureRange(2)), ...
    'AxisFV:Temperature','内能/温度超出明确接受域，未裁剪状态。');
p=Q(:,1).*gas.R.*T;a=sqrt(gas.gamma*gas.R*T);P=[Q(:,1),u,T,p,a];
assert(all(isfinite(P),'all')&&all(p>0),'AxisFV:State','导出压力/声速非有限或非正。');
end
