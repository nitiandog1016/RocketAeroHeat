function Q=axisymmetricFVConservative(P,gas)
%AXISYMMETRICFVCONSERVATIVE 将明确[rho ux ur T]输入变为四个守恒量。
gas=axisymmetricFVGas(gas);
assert(isnumeric(P)&&isreal(P)&&ismatrix(P)&&size(P,2)==4&&all(isfinite(P),'all')&& ...
    all(P(:,1)>0),'AxisFV:State','原始量须为有限实数nc×4。');
P=double(P);
Q=[P(:,1),P(:,1).*P(:,2:3),P(:,1).*(gas.cv*P(:,4)+sum(P(:,2:3).^2,2)/2)];
axisymmetricFVState(Q,gas);
end
