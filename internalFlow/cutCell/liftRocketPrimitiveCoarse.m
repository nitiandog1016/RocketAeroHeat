function direction=liftRocketPrimitiveCoarse(Q,gas,scale,space,coefficients)
%LIFTROCKETPRIMITIVECOARSE 将分组常数[rho,u,v,w,T]扰动按当前Q的精确切线提升为缩放五量。
% 模块一：保留密度、动量及总能量的一致耦合，局部速度和温度可随基态变化。
% 模块二：本函数是固定基态线性映射，不生成非线性候选、不裁剪状态。
assert(isequal(size(coefficients),[space.groupCount 5])&&all(isfinite(coefficients),'all'),'PrimitiveCoarse:Coefficients','原始变量粗系数非法。');
fine=space.prolongation*coefficients;rho=Q(:,1);u=Q(:,2:4)./rho;cv=gas.R/(gas.gamma-1);
specificEnergy=Q(:,5)./rho;drho=fine(:,1);du=fine(:,2:4);dT=fine(:,5);
dQ=[drho,u.*drho+rho.*du,specificEnergy.*drho+rho.*sum(u.*du,2)+rho*cv.*dT];
direction=dQ./scale;assert(all(isfinite(direction),'all'),'PrimitiveCoarse:Lift','原始变量切线提升非有限。');
end
