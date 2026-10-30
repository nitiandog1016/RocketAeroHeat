function sample=ellipticConicShockSample(shape,q)
%ELLIPTICCONICSHOCKSAMPLE 椭圆截面的纵向二次曲线共同激波，B=0为原抛物面。
% 模块一：F+B*F^2/2=s，s=q'*K*q/2；稳定根保留原量纲B[1/m]。
% 模块二：显式求F'、F''及二阶几何；不以多项式补丁拟合独立面元压力。
if ~isfield(shape,'profileB')||shape.profileB==0
    sample=parabolicShockSample(shape,q);return;
end
base=parabolicShockSample(shape,q);
validateattributes(shape.profileB,{'double'},{'real','finite','scalar'});
g=shape.K*q;s=.5*q'*shape.K*q;D=1+2*shape.profileB*s;
assert(D>0,'ShockLayer:Shape','纵向二次曲线在当前半径不再单值，不能裁剪根。');
root=sqrt(D);value=2*s/(1+root);first=1/root;second=-shape.profileB/root^3;
H=first*shape.K+second*(g*g');
sample=struct('position',shape.origin+shape.frame*[value;q], ...
    'firstDerivatives',shape.frame*[first*g';eye(2)], ...
    'secondDerivatives',shape.frame*[H(1,1) H(1,2) H(2,2);zeros(2,3)]);
assert(all(isfinite(base.position)),'ShockLayer:Shape','基本几何验证失败。');
end
