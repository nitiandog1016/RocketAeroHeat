function sample=parabolicShockSample(shape,q)
%PARABOLICSHOCKSAMPLE 同一迎风坐标抛物激波的三维位置和几何导数。
% 模块一：frame第一列沿均匀来流，另外两列为右手横向基。
% 模块二：x=0.5*q'*K*q，K为正定曲率参数；几何仍待与箭体反算匹配。
validateattributes(q,{'double'},{'real','finite','size',[2 1]});
validateattributes(shape.origin,{'double'},{'real','finite','size',[3 1]});
validateattributes(shape.frame,{'double'},{'real','finite','size',[3 3]});
validateattributes(shape.K,{'double'},{'real','finite','size',[2 2]});
assert(norm(shape.frame'*shape.frame-eye(3),'fro')<1e-11&&det(shape.frame)>0, ...
    'ShockLayer:Frame','需要正交右手迎风坐标，不自动修正。');
assert(norm(shape.K-shape.K','fro')<1e-12*norm(shape.K,'fro')&&min(eig(shape.K))>0, ...
    'ShockLayer:Shape','抛物激波K须对称正定。');
g=shape.K*q;
sample=struct('position',shape.origin+shape.frame*[.5*q'*shape.K*q;q], ...
    'firstDerivatives',shape.frame*[g';eye(2)], ...
    'secondDerivatives',shape.frame*[shape.K(1,1) shape.K(1,2) shape.K(2,2);zeros(2,3)]);
end
