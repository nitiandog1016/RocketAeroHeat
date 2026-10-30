function sample=deformedConicShockSample(shape,q)
%DEFORMEDCONICSHOCKSAMPLE 二次曲线共同激波上的全局四次形变及精确一/二阶几何。
% 模块一：原二次曲线几何保持；没有形变时位置和导数逐值沿用原结果。
% 模块二：sum(w_j*(d_j'*q)^4/(4*L^3))是整体几何场，不按面元拟合压力。
% 模块三：返回迎风图的梯度/Hessian，供一般投影质量连续式使用。
base=ellipticConicShockSample(shape,q);sample=base;
localA=shape.frame'*base.firstDerivatives;localB=shape.frame'*base.secondDerivatives;
gradient=localA(1,:)';H=[localB(1,1) localB(1,2);localB(1,2) localB(1,3)];
if isfield(shape,'deformationWeights')&&~isempty(shape.deformationWeights)
    weights=shape.deformationWeights(:);directions=shape.deformationDirections;L=shape.deformationReferenceLengthM;
    validateattributes(weights,{'double'},{'real','finite','vector'});
    validateattributes(directions,{'double'},{'real','finite','size',[2 numel(weights)]});
    validateattributes(L,{'double'},{'real','finite','scalar','positive'});
    assert(all(abs(vecnorm(directions,2,1)-1)<1e-11),'ShockLayer:Shape','全局形变方向须明确为单位向量。');
    r=directions'*q;value=sum(weights.*r.^4)/(4*L^3);
    extraGradient=directions*(weights.*r.^3)/L^3;
    extraHessian=directions*diag(3*weights.*r.^2/L^3)*directions';
    sample.position=base.position+shape.frame(:,1)*value;
    sample.firstDerivatives=base.firstDerivatives+shape.frame(:,1)*extraGradient';
    sample.secondDerivatives=base.secondDerivatives+shape.frame(:,1)* ...
        [extraHessian(1,1) extraHessian(1,2) extraHessian(2,2)];
    gradient=gradient+extraGradient;H=H+extraHessian;
end
sample.graphGradient=gradient;sample.graphHessian=H;
sample.graphMethod='共同二次曲线 + 全局方向四次形变；显式一/二阶导数';
end
