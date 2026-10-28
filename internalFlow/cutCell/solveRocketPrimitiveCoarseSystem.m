function [correction,report]=solveRocketPrimitiveCoarseSystem(coarse,remaining)
%SOLVEROCKETPRIMITIVECOARSESYSTEM 固定原始变量粗像的最小二乘投影，不更新实际流场。
% 模块一：小Gram直接解后按原局部切线提升，检查细模型残差不增。
assert(strcmp(coarse.kind,'volume-primitive-tangent-fine-L2-minimum-residual-coarse-correction')&& ...
    strcmp(coarse.primitiveLiftSourceSHA256,rocketPrimitiveCoarseLiftIdentity),'PrimitiveCoarse:Source','原始变量粗基或提升来源不符。');
vector=remaining(:);rhs=coarse.image'*vector;c=coarse.upper\(coarse.upper'\rhs);
correction=liftRocketPrimitiveCoarse(coarse.baseQ,coarse.gas,coarse.scale,coarse.space,reshape(c.*coarse.columnLiftScale,coarse.space.groupCount,5));
predicted=vector-coarse.image*c;algebraic=norm(coarse.normalMatrix*c-rhs)/max(1,norm(rhs));
assert(algebraic<1e-8&&norm(predicted)<=norm(vector)+1e-10*max(1,norm(vector)),'PrimitiveCoarse:Projection','粗代数误差或模型不增检查失败。');
report=struct('coarseAlgebraicRelativeResidual',algebraic,'predictedRemainingNorm',norm(predicted),'fixedLinearMap',true, ...
    'actualOriginalNSEvaluations',0,'physicalSolution',false,'newProductionSolvedAreaM2',0);
end
