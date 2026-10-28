function [correction,report]=solveRocketLeastSquaresCoarseSystem(coarse,remaining)
%SOLVEROCKETLEASTSQUARESCOARSESYSTEM 在当前细层原收敛范数下作固定线性最小二乘校正。
% 模块一：小Gram解只决定粗基系数，细网格和物理边界不变。
% 模块二：检查投影后模型残差不增，不把该检查当非线性NS或热流验收。
ng=coarse.space.groupCount;n=numel(coarse.space.fineVolume);
assert(isequal(size(remaining),[n 5])&&all(isfinite(remaining),'all'),'RocketCoarse:RHS','粗像校正右端无效。');
vector=remaining(:);rhs=coarse.image'*vector;coefficients=coarse.upper\(coarse.upper'\rhs);
correction=coarse.space.prolongation*reshape(coefficients.*coarse.columnLiftScale,ng,5);
projected=vector-coarse.image*coefficients;before=norm(vector);after=norm(projected);
algebraic=norm(coarse.normalMatrix*coefficients-rhs)/max(1,norm(rhs));
assert(all(isfinite(correction),'all')&&algebraic<1e-8&&after<=before+1e-10*max(1,before), ...
    'RocketCoarse:LeastSquares','粗像代数残差或细层模型不增检查失败，不裁结果。');
report=struct('coarseAlgebraicRelativeResidual',algebraic,'coarsePredictedFineRelativeResidual',after/max(1,before), ...
    'coarseUnknowns',5*ng,'coarseGroupCount',ng,'fineModelResidualNonincreaseVerified',true, ...
    'originalNSLinearAccuracyClaimed',false,'actualOriginalNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
