function [update,report]=solveRocketGalerkinCoarseSystem(coarse,rhs)
%SOLVEROCKETGALERKINCOARSESYSTEM 用固定粗层LU解限制后的五量右端并核对代数残差。
% 模块一：只解C系统，不把粗层小残差当原NS线性精度或流场稳态。
ng=coarse.space.groupCount;assert(isequal(size(rhs),[ng 5])&&all(isfinite(rhs),'all'),'RocketCoarse:RHS','粗层右端无效。');
vector=rhs(:);solution=coarse.upper\(coarse.lower\vector(coarse.rowPermutation));update=reshape(solution,ng,5);
error=coarse.matrix*solution-vector;relative=norm(error)/max(1,norm(vector));
assert(all(isfinite(update),'all')&&relative<1e-8,'RocketCoarse:LinearDefect','粗层代数解未通过明确残差检查。');
report=struct('coarseAlgebraicRelativeResidual',relative,'coarseUnknowns',5*ng,'coarseGroupCount',ng, ...
    'originalNSLinearAccuracyClaimed',false,'physicalTimeAdvanced',0,'newProductionSolvedAreaM2',0);
end
