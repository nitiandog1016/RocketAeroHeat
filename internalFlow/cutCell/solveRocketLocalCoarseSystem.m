function [correction,report]=solveRocketLocalCoarseSystem(coarse,remaining)
%SOLVEROCKETLOCALCOARSESYSTEM 在冻结局部像中最小化细层剩余，不按右端重选方向。
% 模块一：小Gram解与逐列提升组成固定线性预条件，不生成非线性流场解。
% 模块二：模型残差不增不能代替原NS方向、温域、下降或稳态验收。
assert(isequaln(coarse.liftSourceIdentity,rocketLocalCoarseLiftIdentity),'LocalCoarse:Source','局部提升来源已改变。');
assert(isequal(size(remaining),size(coarse.baseQ))&&all(isfinite(remaining),'all'),'LocalCoarse:RHS','剩余形状或数值非法。');
right=coarse.image'*remaining(:);weights=coarse.upper\(coarse.upper'\right);
correction=liftRocketLocalCoarse(coarse,weights);before=norm(remaining(:));after=norm(remaining(:)-coarse.image*weights);
error=norm(coarse.normalMatrix*weights-right)/max(1,norm(right));
assert(error<1e-8&&all(isfinite(correction),'all')&&after<=before+1e-10*max(1,before),'LocalCoarse:Projection','局部投影误差或不增检查失败。');
report=struct('coarseAlgebraicRelativeResidual',error,'coarsePredictedFineRelativeResidual',after/max(1,before), ...
    'coarseUnknowns',size(coarse.image,2),'coarseGroupCount',coarse.space.groupCount, ...
    'coarseRetainedMixedColumns',numel(coarse.retainedMixedColumns),'coarseLocalConservativeColumns',numel(coarse.localCellIDs), ...
    'fineModelResidualNonincreaseVerified',true,'originalNSLinearAccuracyClaimed',false,'actualOriginalNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
