function [correction,report]=solveRocketMixedCoarseSystem(coarse,remaining)
%SOLVEROCKETMIXEDCOARSESYSTEM 以冻结混合粗基作细层范数最小二乘校正，不重新筛选右端。
% 模块一：逐列按保存身份提升，保留原保守量与原始变量方向的不同单位/来源。
% 模块二：粗代数误差和模型残差不增只是预条件检查，不能替代原NS验收。
assert(isequaln(coarse.liftSourceIdentity,rocketMixedCoarseLiftIdentity),'MixedCoarse:LiftSource','提升来源已改变，不使用旧混合像。');
right=coarse.image'*remaining(:);coefficients=coarse.upper\(coarse.upper'\right);correction=liftRocketMixedCoarse(coarse,coefficients);
predicted=remaining(:)-coarse.image*coefficients;before=norm(remaining(:));after=norm(predicted);
algebraic=norm(coarse.normalMatrix*coefficients-right)/max(1,norm(right));
assert(all(isfinite(correction),'all')&&algebraic<1e-8&&after<=before+1e-10*max(1,before),'MixedCoarse:Projection','混合代数误差或模型不增检查失败。');
report=struct('coarseAlgebraicRelativeResidual',algebraic,'coarsePredictedFineRelativeResidual',after/max(1,before), ...
    'coarseUnknowns',size(coarse.image,2),'coarseGroupCount',coarse.space.groupCount,'coarseAdditionalPrimitiveColumns',numel(coarse.selectedPrimitiveColumns), ...
    'fineModelResidualNonincreaseVerified',true,'originalNSLinearAccuracyClaimed',false,'actualOriginalNSEvaluations',0,'newProductionSolvedAreaM2',0);
end
