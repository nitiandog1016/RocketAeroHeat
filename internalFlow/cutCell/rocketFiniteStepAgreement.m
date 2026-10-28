function evidence=rocketFiniteStepAgreement(baseRate,candidateRate,originalJv,scale,alpha)
%ROCKETFINITESTEPAGREEMENT 原NS有限步与同基态一阶预测的一致性保护。
% 不改物理残差或任何温域/正性/方向/稳态门槛；仅更严格地拒绝非线性候选。
% 固定半预测变化的误差预算，兼顾原最大局部残差不升；不裁剪或冻结波系。
validateattributes(alpha,{'numeric'},{'real','finite','positive','scalar'});
assert(isequal(size(baseRate),size(candidateRate),size(originalJv))&&numel(scale)==5, ...
    'FiniteStep:Shape','原NS率和方向作用尺寸必须一致。');
base=baseRate./scale;candidate=candidateRate./scale;prediction=alpha*originalJv;
remainder=candidate-base-prediction;denominator=norm(prediction,'fro');
ratio=norm(remainder,'fro')/max(denominator,realmin);
baseMax=max(abs(base),[],'all');candidateMax=max(abs(candidate),[],'all');
evidence=struct('kind','original-NS-finite-step-model-agreement','relativeNonlinearRemainder',ratio, ...
    'maximumRemainderRatio',.5,'baseMaximumResidual',baseMax,'candidateMaximumResidual',candidateMax, ...
    'maximumResidualNonIncreasing',candidateMax<=baseMax,'modelAgreementMet',isfinite(ratio)&&ratio<=.5, ...
    'acceptable',isfinite(ratio)&&ratio<=.5&&candidateMax<=baseMax, ...
    'physicalEquationChanged',false,'stateClipped',false);
end
