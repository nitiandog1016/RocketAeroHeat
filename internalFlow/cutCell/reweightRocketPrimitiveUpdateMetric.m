function [metric,report]=reweightRocketPrimitiveUpdateMetric(base,reference,options)
%REWEIGHTROCKETPRIMITIVEUPDATEMETRIC 用冻结参考方向加权少数很大的局部原始变量更新。
% 模块一：w=sqrt(1+(W*v_ref/transition)^2)，逐分量放大局部峰值的软惩罚。
% 模块二：权重在一次LSQR内固定，完整正向/转置保持线性，不随右端或迭代重筛。
% 模块三：最大权重只是数值权重上界，不是状态裁剪或更新硬上限；实际状态仍须原NS验收。
if nargin<3,options=struct;end
defaults=struct('transitionFraction',10,'maximumRowWeight',256);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'PrimitivePeak:Options','未知局部峰值权重选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
validateattributes(options.transitionFraction,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(options.maximumRowWeight,{'numeric'},{'real','finite','scalar','>=',1});
z=applyRocketPrimitiveUpdateMetric(base,reference);weights=min(options.maximumRowWeight,hypot(1,z/options.transitionFraction));
metric=base;metric.diagonal=base.diagonal.*weights;
metric.densityToVelocity=base.densityToVelocity.*weights(:,2:4);
metric.densityToTemperature=base.densityToTemperature.*weights(:,5);
metric.momentumToTemperature=base.momentumToTemperature.*weights(:,5);
report=struct('method','frozen-local-primitive-peak-reweighted-quadratic-penalty','options',options, ...
    'weightRange',[min(weights,[],'all') max(weights,[],'all')],'weightedComponents',sum(weights>sqrt(2),'all'), ...
    'upperWeightComponents',sum(weights==options.maximumRowWeight,'all'), ...
    'fixedWithinLinearSolve',true,'hardStateUpdateBoundEnforced',false,'physicalStateClipping',false, ...
    'referenceMaximumRelativeComponents',max(abs(z),[],1),'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0);
metric.peakWeighting=report;metric.peakRowWeights=weights;
end
