function out=invertRocketPrimitiveUpdateMetric(metric,value)
%INVERTROCKETPRIMITIVEUPDATEMETRIC 精确反解原始变量切向度量，不改变非线性候选坐标。
assert(~isfield(metric,'peakRowWeights')&&isequal(size(value),[metric.cellCount 5]), ...
    'PrimitiveMetric:Inverse','仅接受未加峰值权重的真实五量度量。');
out=zeros(size(value));out(:,1)=value(:,1)./metric.diagonal(:,1);
out(:,2:4)=(value(:,2:4)-metric.densityToVelocity.*out(:,1))./metric.diagonal(:,2:4);
out(:,5)=(value(:,5)-metric.densityToTemperature.*out(:,1)-sum(metric.momentumToTemperature.*out(:,2:4),2))./metric.diagonal(:,5);
assert(all(isfinite(out),'all'),'PrimitiveMetric:Inverse','度量逆非有限。');
end
