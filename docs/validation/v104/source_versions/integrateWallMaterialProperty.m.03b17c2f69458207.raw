function value = integrateWallMaterialProperty(material,property,lower,upper)
%INTEGRATEWALLMATERIALPROPERTY 精确积分分段线性材料物性。
% specificHeat 对温度积分给单位质量储能变化 [J/kg]；
% thermalConductivity 积分给 Kirchhoff 势差 [W/m]。
% 逐段积分使用温差乘区间平均物性，避免近温度大焓相减的消减误差。
% 上下限可反向、可含 NaN；其他非法值和表外温度由统一模块拒绝。
if ~ismember(property,{'specificHeat','thermalConductivity'})
    error('RocketAeroHeat:Material:Property','只支持比热或导热系数积分。');
end
if isscalar(lower), lower=lower+zeros(size(upper)); end
if isscalar(upper), upper=upper+zeros(size(lower)); end
if ~isequal(size(lower),size(upper))
    error('RocketAeroHeat:Material:Size','材料物性积分上下限尺寸不匹配。');
end
wallMaterialProperties(material,lower); wallMaterialProperties(material,upper);
if ~isfield(material,'propertyTable') || isempty(material.propertyTable)
    value=material.(property)*(upper-lower); return;
end
t=material.propertyTable.temperature(:);
y=material.propertyTable.(property)(:);
lo=min(lower,upper); hi=max(lower,upper);
value=zeros(size(upper));
for i=1:numel(t)-1
    a=max(lo,t(i)); b=min(hi,t(i+1));
    inside=b>a;
    slope=(y(i+1)-y(i))/(t(i+1)-t(i));
    value(inside)=value(inside)+(b(inside)-a(inside)).* ...
        (y(i)+slope*((a(inside)+b(inside))/2-t(i)));
end
value=value.*sign(upper-lower);
value(isnan(lower)|isnan(upper))=NaN;
end
