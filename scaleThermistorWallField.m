function [rows,scale,deltaT,wallK]=scaleThermistorWallField(rows,b,airK,referenceWallK,referenceOutW)
% 模块一：以完整等温试件的空气侧热导求解集总壁温，或读取给定壁温。
referenceDelta=referenceWallK-airK;
assert(referenceDelta>0 && referenceOutW>0 && isfinite(referenceOutW));
if isfield(b,'netHeatingPowerW')
    power=b.netHeatingPowerW;
    assert(isnumeric(power) && isreal(power) && isscalar(power) && isfinite(power) && ...
        power>=0 && power<=referenceOutW,'RAH:Bead:Power','净对流功率超出同网格参照场范围。');
    scale=power/referenceOutW;deltaT=scale*referenceDelta;wallK=airK+deltaT;
else
    wallK=b.wallTemperatureK;deltaT=wallK-airK;scale=deltaT/referenceDelta;
end
% 模块二：恒物性及给定流场下，壁反力和首层超温随温差线性变化。
if scale==0
    rows.heatFluxWm2(:)=0;rows.nearWallCellTemperatureK(:)=airK;
elseif scale~=1
    rows.heatFluxWm2=scale*rows.heatFluxWm2;
    rows.nearWallCellTemperatureK=airK+scale*(rows.nearWallCellTemperatureK-airK);
end
end
