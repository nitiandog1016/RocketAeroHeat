function result=thermistorBeadMeanConvection(c,wallK)
%THERMISTORBEADMEANCONVECTION 按坪内、佐藤式（28）估算 B 型珠部的平均强迫对流。
% 输入是带线 B 珠的已确认工况和给定壁温，输出功率以珠部向空气放热为正。

% 模块一：绑定实物表 1 的几何与公式口径。L 为 2V/F，F 为论文名义珠部面积。
result=struct('available',false,'reason','','method','坪内、佐藤（1959）式（28），珠部平均强迫对流', ...
    'doi','10.1299/kikai1938.25.810','scope','珠部含接头的积分比较，引线及线端单列', ...
    'lengthM',0.3430e-3,'areaM2',3.094e-6,'eccentricity',0.53,'tableNu0',0.664, ...
    'speedMps',c.speedMps,'airTemperatureK',c.airK,'wallTemperatureK',wallK, ...
    'meanHWm2K',NaN,'outwardW',NaN);
if ~strcmp(c.caseKind,'thermistor-B-axial-two-leads') || c.sourceAreaM2~=result.areaM2
    result.reason='试件与表 1 的 B 型热敏珠不匹配。';return;
end
deltaT=wallK-c.airK;
if ~isscalar(deltaT) || ~isreal(deltaT) || ~isfinite(deltaT) || deltaT<=0
    result.reason='该加热试件的实验范围要求正的珠面与来流温差。';return;
end

% 模块二：按论文规定在膜温取空气物性，使用本项目已核验的解析输运入口。
result.filmTemperatureK=c.airK+deltaT/2;
try
    air=phaseAwareAirState('PT',c.airPa,result.filmTemperatureK,'lemmonAnalytic');
catch ex
    if ~startsWith(ex.identifier,'RAH:PhaseAir:'),rethrow(ex);end
    result.reason=['空气物性不可用。' ex.message];return;
end
nu=air.dynamicViscosityPas/air.densityKgm3;
result.air=struct('pressurePa',c.airPa,'densityKgm3',air.densityKgm3,'dynamicViscosityPas',air.dynamicViscosityPas, ...
    'kinematicViscosityM2s',nu,'conductivityWmK',air.conductivityWmK,'transportModel',air.transportModel, ...
    'model',air.model,'domain',air.domain,'temperatureRangeK',air.temperatureRangeK,'pressureCeilingPa',air.pressureCeilingPa, ...
    'sourceSHA256',air.sourceSHA256,'transportDOI','10.1023/B:IJOT.0000022327.04529.f3');
L=result.lengthM;result.reynoldsL=c.speedMps*L/nu;
result.grashofL=9.80665*(deltaT/result.filmTemperatureK)*L^3/nu^2;
if ~(result.reynoldsL>20 && result.reynoldsL<800 && result.grashofL>1e-4 && result.grashofL<1.4)
    result.reason='工况超出本文采用的强迫对流范围（20 < Re_L < 800，0.0001 < Gr_L < 1.4）。';return;
end

% 模块三：式（16）的形状项与式（28）的速度项分别保留；表列值另算敏感性。
e=result.eccentricity;
shape=sqrt(1-e^2)+asin(e)/e;
result.shapeNu0=16*e/(3*shape^2*log((1+e)/(1-e)));
result.nuL=result.shapeNu0+0.29*sqrt(result.reynoldsL);
result.meanHWm2K=air.conductivityWmK/L*result.nuL;
result.outwardW=result.meanHWm2K*result.areaM2*deltaT;
result.tableMinusFormulaW=(result.tableNu0-result.shapeNu0)*air.conductivityWmK/L*result.areaM2*deltaT;
result.available=true;
end
