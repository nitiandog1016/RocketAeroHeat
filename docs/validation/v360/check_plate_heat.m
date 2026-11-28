function rows = check_plate_heat()
%CHECK_PLATE_HEAT 在公开加热平板的实测站位上复核工程换热公式的量纲和偏差。
% 输入固定为表 I(a) 的 127 ft/s 档。输出为铜带记录，热流采用正入壁。

%% 1. 按印刷表头读取质量通量、温差与热流，以国际表 Btu 换算 SI。
folder = fileparts(mfilename('fullpath'));
data = jsondecode(fileread(fullfile(folder,'table_ia_127.json')));
strip = data.rows(:,1);
assert(isequal(strip,(2:23)') && size(data.rows,2)==8);
massFlux = data.rows(:,2)*1000*0.45359237/(3600*0.3048^2);
density = data.freeDensityLbPerFt3*0.45359237/0.3048^3;
speed = massFlux/density;
edgeK = repmat((data.freeTemperatureF-32)*5/9+273.15,size(strip));
wallK = edgeK + data.rows(:,3)*5/9;
measuredQ = -data.rows(:,4)*1055.05585262/(3600*0.3048^2);
reEdge = data.rows(:,8)*1e6;

%% 2. 用项目实现求恢复温度、参考温度与等效比热。
% 采用干空气近似，气速来自各站质量通量与该速度档的自由流密度。
gas = struct('R',287.05287,'gamma',1.4);
pr = 0.72;
recoveryK = airRecoveryTemperature(edgeK,speed,pr^(1/3));
[referenceK,reference] = airReferenceTemperature(edgeK,wallK,recoveryK,gas,struct());
cpMean = airEnthalpySecantCp(recoveryK,wallK);
assert(all(reference.validFace & edgeK < recoveryK & recoveryK < wallK));

%% 3. 沿用局部平板分支的黏度比及湍流公式；实验在第一条铜带已触发转捩。
% Re_x 采用作者给出的实测列。对同一站，参考态与外缘态取相同静压。
viscosity = @(temp) 1.458e-6*temp.^1.5./(temp+110.4);
reRef = reEdge .* (edgeK./referenceK) .* (viscosity(edgeK)./viscosity(referenceK));
stanton = 0.0592 ./ reRef.^0.2 ./ (2*pr^(2/3));
estimatedQ = stanton .* massFlux .* cpMean .* (recoveryK-wallK);
errorPercent = 100*(estimatedQ./measuredQ-1);

%% 4. 保留局部输入与逐站有符号读数，用直接焓差独立检查热流方向和量纲。
[~,hRecovery] = airThermalProperties(recoveryK);
[~,hWall] = airThermalProperties(wallK);
directQ = stanton .* massFlux .* (hRecovery-hWall);
assert(all(isfinite(estimatedQ) & estimatedQ<0 & abs(estimatedQ-directQ)<1e-8));
rows = table(strip,speed,massFlux,edgeK,wallK,recoveryK,referenceK,reEdge,reRef, ...
    measuredQ,estimatedQ,errorPercent,'VariableNames',{'strip','speedMps','massFluxKgM2s', ...
    'freeK','givenWallK','recoveryK','referenceK','reX','referenceReX', ...
    'measuredQIntoWallWm2','estimatedQIntoWallWm2','relativeErrorPercent'});
end
