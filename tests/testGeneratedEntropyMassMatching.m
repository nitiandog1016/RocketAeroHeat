function report = testGeneratedEntropyMassMatching(outputFolder)
%TESTGENERATEDENTROPYMASSMATCHING 内部激波表到质量匹配核心的独立接口校核。
% 模块作用：实际调用球锥激波表生成器，由内部压力、度量及Gamma/Psi
% 结果检查质量匹配链路。本测试只按已知穿越质量坐标制造质量厚度；
% 这没有求解物理边界层，不证明真实鼻部/后体闭合或实验热流精度。
% 不造Gamma表，不输入外部CFD，不修改默认外流或网格。
projectRoot=fileparts(fileparts(mfilename('fullpath'))); addpath(projectRoot);
if nargin<1, outputFolder=fullfile(projectRoot,'outputs','v52_validation'); end
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
fs=struct('Mach',4,'gamma',1.4,'R',287.05,'pressure',5000, ...
    'temperature',220,'alpha',0,'beta',0);
fs.density=fs.pressure/(fs.R*fs.temperature);
fs.speed=fs.Mach*sqrt(fs.gamma*fs.R*fs.temperature);
geometry=struct('noseRadius',.1,'baseRadius',.5,'coneAngleDeg',25,'tipX',.8);
report=struct('validationKind','内部生成激波表/质量恢复接口与数值协变校核', ...
    'physicalBoundaryLayerClosed',false,'experimentallyValidated',false, ...
    'defaultSolverEnabled',false,'models',struct([]));

%% 1. 同一解析球锥在定比热、NASA7模式下真实生成激波质量表
for model={'caloricallyPerfect','thermallyPerfect'}
    mode=model{1};
    cfg=struct('heating',struct('airThermalModel',mode), ...
        'externalFlow',struct('entropyShock',struct('stationCount',41, ...
        'pressureTableNodes',65,'relativeTolerance',1e-9,'absoluteTolerance',1e-12)), ...
        'entropyLayer',struct('scanSubdivisions',4,'relativeTolerance',1e-9));
    generated=solveAxisymmetricEntropyShockTable(fs,geometry,cfg);
    assert(generated.available,generated.status);
    assert(all(diff(generated.massCoordinate)>0) && all(diff(generated.shockCrossingAngleDeg)<=0));
    assert(generated.massCoordinate(1)==0 && generated.shockCrossingAngleDeg(1)==90);
    before=generated;
    % 目标质量坐标来自实际生成表；线性插值与核心声明的表合同一致。
    knownPsi=generated.massCoordinate(end)*[.2 .5 .8];
    knownGamma=interp1(generated.massCoordinate,generated.shockCrossingAngleDeg,knownPsi,'linear');
    stationIndex=find(generated.surfaceRadius>=.49,1,'first');
    assert(~isempty(stationIndex));
    physicalPressure=generated.physicalPressure(stationIndex);
    physicalMetric=generated.streamlineMetric(stationIndex);
    knownEdge=recoverShockEntropyEdgeState(fs,knownGamma,physicalPressure,cfg);
    assert(all(knownEdge.algebraicallyValid));
    manufacturedThickness=knownPsi./(knownEdge.density.*knownEdge.speed*physicalMetric);
    station=struct('physicalPressure',physicalPressure,'streamlineMetric',physicalMetric, ...
        'massThickness',manufacturedThickness,'separated',false, ...
        'physicalPressureSource',generated.physicalPressureSource, ...
        'streamlineMetricSource',generated.streamlineMetricSource, ...
        'massThicknessSource','接口校核用已知根制造质量厚度；未解物理边界层', ...
        'separationSource','接口校核明确假设附着；未进行物理分离判定');

    %% 2. 直接消费完整生成表，检查根、实际输入及逐量来源
    matched=recoverEntropyLayerEdgeState(fs,station,generated,cfg);
    assert(isequaln(before,generated),'消费生成表不能改写其物理或数值数据');
    assert(all(matched.algebraicallyValid),matched.status);
    relativePsi=abs(matched.massCoordinate./knownPsi-1);
    assert(max(relativePsi)<2e-9);
    assert(max(abs(matched.shockCrossingAngleDeg-knownGamma))<1e-7);
    assert(max(abs(matched.speed./knownEdge.speed-1))<1e-8);
    assert(max(abs(matched.temperature./knownEdge.temperature-1))<1e-8);
    assert(all(matched.pressure==physicalPressure));
    assert(all(matched.streamlineMetricInput==physicalMetric));
    assert(strcmp(matched.massCoordinateConvention,generated.massCoordinateUnit));
    assert(matched.inputTableAvailable);
    assert(strcmp(matched.inputMetadata.physicalPressure.source{1},generated.physicalPressureSource));
    assert(strcmp(matched.inputMetadata.streamlineMetric.source{1},generated.streamlineMetricSource));
    assert(contains(matched.inputMetadata.massThickness.source{1},'未解物理边界层'));
    assert(strcmp(matched.fieldMetadata.pressure.valueKind,'input'));
    assert(strcmp(matched.fieldMetadata.temperature.valueKind,'computed'));
    assert(contains(matched.massCoordinateSource,'内部ODE'));
    assert(~matched.geometryVerified && ~matched.defaultFlowIntegrated);

    %% 3. 几何全尺度协变：r/S/h/delta乘k，真实生成Psi应乘k^2
    factor=1.7;
    largerGeometry=geometry;
    largerGeometry.noseRadius=geometry.noseRadius*factor;
    largerGeometry.baseRadius=geometry.baseRadius*factor;
    largerGeometry.tipX=geometry.tipX*factor;
    larger=solveAxisymmetricEntropyShockTable(fs,largerGeometry,cfg);
    assert(larger.available,larger.status);
    assert(numel(larger.massCoordinate)==numel(generated.massCoordinate));
    assert(max(abs(larger.massCoordinate/factor^2-generated.massCoordinate))<1e-10);
    assert(max(abs(larger.shockCrossingAngleDeg-generated.shockCrossingAngleDeg))<1e-8);
    scaledStation=station;
    scaledStation.streamlineMetric=larger.streamlineMetric(stationIndex);
    scaledStation.physicalPressure=larger.physicalPressure(stationIndex);
    scaledStation.massThickness=station.massThickness*factor;
    scaled=recoverEntropyLayerEdgeState(fs,scaledStation,larger,cfg);
    assert(all(scaled.algebraicallyValid),scaled.status);
    assert(max(abs(scaled.massCoordinate/factor^2-matched.massCoordinate))<1e-9);
    assert(max(abs(scaled.speed-matched.speed))<1e-7);
    assert(max(abs(scaled.temperature-matched.temperature))<1e-7);
    assert(max(abs(scaled.totalPressure./matched.totalPressure-1))<1e-8);

    %% 4. 原表的质量域外、不可用状态、变攻角工况与物性错误不能误消费
    outside=station; outside.massThickness=1e6*station.massThickness;
    refused=recoverEntropyLayerEdgeState(fs,outside,generated,cfg);
    assert(all(refused.statusCode==6)); assertPhysicalMissing(refused);
    badTable=generated; badTable.available=false;
    refused=recoverEntropyLayerEdgeState(fs,station,badTable,cfg);
    assert(all(refused.statusCode==3)); assertPhysicalMissing(refused);
    angled=fs; angled.alpha=1;
    refused=recoverEntropyLayerEdgeState(angled,station,generated,cfg);
    assert(all(refused.statusCode==3)); assertPhysicalMissing(refused);
    otherModel=cfg;
    if strcmp(mode,'caloricallyPerfect'), otherModel.heating.airThermalModel='thermallyPerfect';
    else, otherModel.heating.airThermalModel='caloricallyPerfect'; end
    refused=recoverEntropyLayerEdgeState(fs,station,generated,otherModel);
    assert(all(refused.statusCode==3)); assertPhysicalMissing(refused);
    tooHigh=station; tooHigh.physicalPressure=1e99;
    refused=recoverEntropyLayerEdgeState(fs,tooHigh,generated,cfg);
    assert(all(refused.statusCode==4)); assertPhysicalMissing(refused);

    entry=struct('thermalPropertyMode',mode,'table',generated,'knownMassCoordinate',knownPsi, ...
        'manufacturedMassThickness',manufacturedThickness,'station',station,'matched',matched, ...
        'geometryScaleFactor',factor,'maximumRelativeMassError',max(relativePsi), ...
        'maximumScaledMassError',max(abs(scaled.massCoordinate/factor^2-matched.massCoordinate)), ...
        'physicalBoundaryLayerClosed',false);
    if isempty(report.models), report.models=entry; else, report.models(end+1)=entry; end
    fprintf('Generated entropy-mass interface %s: %d real ODE table nodes, %d scan samples, mass relative error %.3g, scaled mass error %.3g kg/s.\n', ...
        mode,numel(generated.massCoordinate),matched.scanSampleCount,max(relativePsi),entry.maximumScaledMassError);
end
save(fullfile(outputFolder,'generated_entropy_mass_matching.mat'),'report','-v7');
fprintf('Generated entropy mass matching passed: actual internal Gamma/Psi tables, input provenance, caloric/NASA7 roots, geometry-scale covariance and domain/state refusals.\n');
fprintf('Manufactured mass thickness tests interface closure only; physical boundary-layer mass thickness, heat flux, nose/afterbody validation and experimental comparison remain incomplete.\n');
end

function assertPhysicalMissing(edge)
assert(~edge.available && ~any(edge.algebraicallyValid(:)) && ~any(edge.massMatchingSolved(:)));
for name={'pressure','temperature','density','speed','totalPressure','massCoordinate','shockCrossingAngleDeg'}
    assert(all(isnan(edge.(name{1})(:))));
end
end
