function testRunAeroHeatSingleCase(out)
%TESTRUNAEROHEATSINGLECASE 校核显式单工况入口的真实计算与失败通知。
% 模块一：经典锥柱真实后端结果、当前非均匀/未知壁温及原输入不变。
% 模块二：改攻角使压力/载荷更新，已解快照不能继承旧工况或旧历史。
% 模块三：早期物理条件、旧几何、晚期模型失败保留原错误及真实阶段；
% 不发假完成，不返回中间结果。本测试不是实验精度验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(root,'outputs','v55_validation');end
if ~exist(out,'dir'),mkdir(out);end
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=5;
cfg.flight.beta=0;cfg.flight.Twall=300;
cfg.heating.airThermalModel='caloricallyPerfect';
mesh=buildCanonicalBody('coneCylinder',struct('resolution',24));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);original=mesh;events=struct([]);

%% 一、真实工况链与计算面积，给定量不计入求解面积。
[baseline,fs,report]=runAeroHeatSingleCase(mesh,cfg,@record);
assert(strcmp(report.status,'completed') && isequaln(mesh,original));
assert(isequaln(report.cfg,cfg) && report.faceCount==size(mesh.F,1));
assert(baseline.aerodynamicLoads.available);
assert(report.pressurePrescribedAreaFraction>0 && report.pressureComputedAreaFraction<1);
assert(report.heatComputedAreaFraction>0 && report.heatComputedAreaFraction<1);
assert(abs(report.pressureComputedAreaFraction+report.pressurePrescribedAreaFraction+ ...
    report.pressureMissingAreaFraction-1)<1e-14);
assert(isequaln(fs,baseline.pressure.freestream) && isequaln(events,report.events));
assert(numel(events)==22 && all(diff([events.fraction])>=0) && events(end).fraction==1);
assert(all(strcmp({events(1:2:end).status},'started')));
assert(all(strcmp({events(2:2:end).status},'completed')));

%% 二、当前面元壁温保留；一个未知面元不得用标量壁温补齐。
Tw=300+70*(mesh.localGeometry.center(:,2)>0);
known=find(baseline.heating.eligibleFace & Tw>300);
assert(numel(known)>2);unknown=known(1);Tw(unknown)=NaN;
nonuniform=cfg;nonuniform.heating.wallTemperatureByFace=Tw;
events=struct([]);
[current,~,currentReport]=runAeroHeatSingleCase(baseline,nonuniform,@record);
assert(isequaln(current.heating.wallTemperatureByFace,Tw));
assert(isnan(current.heating.heatFlux(unknown)) && ~current.heating.eligibleFace(unknown));
assert(isequaln(currentReport.cfg.heating.wallTemperatureByFace,Tw));
changed=known(2:end);changed=changed(isfinite(current.heating.heatFlux(changed)));
assert(~isempty(changed));
assert(any(abs(current.heating.heatFlux(changed)-baseline.heating.heatFlux(changed))>1e-3));
assert(isequaln(current.pressure.staticPressure,baseline.pressure.staticPressure));
assert(isequaln(baseline.heating.wallTemperatureByFace,repmat(300,size(mesh.F,1),1)));

%% 三、从旧快照再算新攻角，压力/总载荷同步，旧历史字段清除。
changedCfg=cfg;changedCfg.flight.alpha=15;seed=baseline;
seed.wallThermal=struct('stale',true);seed.flightHeatingHistory=struct('stale',true);
seed.aerothermalExposure=struct('stale',true);events=struct([]);
[changedCase,~,~]=runAeroHeatSingleCase(seed,changedCfg,@record);
assert(changedCase.pressure.freestream.alpha==15);
assert(any(abs(changedCase.pressure.staticPressure-baseline.pressure.staticPressure)>1));
assert(norm(changedCase.aerodynamicLoads.forceBody-baseline.aerodynamicLoads.forceBody)>1e-3);
assert(~any(isfield(changedCase,{'wallThermal','flightHeatingHistory','aerothermalExposure'})));

%% 四、受限底压在锥柱上被真实拒绝，链完成也不能宣称全表面求解。
restricted=cfg;restricted.basePressure=struct('model','sand93Laminar', ...
    'hasRearSupport',false,'hasGasInjection',false,'operationSource','软件验证指定孤立模型');
events=struct([]);
[partial,~,partialReport]=runAeroHeatSingleCase(mesh,restricted,@record);
base=mesh.localGeometry.isAftCap;
assert(any(base) && all(isnan(partial.pressure.staticPressure(base))));
assert(~partialReport.loadsAvailable && partialReport.pressureMissingAreaFraction>0);
assert(strcmp(partialReport.status,'completed') && events(end).fraction==1);
assert(all(isnan(partial.heating.heatFlux(base))));
assert(all(isnan(partial.boundaryLayerEdge.temperature(base))));

%% 五、真实失败阶段与错误；通知自身失败也不能遮盖后端错误。
bad=cfg;bad.flight.Mach=1;events=struct([]);
expectFailure(@()runAeroHeatSingleCase(mesh,bad,@record),'RocketAeroHeat:SingleCase:Mach','freestream');
stale=mesh;stale.V(1,2)=stale.V(1,2)+.001;events=struct([]);
expectFailure(@()runAeroHeatSingleCase(stale,cfg,@record), ...
    'RocketAeroHeat:Region:StaleGeometry','regions');
bad=cfg;bad.heating.noseModel='unsupported-test-model';events=struct([]);
expectFailure(@()runAeroHeatSingleCase(mesh,bad,@record), ...
    'RocketAeroHeat:Heating:NoseModel','heating');
events=struct([]);
expectFailure(@()runAeroHeatSingleCase(mesh,bad,@failedNotificationThrows), ...
    'RocketAeroHeat:Heating:NoseModel','heating');
bad=cfg;bad.heating.wallTemperatureByFace=[300;310];events=struct([]);
try
    runAeroHeatSingleCase(mesh,bad,@record);error('Test:UnexpectedSuccess','壁温长度错误没有被拒绝。');
catch ME
    assert(~strcmp(ME.identifier,'Test:UnexpectedSuccess'));
end
assert(strcmp(events(end).phase,'inputs') && strcmp(events(end).status,'failed'));
assert(~any(strcmp({events.status},'completed')));
save(fullfile(out,'single_case_entry.mat'),'baseline','current','changedCase','partial', ...
    'cfg','nonuniform','restricted','report','partialReport');
fprintf('Single-case entry passed: %d faces; real pressure/loads/heat, current nonuniform/unknown wall temperature, changed-condition refresh and early/late failure progress.\n',size(mesh.F,1));

    function record(event)
        if isempty(events),events=event;else,events(end+1,1)=event;end %#ok<AGROW>
    end
    function failedNotificationThrows(event)
        record(event);
        if strcmp(event.status,'failed')
            error('Test:NotificationFailure','故意破坏失败通知。');
        end
    end
    function expectFailure(action,id,phase)
        try
            action();error('Test:UnexpectedSuccess','应当失败的真实工况没有被拒绝。');
        catch ME
            assert(strcmp(ME.identifier,id),ME.identifier);
        end
        assert(strcmp(events(end).phase,phase) && strcmp(events(end).status,'failed'));
        assert(strcmp(events(end).errorIdentifier,id) && ~isempty(events(end).errorMessage));
        assert(events(end).fraction<1);
        assert(~any(strcmp({events.phase},'coverage') & strcmp({events.status},'completed')));
    end
end
