function testSphericalNoseWorkbench(out)
%TESTSPHERICALNOSEWORKBENCH 单工况球鼻方法在真实计算链路及统一GUI中的契约验证。
% 模块作用：校核同网格/工况切换、独立场量来源、拒绝区无回退和后体保留。
% 数值模型精度由独立球鼻剖面测试验证；本测试不等同接口或实验验证。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
if nargin<1,out=fullfile(root,'outputs','v53_validation');end
if ~exist(out,'dir'),mkdir(out);end
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=0;cfg.flight.beta=0;
cfg.flight.Twall=300;cfg.heating.sphericalNoseAngleStepDeg=5;
cfg.externalFlow.sphericalNoseMaxPolarAngleDeg=60;
[mesh,meta]=buildCanonicalBody('hemisphereCylinder',struct('resolution',24));
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
base=calcFlightStationState(mesh,cfg.flight,cfg);
assert(strcmp(base.heating.noseModel,'cosineExtension'));
fs=calcFreestream(cfg.flight);

%% 1. 实际GUI回调重算；计算与网页必须对应同一个工况/网格
fig=pressureGUI(base,cfg);cleanup=onCleanup(@()closeFigure(fig)); %#ok<NASGU>
selector=fig.UserData.noseSelector;
selector.Value='sphericalLocalSimilarity';selector.ValueChangedFcn(selector,[]);
drawnow;current=fig.UserData.currentMesh;selectedCfg=fig.UserData.currentCfg;
assert(strcmp(current.heating.noseModel,'sphericalLocalSimilarity'));
assert(strcmp(selectedCfg.heating.noseModel,'sphericalLocalSimilarity'));
assert(isequaln(base.V,current.V) && isequaln(base.F,current.F));
assert(isequaln(base.pressure.staticPressure,current.pressure.staticPressure));
nose=current.heating.noseDistribution;sphere=nose.edgeState.geometryFace(:);
computed=sphere & current.heating.eligibleFace(:);
assert(nnz(computed)>10,'正式GUI切换后必须实际得到非零球鼻区域解。');
high=sphere & nose.edgeState.algebraicallyValid & nose.edgeState.polarAngleDeg>=45;
assert(nnz(high)>10 && all(computed(high)), ...
    '实际GUI必须覆盖修正后45度以上的有效热边界，不能只验证低角度区域。');
assert(all(current.heating.faceMethodID(computed)==6));
assert(max(abs(current.heating.heatFlux(computed)-nose.heatFlux(computed)))<1e-9);
assert(max(abs(current.heating.recoveryTemperature(computed)-nose.recoveryTemperature(computed)))<1e-9);
assert(all(isnan(current.heating.heatFlux(sphere & ~computed))));
assert(~any(current.heating.faceMethodID(sphere)==3),'拒绝球鼻不得退回余弦补齐。');
assert(isequaln(base.heating.heatFlux(~sphere),current.heating.heatFlux(~sphere)), ...
    '独立鼻区选择不得改写后体热流。');
assert(~current.heating.historyCompatible);

%% 2. 同一区域模型分别记录热流、恢复温度及外缘压力/温度来源
data=buildAerothermalViewerData(current);ids=find(computed);
heat=field(data,'heatFlux');recovery=field(data,'recoveryTemperature');
edgeP=field(data,'edgePressure');edgeT=field(data,'edgeTemperature');
wall=field(data,'wallInput');solved=field(data,'wallSolved');
assert(max(abs(heat.values(computed)*1000-current.heating.heatFlux(computed)))<1e-8);
assert(all(wall.status==2) && all(solved.status==0) && all(isnan(solved.values)));
for id=ids(:)'
    assert(contains(heat.methodNames{heat.methodIndex(id)+1},'壁面导数'));
    assert(contains(recovery.methodNames{recovery.methodIndex(id)+1},'q(Twall)=0'));
    assert(contains(edgeP.scopeNames{edgeP.scopeIndex(id)+1},'接口尚未匹配'));
    assert(contains(edgeT.scopeNames{edgeT.scopeIndex(id)+1},'接口尚未匹配'));
end
keep=sphere & nose.edgeState.algebraicallyValid;
assert(max(abs(edgeP.values(keep)*1000-nose.edgeState.pressure(keep)))<1e-8);
ui=fig.UserData.webUI;start=tic;
while (~ui.html.UserData.ready || isempty(ui.html.UserData.rendered)) && toc(start)<40
    drawnow;pause(.1);
end
assert(ui.html.UserData.ready && ui.html.UserData.rendered.webglError==0);
assert(ui.html.UserData.rendered.visiblePixelSamples>0);
assert(strcmp(fig.UserData.noseSelector.Value,'sphericalLocalSimilarity'));
exportAerothermalWeb(current,fullfile(out,'spherical_nose'));
save(fullfile(out,'spherical_nose_workbench.mat'),'current','selectedCfg','base','cfg','fs','meta');

%% 3. 单时刻模型拒绝未验证的历史反馈；不会调用旧割线反馈冒充BVP耦合
expect(@()calcWallThermalResponse(current,[0 1],struct,selectedCfg), ...
    'RocketAeroHeat:WallThermal:SingleConditionModel');
expect(@()calcFlightHeatingHistory(current,struct,selectedCfg), ...
    'RocketAeroHeat:FlightHistory:SingleConditionModel');

%% 4. 非零攻角整链路无球鼻热流回退，并保存自身拒绝原因
obliqueCfg=selectedCfg;obliqueCfg.flight.alpha=10;
oblique=calcFlightStationState(mesh,obliqueCfg.flight,obliqueCfg);
blocked=oblique.heating.noseDistribution.edgeState.geometryFace(:);
assert(any(blocked));
assert(all(isnan(oblique.heating.heatFlux(blocked))));
assert(all(oblique.heating.faceMethodID(blocked)==0));
obliqueData=buildAerothermalViewerData(oblique);
assert(all(cellfun(@(s)contains(s,'零攻角')||contains(s,'轴对称'),obliqueData.missingReason(blocked))));

%% 5. 实际回调切回默认；默认模型与已保存生产算例一致
selector=fig.UserData.noseSelector;selector.Value='cosineExtension';
selector.ValueChangedFcn(selector,[]);drawnow;
assert(strcmp(fig.UserData.currentMesh.heating.noseModel,'cosineExtension'));
assert(isequaln(base.heating.heatFlux,fig.UserData.currentMesh.heating.heatFlux));
fprintf('Spherical nose workbench passed: %d/%d sphere faces computed; actual GUI callback, WebGL, per-field provenance, no fallback, unchanged afterbody and history guards.\n',nnz(computed),nnz(sphere));
end
function f=field(data,key)
f=data.fields(strcmp({data.fields.key},key));
end
function expect(fun,id)
try,fun();catch ex,assert(strcmp(ex.identifier,id),ex.identifier);return;end
error('Expected error not raised: %s',id);
end
function closeFigure(fig)
if isgraphics(fig),delete(fig);end
end
