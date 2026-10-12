function testAxisymmetricShockExpansion(outputFolder)
% 真实箭体：零攻角直鼻锥之后，总压沿轴向不能恢复到自由来流。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);
[raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
fs=calcFreestream(cfg.flight);
mesh=modifiedNewtonianPressure(mesh,fs,cfg);
mesh=findStagnationRegion(mesh,fs,cfg);
local=solveExternalFlow(mesh,fs,cfg);
advanced=marchAxisymmetricShockExpansion(local,fs,cfg);
a=advanced.axisymmetricMarch;
e=advanced.boundaryLayerEdge;
assert(a.available && a.applied && a.summary.marchedFaces>10000);
assert(strcmp(a.thermalPropertyMode,'thermallyPerfect'));
assert(all(ismember(e.thermalPropertyModeByFace(a.faceMask),[2 3])));
assert(any(e.thermalPropertyModeByFace(a.faceMask)==3));
assert(min(e.temperature(a.faceMask))>=80);
assert(a.summary.expansionBands>=2 && a.summary.compressionBands>=1);
assert(a.summary.unresolvedBands==0);
assert(all(e.faceMethodID(a.faceMask)==6));
assert(all(e.algebraicallyValid(a.faceMask)));
assert(all(abs(e.pressure(a.faceMask)- ...
    e.density(a.faceMask).*fs.R.*e.temperature(a.faceMask))<1e-8));
assert(all(abs(e.totalTemperature(a.faceMask)- ...
    a.anchorTotalTemperature)<1e-5));
assert(all(e.totalPressure(a.faceMask)<= ...
    a.anchorTotalPressure*(1+1e-8)));
carried=flipud(a.bandTotalPressure(1:a.anchorBand-1));
assert(all(diff(carried)<=1e-8));
assert(all(isfinite(e.shockAngleDeg(e.axialWaveType==3))));
assert(all(isnan(e.shockAngleDeg( ...
    a.faceMask & e.axialWaveType~=3))));
x=advanced.localGeometry.center(:,1);
large=x<-.22 & x>-.35 & a.faceMask;
taper=x<-.45 & x>-.48 & a.faceMask;
small=x<-1 & x>-1.3 & a.faceMask;
assert(any(large) && any(taper) && any(small));
assert(mean(e.totalPressure(large))<0.5*fs.totalPressure);
assert(mean(e.pressure(large))<fs.pressure);
assert(mean(e.pressure(taper))<mean(e.pressure(large)));
assert(mean(e.totalPressure(small))<mean(e.totalPressure(large)));
assert(mean(e.pressure(small))>mean(e.pressure(taper)));
assert(mean(local.boundaryLayerEdge.totalPressure(large))> ...
    mean(e.totalPressure(large)));
% 人工温暖干空气状态：整个肩部膨胀仍处于热物性域，验证 NASA7
% 切锥锚点与后续膨胀/压缩确实作为同一条轴向流管运行。
warm=cfg; warm.flight.altitude=0; warm.flight.Mach=4;
air=standardAtmosphere(0);
air.temperature=400;
air.density=air.pressure/(air.R*air.temperature);
air.speedOfSound=sqrt(air.gamma*air.R*air.temperature);
air.dynamicViscosity=1.458e-6*air.temperature^(3/2)/ ...
    (air.temperature+110.4);
warm.flight.atmosphere=air;
warmStream=calcFreestream(warm.flight);
warmInput=modifiedNewtonianPressure(mesh,warmStream,warm);
warmInput=findStagnationRegion(warmInput,warmStream,warm);
warmLocal=solveExternalFlow(warmInput,warmStream,warm);
warmAdvanced=marchAxisymmetricShockExpansion(warmLocal,warmStream,warm);
warmMarch=warmAdvanced.axisymmetricMarch;
warmEdge=warmAdvanced.boundaryLayerEdge;
assert(warmMarch.applied && ...
    strcmp(warmMarch.thermalPropertyMode,'thermallyPerfect'));
assert(warmMarch.summary.unresolvedBands==0);
assert(all(warmEdge.thermalPropertyModeByFace(warmMarch.faceMask)==2));
assert(all(warmEdge.totalPressure(warmMarch.faceMask)<= ...
    warmMarch.anchorTotalPressure*(1+1e-8)));
assert(all(abs(warmEdge.pressure(warmMarch.faceMask)- ...
    warmEdge.density(warmMarch.faceMask).*warmStream.R.* ...
    warmEdge.temperature(warmMarch.faceMask))<1e-6));
carriedWarm=flipud(warmMarch.bandTotalPressure(1:warmMarch.anchorBand-1));
assert(all(diff(carriedWarm)<=1e-6));
% 人工冷来流跨越 80 K 下限时，整条轴向流管必须统一回退；
% 不能把低温定比热状态拼接在 NASA7 锥流状态的后面。
cold=warm;
cold.flight.atmosphere.temperature=160;
cold.flight.atmosphere.density= ...
    air.pressure/(air.R*cold.flight.atmosphere.temperature);
cold.flight.atmosphere.speedOfSound=sqrt( ...
    air.gamma*air.R*cold.flight.atmosphere.temperature);
cold.flight.atmosphere.dynamicViscosity= ...
    1.458e-6*160^(3/2)/(160+110.4);
coldStream=calcFreestream(cold.flight);
coldInput=modifiedNewtonianPressure(mesh,coldStream,cold);
coldInput=findStagnationRegion(coldInput,coldStream,cold);
coldLocal=solveExternalFlow(coldInput,coldStream,cold);
coldAdvanced=marchAxisymmetricShockExpansion(coldLocal,coldStream,cold);
coldMarch=coldAdvanced.axisymmetricMarch;
assert(coldMarch.applied && ...
    strcmp(coldMarch.thermalPropertyMode,'caloricallyPerfect'));
assert(strcmp(coldMarch.requestedThermalPropertyMode,'thermallyPerfect'));
assert(contains(coldMarch.thermalFallbackReason,'80 K'));
assert(all(coldAdvanced.boundaryLayerEdge.thermalPropertyModeByFace( ...
    coldMarch.faceMask)==1));
% 关闭模型或使用非零攻角时，原逐面元外流保持原值。
off=cfg; off.externalFlow.axisymmetricMarchEnabled=false;
unchanged=marchAxisymmetricShockExpansion(local,fs,off);
assert(~unchanged.axisymmetricMarch.applied);
assert(isequaln(unchanged.boundaryLayerEdge.pressure, ...
    local.boundaryLayerEdge.pressure));
tilt=cfg.flight; tilt.alpha=10;
tilted=marchAxisymmetricShockExpansion(local,calcFreestream(tilt),cfg);
assert(~tilted.axisymmetricMarch.applied);
advanced=traceSurfaceStreamlines(advanced,cfg);
advanced=calcSurfaceHeating(advanced,fs,cfg);
b=advanced.boundaryLayer; h=advanced.heating;
consistent=b.validFace & h.laminarFace;
assert(any(consistent));
assert(max(abs(b.recoveryTemperature(consistent)- ...
    h.recoveryTemperatureLaminar(consistent)))<1e-8);
expectedReference=advanced.boundaryLayerEdge.temperature(consistent)+ ...
    0.5*(cfg.flight.Twall- ...
    advanced.boundaryLayerEdge.temperature(consistent))+ ...
    0.22*(h.recoveryTemperatureLaminar(consistent)- ...
    advanced.boundaryLayerEdge.temperature(consistent));
assert(max(abs(b.referenceTemperature(consistent)-expectedReference))<1e-8);
assert(advanced.heatingCoverage.summary.computedAreaFraction>0.9);
assert(advanced.heatingCoverage.summary.computedAreaFraction>0.97);
assert(advanced.heating.summary.calculatedFaces>1000);
if nargin>0
    cfg.pressureGUI.visible='off';
    fig=pressureGUI(advanced,cfg);
    ui=fig.UserData;
    ui.heatChoice.Value=6;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
    ui.tabs.SelectedTab=ui.heatTab;
    view(ui.heatAxis,[0 0]);
    drawnow;
    assert(any(contains(string(ui.heatDescription.Value),'轴向波系')));
    assert(any(contains(string(ui.heatDescription.Value),'低温振动冻结')));
    assert(any(ui.heatSurface.FaceVertexCData==6));
    exportapp(fig,fullfile(outputFolder,'axisymmetric_march_gui.png'));
    ui.heatChoice.Value=12;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
    drawnow;
    assert(any(ui.heatSurface.FaceVertexCData==3));
    assert(isequal(ui.heatColorbar.Ticks,0:3));
    exportapp(fig,fullfile(outputFolder,'axisymmetric_thermal_modes_gui.png'));
    ui.tabs.SelectedTab=ui.mapTab;
    ui.selector.Value=5;
    ui.colorScale.Value=2;
    ui.selector.ValueChangedFcn(ui.selector,[]);
    view(ui.axes,[0 0]);
    drawnow;
    assert(isequaln(ui.surface.FaceVertexCData,e.pressure/1000));
    assert(strcmp(ui.colorScale.Enable,'on'));
    exportapp(fig,fullfile(outputFolder,'axisymmetric_pressure_gui.png'));
    delete(fig);
end
fprintf(['Axisymmetric march passed: anchor %.2f deg, %d faces, ' ...
    'p0 large/small cylinder %.3f/%.3f of freestream.\n'], ...
    a.anchorAngleDeg,a.summary.marchedFaces, ...
    mean(e.totalPressure(large))/fs.totalPressure, ...
    mean(e.totalPressure(small))/fs.totalPressure);
end
