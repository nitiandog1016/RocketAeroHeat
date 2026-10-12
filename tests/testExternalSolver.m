function testExternalSolver(outputFolder)
% 真实曲面箭体：零攻角鼻部切锥、柱段自由来流、攻角柱段压差。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);
[raw,~]=alignGeometryAxis(raw,'pca');
evalc('mesh=autoRemeshSurface(raw,cfg); mesh=fixNormalDirection(mesh,cfg);');
mesh=buildSurfaceCoordinates(mesh,cfg); mesh=analyzeLocalGeometry(mesh,cfg);
fs=calcFreestream(cfg.flight);
mesh=modifiedNewtonianPressure(mesh,fs,cfg);
mesh=findStagnationRegion(mesh,fs,cfg);
tic;
zero=solveExternalFlow(mesh,fs,cfg);
elapsed=toc;
e=zero.boundaryLayerEdge;
assert(zero.stagnationShock.available);
assert(strcmp(zero.stagnationShock.thermalPropertyMode,'thermallyPerfect'));
caloricShock=solveNormalShock(fs);
assert(zero.stagnationShock.postShockTemperature< ...
    caloricShock.postShockTemperature);
assert(zero.stagnationShock.postShockMach<1);
assert(isnan(zero.stagnationShock.shockStandOffDistance));
assert(zero.bowShockPreview.available);
assert(zero.bowShockPreview.standOffDistance>0);
assert(e.summary.validFaces>1000 && e.summary.coneFaces>0);
assert(all(e.thermalPropertyModeByFace(e.faceMethodID==2)==2));
assert(all(e.thermalPropertyModeByFace(e.faceMethodID==1)==2));
assert(e.summary.grazingFaces>1000);
assert(all(isfinite(e.shockAngleDeg(e.faceMethodID==2))));
assert(all(isnan(e.shockAngleDeg(e.faceMethodID~=2))));
assert(all(e.shockAngleDeg(e.faceMethodID==2)>asind(1/fs.Mach)));
assert(all(isnan(e.pressure(~e.algebraicallyValid))));
thermalFree=solveThermallyPerfectConicalShock(fs,0);
assert(all(e.totalPressure(e.algebraicallyValid)<= ...
    thermalFree.totalPressure*(1+1e-9)));
assert(any(e.pressure(e.faceMethodID==2)>fs.pressure));
body=zero.localGeometry.center(:,1)<-1 & zero.localGeometry.center(:,1)>-2;
assert(all(e.faceMethodID(body)==1));
zero=traceSurfaceStreamlines(zero,cfg);
zero=calcSurfaceHeating(zero,fs,cfg);
assert(isstruct(zero.heating.turbulentBoundaryLayer));
assert(zero.heating.turbulentBoundaryLayer.summary.integratedFaces>100);
assert(contains(zero.stagnationHeating.shockStateMethod,'Rankine-Hugoniot'));
assert(abs(zero.stagnationHeating.adiabaticWallTemperature- ...
    zero.stagnationShock.stagnationTemperature)<1e-9);
assert(zero.surfaceStreamlines.summary.connectedFaces>1000);
assert(zero.heating.summary.calculatedFaces>1000);
assert(zero.heating.summary.noseFaces>10);
assert(max(zero.heating.heatFlux(zero.heating.noseFace))<= ...
    zero.stagnationHeating.heatFlux*(1+1e-12));
assert(max(zero.heating.heatFlux(zero.heating.noseFace))> ...
    0.9*zero.stagnationHeating.heatFlux);
assert(isfield(zero,'boundaryLayer') && zero.boundaryLayer.summary.validFaces>1000);
material=struct('density',2700,'specificHeat',900,'thickness',0.003, ...
    'thermalConductivity',180,'emissivity',0.8, ...
    'radiationTemperature',cfg.flight.Twall);
heatedWall=calcWallThermalResponse(zero,[0 0.2 0.4],material,cfg);
w=heatedWall.wallThermal;
assert(w.wallTemperatureFeedback && w.summary.calculatedFaces>1000);
assert(any(abs(w.heatTransferCoefficientHistory( ...
    zero.heating.noseFace,2)-w.heatTransferCoefficientHistory( ...
    zero.heating.noseFace,1))>1e-8));
assert(w.stagnationPoint.available);
assert(abs(w.stagnationPoint.convectiveHeatFlux(1)- ...
    zero.stagnationHeating.heatFlux)<1e-9*zero.stagnationHeating.heatFlux);
assert(w.stagnationPoint.temperature(end)>cfg.flight.Twall);
assert(w.stagnationPoint.convectiveHeatFlux(end)< ...
    w.stagnationPoint.convectiveHeatFlux(1));
assert(~w.stagnationPoint.oneDimensional);
pointCapacity=material.density*material.specificHeat*material.thickness;
pointResidual=pointCapacity/0.2*( ...
    w.stagnationPoint.temperature(2)-w.stagnationPoint.temperature(1)) ...
    -w.stagnationPoint.netHeatFlux(2);
assert(abs(pointResidual)<1e-8*max(1,w.stagnationPoint.netHeatFlux(2)));
% 新默认耦合模型统一求一维导热；必须用厚度储热核对能量，
% 不能把表面温度乘整层热容当成材料储热。
ids=find(w.validFace);
assert(all(w.oneDimensionalFace(ids)) && ~any(w.lumpedFace));
assert(any(abs(w.heatTransferCoefficientHistory(ids,2)- ...
    w.heatTransferCoefficientHistory(ids,1))>1e-6));
stored=material.density*material.specificHeat*trapz(w.oneDimensional.depth, ...
    w.oneDimensional.finalDepthTemperature-cfg.flight.Twall,2);
energyResidual=stored-sum(w.netHeatFlux(ids,2:end).*diff(w.time),2);
assert(max(abs(energyResidual))<1e-3);
assert(all(isnan(w.temperature(~w.validFace,:)),'all'));
tilt=cfg.flight; tilt.alpha=10;
fsTilt=calcFreestream(tilt);
heatedWall=calcAerothermalExposure(heatedWall);
heatedWall.flightHeatingHistory=struct('time',[0 0.2 0.4]);
retrace=traceSurfaceStreamlines(heatedWall,cfg);
assert(~isfield(retrace,'boundaryLayer') && ...
    ~isfield(retrace,'heating') && ...
    ~isfield(retrace,'heatingCoverage') && ...
    ~isfield(retrace,'wallThermal') && ...
    ~isfield(retrace,'aerothermalExposure') && ...
    ~isfield(retrace,'flightHeatingHistory'));
reset=modifiedNewtonianPressure(heatedWall,fsTilt,cfg);
assert(~isfield(reset,'wallThermal') && ...
    ~isfield(reset,'aerothermalExposure') && ...
    ~isfield(reset,'flightHeatingHistory'));
mesh=modifiedNewtonianPressure(mesh,fsTilt,cfg);
mesh=findStagnationRegion(mesh,fsTilt,cfg);
angled=solveExternalFlow(mesh,fsTilt,cfg);
assert(angled.stagnationShock.available);
assert(~angled.bowShockPreview.available);
et=angled.boundaryLayerEdge;
assert(et.summary.detachedClosureFaces>0);
assert(et.summary.thermalWedgeFaces>0);
assert(all(et.thermalPropertyModeByFace(et.faceMethodID==3)==2));
assert(et.summary.thermalExpansionFaces>0);
assert(et.summary.thermalGrazingFaces>0);
thermalFree=solveThermallyPerfectPrandtlMeyerExpansion(fsTilt,0);
thermalGrazing=et.faceMethodID==1 & ...
    et.thermalPropertyModeByFace==2;
assert(max(abs(et.totalPressure(thermalGrazing)- ...
    thermalFree.totalPressure))<1e-8*thermalFree.totalPressure);
thermalExpansion=et.faceMethodID==4 & ...
    et.thermalPropertyModeByFace>=2;
assert(max(abs(et.totalPressure(thermalExpansion)- ...
    thermalFree.totalPressure))<1e-8*thermalFree.totalPressure);
assert(all(ismember(et.thermalPropertyModeByFace( ...
    et.faceMethodID==4 & ~et.thermalExpansionOutOfRangeFace),[2 3])));
assert(all(et.thermalPropertyModeByFace( ...
    et.faceMethodID==4 & et.temperature<150 & ...
    ~et.thermalExpansionOutOfRangeFace)==3));
assert(all(et.thermalPropertyModeByFace( ...
    et.faceMethodID==4 & et.thermalExpansionOutOfRangeFace)==1));
detached=et.faceMethodID==5;
assert(all(et.thermalPropertyModeByFace(detached)==2));
normal=solveThermallyPerfectNormalShock(fsTilt);
limit=solveThermallyPerfectObliqueShock(fsTilt,0);
limit=solveThermallyPerfectObliqueShock( ...
    fsTilt,limit.maximumAttachedDeflectionDeg);
assert(all(isnan(et.shockAngleDeg(detached))));
assert(all(et.totalPressure(detached)<= ...
    limit.totalPressure*(1+1e-9)));
assert(all(et.totalPressure(detached)>= ...
    normal.stagnationPressure*(1-1e-9)));
assert(all(et.pressure(detached)>= ...
    limit.pressure*(1-1e-9)));
assert(all(et.pressure(detached)<et.totalPressure(detached)));
assert(all(abs(et.pressure(detached)- ...
    et.density(detached).*fsTilt.R.*et.temperature(detached))<1e-8));
nearLimit=detached & et.incidenceAngleDeg< ...
    limit.maximumAttachedDeflectionDeg+1;
assert(any(nearLimit));
assert(max(abs(et.pressure(nearLimit)/limit.pressure-1))<0.03);
z=angled.localGeometry.normal(:,3);
windward=body & z<-0.5 & et.algebraicallyValid;
leeward=body & z>0.5 & et.algebraicallyValid;
assert(any(windward) && any(leeward));
assert(mean(et.pressure(windward))>fsTilt.pressure);
assert(mean(et.pressure(leeward))<fsTilt.pressure);
assert(et.summary.wedgeFaces>0 && et.summary.expansionFaces>0);
assert(all(isfinite(et.shockAngleDeg(et.faceMethodID==3))));
assert(all(isnan(et.shockAngleDeg(et.faceMethodID~=3))));
caloricCfg=cfg;
caloricCfg.heating.airThermalModel='caloricallyPerfect';
caloric=solveExternalFlow(mesh,fsTilt,caloricCfg);
assert(caloric.boundaryLayerEdge.summary.thermalWedgeFaces==0);
assert(caloric.boundaryLayerEdge.summary.thermalExpansionFaces==0);
assert(all(caloric.boundaryLayerEdge.thermalPropertyModeByFace<=1));
assert(contains(caloric.boundaryLayerEdge.wedgeThermalStatus,'定比热'));
angled=traceSurfaceStreamlines(angled,cfg);
angled=calcSurfaceHeating(angled,fsTilt,cfg);
assert(isstruct(angled.heating.turbulentBoundaryLayer));
assert(angled.heating.turbulentBoundaryLayer.summary.integratedFaces>100);
assert(angled.heating.summary.detachedClosureHeatFaces>0);
assert(all(angled.heating.faceMethodID( ...
    angled.heating.detachedClosureHeatFace)>=4));
coverage=angled.heatingCoverage;
assert(coverage.summary.computedAreaFraction>0.8);
% 攻角下热流页默认从迎风侧观察，避免把偏置热点藏在箭体背面。
viewCfg=cfg; viewCfg.flight=tilt; viewCfg.pressureGUI.visible='off';
viewFig=pressureGUI(angled,viewCfg);
viewUI=viewFig.UserData;
camera=viewUI.heatAxis.CameraPosition-viewUI.heatAxis.CameraTarget;
camera=camera/norm(camera);
windwardDirection=[0 -fsTilt.unitFlowDirection(2:3)];
windwardDirection=windwardDirection/norm(windwardDirection);
assert(dot(camera,windwardDirection)>0.9);
assert(strcmp(viewUI.heatDescription.Visible,'off'));
delete(viewFig);
assert(angled.continuumRegime.summary.directUnresolvedFaces>0);
assert(coverage.areaByReason(11)>0);
assert(~any(angled.heating.eligibleFace( ...
    angled.continuumRegime.unresolvedFace)));
assert(all(isnan(angled.heating.heatFlux( ...
    angled.continuumRegime.unresolvedFace))));
bridge=et.thermalExpansionOutOfRangeFace & ...
    ~angled.heating.noseFace;
if any(bridge)
    assert(~any(angled.heating.eligibleFace(bridge)));
    assert(any(angled.heating.thermalPropertyUnavailableFace(bridge)));
    assert(coverage.areaByReason(9)>0);
else
    assert(et.summary.lowTemperatureThermalFaces>0);
end
assert(sum(coverage.countByReason)==size(angled.F,1));
assert(isequal(coverage.reasonCode==0,angled.heating.eligibleFace));
assert(angled.heating.separationTrace.summary.onsetFaces==0);
% 角度节点连续插值后，分箱步长减半不应改变分离判断或热流覆盖。
quarterCfg=cfg; quarterCfg.externalFlow.angleStepDeg=0.25;
quarter=solveExternalFlow(mesh,fsTilt,quarterCfg);
quarter=traceSurfaceStreamlines(quarter,quarterCfg);
quarter=calcSurfaceHeating(quarter,fsTilt,quarterCfg);
assert(quarter.heating.separationTrace.summary.onsetFaces==0);
assert(abs(quarter.heatingCoverage.summary.computedAreaFraction- ...
    coverage.summary.computedAreaFraction)<1e-4);
samePressure=et.algebraicallyValid & ...
    quarter.boundaryLayerEdge.algebraicallyValid;
assert(max(abs(quarter.boundaryLayerEdge.pressure(samePressure)./ ...
    et.pressure(samePressure)-1))<0.02);
% 故意退回单步梯度以重现局部数值梯度敏感性，核对整条下游
% 路径不会被平板湍流公式或球鼻余弦公式悄悄补齐。
sensitiveCfg=cfg; sensitiveCfg.heating.gradientPathSteps=1;
sensitive=calcSurfaceHeating(angled,fsTilt,sensitiveCfg);
assert(sensitive.heating.separationTrace.summary.onsetFaces>0);
blocked=sensitive.heating.separationUnresolvedFace;
assert(all(isnan(sensitive.heating.heatFlux(blocked))));
assert(all(sensitive.heatingCoverage.flags.laminarSeparation(blocked)));
assert(all(sensitive.heatingCoverage.reasonCode(blocked)==6 | ...
    sensitive.heatingCoverage.reasonCode(blocked)==1 | ...
    sensitive.heatingCoverage.reasonCode(blocked)==2 | ...
    sensitive.heatingCoverage.reasonCode(blocked)==3 | ...
    sensitive.heatingCoverage.reasonCode(blocked)==4 | ...
    sensitive.heatingCoverage.reasonCode(blocked)==5));
assert(coverage.summary.computedAreaFraction> ...
    sensitive.heatingCoverage.summary.computedAreaFraction);
assert(angled.stagnationHeating.available);
assert(angled.heating.noseDistribution.available);
assert(angled.heating.summary.noseFaces>1000);
assert(angled.stagnationHeating.position(3)<-0.005);
noseIDs=find(angled.heating.noseFace);
[~,peakIndex]=max(angled.heating.heatFlux(noseIDs));
assert(angled.localGeometry.center(noseIDs(peakIndex),3)<-0.005);
assert(all(abs(angled.heating.heatFlux(noseIDs)- ...
    angled.stagnationHeating.heatFlux.* ...
    cosd(angled.heating.noseDistribution.polarAngleDeg(noseIDs)))<1e-8));
assert(angled.heating.summary.calculatedFaces>1000);
assert(any(angled.heating.eligibleFace(leeward)));
assert(contains(angled.heating.boundaryLayerMethod,'Thwaites'));
disabledCfg=cfg;
disabledCfg.externalFlow.detachedClosureEnabled=false;
unclosed=solveExternalFlow(mesh,fsTilt,disabledCfg);
unclosed=traceSurfaceStreamlines(unclosed,disabledCfg);
unclosed=calcSurfaceHeating(unclosed,fsTilt,disabledCfg);
assert(unclosed.boundaryLayerEdge.summary.detachedClosureFaces==0);
assert(coverage.summary.computedArea> ...
    unclosed.heatingCoverage.summary.computedArea);
% 不同 Mach 下，闭合压力随转角增加、总压随转角降低，且状态方程成立。
for testMach=[1.2 2 5 10]
    trialFlight=tilt; trialFlight.Mach=testMach;
    trialFree=calcFreestream(trialFlight);
    trial=modifiedNewtonianPressure(mesh,trialFree,cfg);
    trial=findStagnationRegion(trial,trialFree,cfg);
    trial=solveExternalFlow(trial,trialFree,cfg);
    trialEdge=trial.boundaryLayerEdge;
    if testMach==10
        assert(trialEdge.summary.thermalWedgeFaces==0);
        assert(contains(trialEdge.wedgeThermalStatus,'退回定比热'));
    end
    trialIDs=trialEdge.faceMethodID==5;
    assert(any(trialIDs));
    if all(trialEdge.thermalPropertyModeByFace(trialIDs)==2)
        trialNormal=solveThermallyPerfectNormalShock(trialFree);
        trialLimit=solveThermallyPerfectObliqueShock(trialFree,0);
        trialLimit=solveThermallyPerfectObliqueShock( ...
            trialFree,trialLimit.maximumAttachedDeflectionDeg);
    else
        trialNormal=solveNormalShock(trialFree);
        trialLimit=solveObliqueShock(trialFree,0);
        trialLimit=solveObliqueShock( ...
            trialFree,trialLimit.maximumAttachedDeflectionDeg);
    end
    assert(all(trialEdge.pressure(trialIDs)>= ...
        trialLimit.pressure*(1-1e-9)));
    assert(all(trialEdge.totalPressure(trialIDs)<= ...
        trialLimit.totalPressure*(1+1e-9)));
    assert(all(trialEdge.totalPressure(trialIDs)>= ...
        trialNormal.stagnationPressure*(1-1e-9)));
    assert(all(trialEdge.pressure(trialIDs)< ...
        trialEdge.totalPressure(trialIDs)));
    assert(all(abs(trialEdge.pressure(trialIDs)- ...
        trialEdge.density(trialIDs).*trialFree.R.* ...
        trialEdge.temperature(trialIDs))<1e-8));
    trialAngle=asind(trial.aerothermalRegion.mu(trialIDs));
    [~,angleOrder]=sort(trialAngle);
    trialP=trialEdge.pressure(trialIDs);
    trialP0=trialEdge.totalPressure(trialIDs);
    assert(all(diff(trialP(angleOrder))>=-1e-8));
    assert(all(diff(trialP0(angleOrder))<=1e-8));
    assert(max(trialP)>min(trialP));
end
% 攻角反号时偏置热点应翻到另一侧；超出适用角不计算驻点热流。
reverse=tilt; reverse.alpha=-10;
fsReverse=calcFreestream(reverse);
reverseMesh=findStagnationRegion(mesh,fsReverse,cfg);
reverseMesh=calcStagnationHeating(reverseMesh,fsReverse,cfg);
reverseNose=calcBluntNoseHeating(reverseMesh,fsReverse);
assert(reverseNose.available && reverseMesh.stagnationHeating.position(3)>0.005);
reverseIDs=find(reverseNose.faceMask);
[~,reversePeak]=max(reverseNose.heatFlux(reverseIDs));
assert(reverseMesh.localGeometry.center(reverseIDs(reversePeak),3)>0.005);
atLimit=tilt; atLimit.alpha=30;
fsLimit=calcFreestream(atLimit);
limitMesh=findStagnationRegion(mesh,fsLimit,cfg);
limitMesh=calcStagnationHeating(limitMesh,fsLimit,cfg);
limitNose=calcBluntNoseHeating(limitMesh,fsLimit);
assert(limitMesh.stagnation.noseRadiusEstimate.outerDepthRadiusFraction>= ...
    1+fsLimit.unitFlowDirection(1));
assert(limitMesh.stagnationHeating.available && limitNose.available);
tooSteep=tilt; tooSteep.alpha=35;
fsSteep=calcFreestream(tooSteep);
steep=findStagnationRegion(mesh,fsSteep,cfg);
steep=calcStagnationHeating(steep,fsSteep,cfg);
assert(~steep.stagnationHeating.available && isnan(steep.stagnationHeating.heatFlux));
if nargin>0
    % 同一 GUI 中核对偏置驻点的算法说明与鼻端热流云图。
    tiltedCfg=cfg; tiltedCfg.flight=tilt;
    tiltedCfg.pressureGUI.visible='off';
    fig=pressureGUI(angled,tiltedCfg);
    ui=fig.UserData;
    assert(any(contains(string(ui.heatDescription.Value),'入射角 10.0')));
    assert(any(contains(string(ui.heatDescription.Value),'偏置驻点')));
    assert(any(contains(string(ui.heatDescription.Value),'未假定自动再附')));
    ui.noseZoomButton.ButtonPushedFcn(ui.noseZoomButton,[]);
    assert(diff(xlim(ui.heatAxis))<4*angled.stagnationHeating.radius);
    ui.tabs.SelectedTab=ui.heatTab; drawnow;
    exportapp(fig,fullfile(outputFolder,'oblique_nose_heating_gui.png'));
    ui.heatChoice.Value=6;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
    drawnow;
    exportapp(fig,fullfile(outputFolder,'oblique_external_methods_gui.png'));
    ui.resetButton.ButtonPushedFcn(ui.resetButton,[]);
    ui.heatChoice.Value=12;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
    assert(isequal(ui.heatSurface.FaceVertexCData, ...
        double(ui.boundaryLayerEdge.thermalPropertyModeByFace)));
    drawnow;
    exportapp(fig,fullfile(outputFolder,'oblique_thermal_modes_gui.png'));
    ui.heatChoice.Value=14;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
    assert(isequal(ui.heatSurface.FaceVertexCData, ...
        double(ui.continuumRegime.applicabilityID)));
    drawnow;
    exportapp(fig,fullfile(outputFolder,'oblique_rarefaction_gui.png'));
    ui.heatChoice.Value=10;
    ui.heatChoice.ValueChangedFcn(ui.heatChoice,[]);
    drawnow;
    exportapp(fig,fullfile(outputFolder,'oblique_heat_methods_gui.png'));
    delete(fig);
    % 单步梯度只作敏感性演示：显示疑似分离及整条下游无解路径。
    separationCfg=tiltedCfg;
    separationCfg.heating.gradientPathSteps=1;
    separationFig=pressureGUI(sensitive,separationCfg);
    separationUI=separationFig.UserData;
    separationUI.heatChoice.Value=11;
    separationUI.heatChoice.ValueChangedFcn(separationUI.heatChoice,[]);
    separationUI.tabs.SelectedTab=separationUI.heatTab;
    view(separationUI.heatAxis,[0 0]);
    xlim(separationUI.heatAxis,[-0.5 0.06]);
    zlim(separationUI.heatAxis,[-0.15 0.15]);
    drawnow;
    assert(any(contains(string(separationUI.heatDescription.Value), ...
        '未假定自动再附')));
    assert(separationUI.heating.separationTrace.summary.onsetFaces>0);
    exportapp(separationFig, ...
        fullfile(outputFolder,'separation_coverage_gui.png'));
    delete(separationFig);
end
fprintf('Built-in external flow passed: zero-alpha valid %d/%d, cone %d, AoA=10 deg windward/leeward pressure %.1f/%.1f Pa; zero-alpha runtime %.2f s.\n', ...
    e.summary.validFaces,numel(e.pressure),e.summary.coneFaces, ...
    mean(et.pressure(windward)),mean(et.pressure(leeward)),elapsed);
fprintf('Stagnation point: Rn=%.5f m, q0=%.3f W/cm2, Tw(0.4 s)=%.2f K (test material).\n', ...
    zero.stagnationHeating.radius,zero.stagnationHeating.heatFlux/1e4, ...
    w.stagnationPoint.temperature(end));
fprintf('AoA=10 heat coverage: %.2f%% area; unresolved reason areas / total: %s.\n', ...
    100*coverage.summary.computedAreaFraction, ...
    mat2str(100*coverage.areaByReason(2:end)/coverage.summary.totalArea,3));
for item={zero,angled}
    h=item{1}.heating; b=h.turbulentBoundaryLayer;
    use=h.turbulentFace & b.validFace & isfinite(b.flatPlateSkinFrictionCoefficient);
    ratio=b.skinFrictionCoefficient(use)./b.flatPlateSkinFrictionCoefficient(use);
    assert(all(isfinite(ratio) & ratio>0));
    fprintf('Turbulent BL: handoff %d, integrated %d, local %d, fallback %d; Cf/flat median %.3f.\n', ...
        b.summary.laminarHandoffFaces,b.summary.integratedFaces, ...
        b.summary.localStartupFaces,b.summary.flatPlateFallbackFaces, ...
        median(ratio));
    fprintf('Head H range %.4f~%.4f; separation %d, numerical %d including downstream.\n', ...
        min(b.shapeFactor),max(b.shapeFactor), ...
        b.summary.separationUnresolvedFaces,b.summary.numericalUnresolvedFaces);
    transitionState=h.transition;
    fprintf('Transition continuity: carrier %d, inherited beyond heat gap %d, hidden onset %d.\n', ...
        transitionState.summary.turbulentCarrierFaces, ...
        transitionState.summary.inheritedThroughUnavailableFaces, ...
        transitionState.summary.onsetsOutsideHeatDomain);
    assert(all(transitionState.turbulentFace( ...
        transitionState.inheritedThroughUnavailableFace)));
    oldPartition=calcTransition(h.edgeReynolds, ...
        transitionState.resultEligibleFace,cfg,item{1}.surfaceStreamlines);
    fprintf('Transition classification changes versus heat-mask reset: %d faces.\n', ...
        nnz(oldPartition.turbulentFace~=transitionState.turbulentFace));
end
end
