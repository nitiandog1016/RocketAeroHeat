function testLayeredWallRoutes()
%TESTLAYEREDWALLROUTES 验证分层在固定、逐站和耦合入口保持同一储热定义。
addpath(fileparts(fileparts(mfilename('fullpath'))));
layer=struct('density',1000,'specificHeat',800,'thickness',.004, ...
    'thermalConductivity',2,'thermalNodes',3);
inner=layer; inner.density=2500; inner.thermalConductivity=5;
T=[200;800;2000];
inner.propertyTable=struct('temperature',T,'specificHeat',600+.5*T, ...
    'thermalConductivity',3+.002*T);
mat=struct('layers',{{layer,inner}},'contactResistance',.002, ...
    'emissivity',.8,'radiationTemperature',300);
heat=struct('heatTransferCoefficient',[300;100;NaN], ...
    'recoveryTemperature',[1200;1000;NaN],'eligibleFace',[true;true;false]);
t=0:.1:1; cfg.flight.Twall=300;
m=struct('F',ones(3,3),'localGeometry',struct('area',ones(3,1)),'heating',heat);
fixed=calcWallThermalResponse(m,t,mat,cfg);
coupled=solveCoupledWallConduction(@(~)heat,heat,t,mat,300);
profile=struct('time',t,'initialWallTemperature',300, ...
    'heatTransferCoefficient',repmat(heat.heatTransferCoefficient,1,numel(t)), ...
    'recoveryTemperature',repmat(heat.recoveryTemperature,1,numel(t)), ...
    'eligibleFace',repmat(heat.eligibleFace,1,numel(t)));
history=calcWallThermalHistory(m,profile,mat,cfg);
assert(max(abs(fixed.wallThermal.temperature(1:2,:)-coupled.temperature(1:2,:)),[],'all')<1e-7);
assert(max(abs(history.wallThermal.temperature(1:2,:)-coupled.temperature(1:2,:)),[],'all')<1e-7);
for wall={fixed.wallThermal,coupled,history.wallThermal}
    w=wall{1}; m.wallThermal=w; e=calcAerothermalExposure(m);
    assert(e.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
    assert(all(isnan(w.Biot),'all') && ~any(w.lumpedFace));
    assert(isequal(size(w.oneDimensional.backTemperature),[2 numel(t)]));
    assert(all(isnan(w.temperature(3,:))));
    assert(max(abs(w.oneDimensional.interfaceFinal.energyResidual),[],'all')<1e-5);
    assert(contains(w.method,'接触热阻') && contains(w.materialMethod,'逐层'));
end
% 随壁温和时间变化的回调：检验真正的耦合残差及途中失效后的行映射。
w=solveCoupledWallConduction(@boundary,boundary(300*ones(3,1),0),t,mat,300);
m.wallThermal=w; e=calcAerothermalExposure(m);
assert(isequal(w.validFace,[true;false;false]));
assert(all(isnan(w.temperature(2:3,:)),'all'));
assert(size(w.oneDimensional.backTemperature,1)==1);
assert(size(w.oneDimensional.interfaceFinal.leftTemperature,1)==1);
assert(e.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
fprintf('Layered wall routes passed: fixed/local history/coupling, storage, back histories and dropped-face mapping.\n');
    function out=boundary(Tw,at)
        out=heat;
        out.heatTransferCoefficient(1:2)=[300;100].*(1+.0001*(Tw(1:2)-300))*(1+.01*at);
        if at>=.5, out.eligibleFace(2)=false; end
    end
end
