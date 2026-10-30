function state=shockLayerStateFromHS(staticEnthalpy,entropy,gas)
%SHOCKLAYERSTATEFROMHS 从统一静焓及流线入口熵批量恢复层内状态。
% 模块一：单调焓括住温度，保护Newton只在声明温域中取值。
% 模块二：已收敛根立即冻结，防止精确根落在括号端点后被误换成中点。
% 模块三：同一熵恢复p/rho/声速；输出能量和熵回代误差，不改压力混合。
assert(isequal(size(staticEnthalpy),size(entropy))&&isreal(staticEnthalpy)&&isreal(entropy)&& ...
    all(isfinite(staticEnthalpy)&isfinite(entropy),'all'),'ShockPatch:ThermalInput','静焓和熵须同形且有限。');
lo=repmat(gas.temperatureRangeK(1),size(staticEnthalpy));hi=repmat(gas.temperatureRangeK(2),size(staticEnthalpy));
[~,hlo]=shockLayerThermalBatch(lo,gas);[~,hhi]=shockLayerThermalBatch(hi,gas);
assert(all(staticEnthalpy>=hlo&staticEnthalpy<=hhi,'all'), ...
    'ShockPatch:Temperature','静焓不能在声明温域内恢复，禁止夹温。');
if strcmp(gas.thermalPropertyMode,'caloricallyPerfect')
    T=staticEnthalpy/(gas.gamma*gas.R/(gas.gamma-1));iterations=0;
else
    T=lo+(hi-lo).*(staticEnthalpy-hlo)./(hhi-hlo);active=true(size(T));
    for iterations=1:40
        [cp,h]=shockLayerThermalBatch(T,gas);residual=h-staticEnthalpy;
        active=active&abs(residual)>16*eps*max(1,abs(staticEnthalpy));
        if ~any(active,'all'),break;end
        lower=active&residual<0;upper=active&~lower;lo(lower)=T(lower);hi(upper)=T(upper);
        candidate=T-residual./cp;bisect=active&(~isfinite(candidate)|candidate<=lo|candidate>=hi);
        candidate(bisect)=(lo(bisect)+hi(bisect))/2;T(active)=candidate(active);
    end
    assert(~any(active,'all'),'ShockPatch:ThermalRoot','批量焓根未收敛，不能返回候选温度。');
end
[cp,h,s0]=shockLayerThermalBatch(T,gas);
p=gas.referencePressurePa*exp((s0-entropy)/gas.R);rho=p./(gas.R*T);
assert(all(isfinite(p)&p>0&isfinite(rho)&rho>0,'all'),'ShockPatch:Pressure','同熵压力或密度失效。');
actualEntropy=s0-gas.R*log(p/gas.referencePressurePa);gamma=cp./(cp-gas.R);
state=struct('temperature',T,'pressure',p,'density',rho,'specificHeat',cp,'enthalpy',h, ...
    'entropy',actualEntropy,'sourceEntropy',entropy,'soundSpeedSquared',gamma.*gas.R.*T, ...
    'enthalpyResidual',h-staticEnthalpy,'entropyResidual',actualEntropy-entropy, ...
    'temperatureRootIterations',iterations);
end
