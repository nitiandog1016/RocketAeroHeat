function state=shockFittedState3D(value,gas,representation)
%SHOCKFITTEDSTATE3D 在明确的同一cp约定下转换三维五量保守态与原始态。
% 模块一：原始态为[rho ux uy uz T]，保守态为[rho rho*u rho*v rho*w rho*E]。
% 模块二：NASA7内能根严格位于声明温域，不使用常gamma反解或裁温。
% 模块三：输出实验室系总焓、压力、声速和熵；不赋予壁温/热流资格。
assert(isa(value,'double')&&isreal(value)&&size(value,2)==5&&size(value,1)>0&& ...
    all(isfinite(value),'all')&&all(value(:,1)>0),'ShockFit:State','需要有限且密度为正的N×5状态。');
rho=value(:,1);
if strcmp(representation,'primitive')
    velocity=value(:,2:4);T=value(:,5);[cp,h,s0]=shockLayerThermalBatch(T,gas);
elseif strcmp(representation,'conservative')
    velocity=value(:,2:4)./rho;internal=value(:,5)./rho-.5*sum(velocity.^2,2);
    [~,hlo]=shockLayerThermalBatch(gas.temperatureRangeK(1),gas);
    [~,hhi]=shockLayerThermalBatch(gas.temperatureRangeK(2),gas);
    elo=hlo-gas.R*gas.temperatureRangeK(1);ehi=hhi-gas.R*gas.temperatureRangeK(2);
    assert(all(internal>=elo&internal<=ehi),'ShockFit:Temperature','保守态内能超出声明温域；禁止裁温。');
    if strcmp(gas.thermalPropertyMode,'caloricallyPerfect')
        T=internal/(gas.R/(gas.gamma-1));
    else
        lo=repmat(gas.temperatureRangeK(1),size(rho));hi=repmat(gas.temperatureRangeK(2),size(rho));
        T=(lo+hi)/2;active=true(size(rho));
        for iteration=1:70
            [cp,h]=shockLayerThermalBatch(T,gas);f=h-gas.R*T-internal;
            active=active&abs(f)>8*eps*max(1,abs(internal));
            if ~any(active),break;end
            hi(active&f>0)=T(active&f>0);lo(active&f<=0)=T(active&f<=0);
            proposal=T-f./(cp-gas.R);bad=proposal<=lo|proposal>=hi|~isfinite(proposal);
            proposal(bad)=(lo(bad)+hi(bad))/2;T(active)=proposal(active);
        end
        assert(~any(active),'ShockFit:TemperatureRoot','同cp内能根未达到原数值门槛。');
    end
    [cp,h,s0]=shockLayerThermalBatch(T,gas);
else
    error('ShockFit:Representation','状态表示必须显式为primitive或conservative。');
end
p=rho.*gas.R.*T;E=h-gas.R*T+.5*sum(velocity.^2,2);Q=[rho rho.*velocity rho.*E];
state=struct('Q',Q,'density',rho,'velocity',velocity,'temperature',T,'pressure',p, ...
    'enthalpy',h,'totalEnthalpy',h+.5*sum(velocity.^2,2), ...
    'soundSpeed',sqrt(cp./(cp-gas.R).*gas.R.*T),'entropy',s0-gas.R*log(p/gas.referencePressurePa));
end
