function c=lewisCylinderCondition(staticPressurePa,airTemperatureK,wallTemperatureK)
%LEWISCYLINDERCONDITION 将新风洞设定换算为与 TM-83557 定义一致的速度与物性。
% 压力为工作段上游参考静压，温度为静温。三项输入均须由算例明确给定。
validateattributes(staticPressurePa,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(airTemperatureK,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(wallTemperatureK,{'numeric'},{'scalar','real','finite','positive'});
staticPressurePa=double(staticPressurePa);
airTemperatureK=double(airTemperatureK);
wallTemperatureK=double(wallTemperatureK);
assert(all([airTemperatureK wallTemperatureK]>=280 & [airTemperatureK wallTemperatureK]<=310), ...
    'RAH:Lewis:GasRange','本次近室温干空气物性核对覆盖 280–310 K。');
diameter=.066; re=177000; gasR=287.05; cp=1005; pr=.7;
density=staticPressurePa/(gasR*airTemperatureK);
% 同时满足质量通量雷诺数与试验的 Eckert 参考温度。46 m/s 为报告中的风洞上限。
mu=@(t) 1.716e-5*(t/273.15).^1.5.*(273.15+110.4)./(t+110.4);
tad=@(u) airTemperatureK+sqrt(pr)*u.^2/(2*cp);
tref=@(u) .5*wallTemperatureK+.28*airTemperatureK+.22*tad(u);
residual=@(u) density*u*diameter/mu(tref(u))-re;
assert(residual(46)>=0,'RAH:Lewis:FlowRange','当前输入达到目标雷诺数所需速度超过风洞的 46 m/s 上限。');
assert(staticPressurePa>=90000 && staticPressurePa<=110000, ...
    'RAH:Lewis:GasRange','本次近常压干空气物性核对覆盖 90–110 kPa。');
speed=fzero(residual,[0 46]);
assert(speed/sqrt(cp/(cp-gasR)*gasR*airTemperatureK)<.2, ...
    'RAH:Lewis:Mach','当前输入超出此低速定密度工况范围。');
c=struct('referencePressurePa',staticPressurePa,'airTemperatureK',airTemperatureK, ...
    'wallTemperatureK',wallTemperatureK,'velocityMps',speed,'densityKgM3',density, ...
    'adiabaticWallK',tad(speed),'referenceTemperatureK',tref(speed), ...
    'dynamicViscosityPaS',mu(tref(speed)),'conductivityWmK',mu(tref(speed))*cp/pr, ...
    'cpJkgK',cp,'prandtl',pr,'diameterM',diameter,'reynolds',re, ...
    'dynamicPressurePa',.5*density*speed^2);
end
