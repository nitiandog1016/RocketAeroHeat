function item=createRocketAxisymmetricFluxCase(options)
%CREATEROCKETAXISYMMETRICFLUXCASE 完整箭体内部通量的显式独立研究工况。
% 模块一：继承真实几何、近壁层、压力源、远场/Tw/SA及均匀五量初猜。
% 模块二：仅新增内部Euler通量选择；开放面及既有物性/接受门槛不改。
% 模块三：新源码新链，原链不续写；不把稳健通量直接当成物理精度验证。
if nargin<1,options=struct;end
assert(isstruct(options)&&isscalar(options),'RocketFluxCase:Options','工况选项须为单一结构。');
base=options;policy='hlle-einfeldt';
if isfield(base,'interiorConvectiveFlux'),policy=base.interiorConvectiveFlux;base=rmfield(base,'interiorConvectiveFlux');end
assert(((ischar(policy)&&isrow(policy))||(isstring(policy)&&isscalar(policy)&&~ismissing(policy)))&& ...
    ismember(char(policy),{'hllc','hlle-davis','hlle-einfeldt'}), ...
    'RocketFluxCase:Policy','内部通量只接受hllc/hlle-davis/hlle-einfeldt。');policy=char(policy);
item=createRocketAxisymmetricPressureCase(base);item.boundary.interiorConvectiveFlux=policy;
item.caseFactory='createRocketAxisymmetricFluxCase';item.options=options;
item.inputContract.interiorConvectiveFluxPolicy=policy;
item.inputContract.interiorFluxIsPrescribedBoundary=false;
item.inputContract.interiorFluxIsPhysicalAccuracyValidation=false;
item.inputContract.changedFromV72=['same full uniform pressure-source case; explicit interior Euler flux=' policy];
end
