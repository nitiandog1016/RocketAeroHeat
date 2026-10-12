function testVariableWallMaterial()
%TESTVARIABLEWALLMATERIAL 检查表内插值、精确积分、守恒和独立 ODE 基准。
% 以下人工线性 cp/k 仅用于验证，不作为任何真实材料推荐数据。
addpath(fileparts(fileparts(mfilename('fullpath'))));
mat=struct('density',1000,'specificHeat',750,'thickness',.01, ...
    'thermalConductivity',1.1,'emissivity',.8,'radiationTemperature',300,'thermalNodes',6);
T=[200;400;900;2000];
mat.propertyTable=struct('temperature',T,'specificHeat',600+.5*T, ...
    'thermalConductivity',.8+.001*T);
p=wallMaterialProperties(mat,[300 NaN 1000]);
assert(abs(p.specificHeat(1)-750)<1e-10 && ...
    abs(p.thermalConductivity(3)-1.8)<1e-12 && isnan(p.specificHeat(2)));
a=[250 800;400 900]; b=[850 300;400 1500];
expected=(b-a).*(600+.25*(a+b));
actual=integrateWallMaterialProperty(mat,'specificHeat',a,b);
assert(max(abs(actual-expected),[],'all')<1e-8);
% 交叉非平滑表节点时仍须积分，不能使用终点比热乘总温差。
kink=mat; kink.propertyTable.specificHeat=[700;900;850;1200];
assert(abs(integrateWallMaterialProperty(kink,'specificHeat',300,650)- ...
    (100*(800+900)/2+250*(900+875)/2))<1e-8);
for bad={199,2001,Inf,-1}
    failed=false; try, wallMaterialProperties(mat,bad{1}); catch, failed=true; end
    assert(failed);
end
bad=mat; bad.propertyTable.temperature(2)=200;
failed=false; try, wallMaterialProperties(bad,300); catch, failed=true; end
assert(failed);
% 求解所需温度超出表上限时应失败，不能把物性固定在端点外推。
narrow=mat; narrow.propertyTable=struct('temperature',[299;301], ...
    'specificHeat',[749.5;750.5],'thermalConductivity',[1.099;1.101]);
failed=false;
try
    calcWallConduction1D(1e5,1200,[0 1],narrow,300,[]);
catch problem
    failed=startsWith(problem.identifier,'RocketAeroHeat:MaterialConduction:');
end
assert(failed);
% 温度表为常数时与旧核一致，包含非均匀初始厚度温度。
flat=mat; flat.propertyTable.specificHeat(:)=750; flat.propertyTable.thermalConductivity(:)=1.1;
old=rmfield(flat,'propertyTable'); init=[300 310 320 320 310 300];
t=[0 .1 .5 1];
tabulated=calcWallConduction1D(100,1000,t,flat,init,[]);
constant=calcWallConduction1D(100,1000,t,old,init,[]);
assert(max(abs(tabulated.finalDepthTemperature-constant.finalDepthTemperature))<1e-6);
% 用解析 cp/k 直接写半离散 ODE；不调用生产物性积分或隐式求解器。
[~,reference]=ode45(@derivative,[0 1.5],300*ones(6,1), ...
    odeset('RelTol',1e-10,'AbsTol',1e-10));
coarse=calcWallConduction1D(300,1200,0:.5:1.5,mat,300,[]);
fine=calcWallConduction1D(300,1200,0:.1:1.5,mat,300,[]);
assert(abs(fine.surfaceTemperature(end)-reference(end,1))< ...
    .4*abs(coarse.surfaceTemperature(end)-reference(end,1)));
% 统一耦合核、局部固定壁温入口与给定逐站边界均不得退回常比热储热。
heat=struct('heatTransferCoefficient',300,'recoveryTemperature',1200,'eligibleFace',true);
w=solveCoupledWallConduction(@(~)heat,heat,0:.1:1.5,mat,300);
m=struct('F',[1 2 3],'localGeometry',struct('area',1),'wallThermal',w);
m=calcAerothermalExposure(m);
assert(m.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(max(abs(w.temperature-fine.surfaceTemperature))<1e-7);
cfg.flight.Twall=300;
m.heating=heat; m.heating.wallTemperature=300;
fixed=calcWallThermalResponse(m,0:.1:1.5,mat,cfg);
assert(fixed.wallThermal.oneDimensionalFace && ~fixed.wallThermal.lumpedFace);
assert(max(abs(fixed.wallThermal.temperature-w.temperature))<1e-7);
profile=struct('time',0:.1:1.5,'initialWallTemperature',300, ...
    'heatTransferCoefficient',300*ones(1,16),'recoveryTemperature',1200*ones(1,16), ...
    'eligibleFace',true(1,16));
history=calcWallThermalHistory(m,profile,mat,cfg);
assert(max(abs(history.wallThermal.temperature-w.temperature))<1e-7);
history=calcAerothermalExposure(history);
assert(history.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
fprintf('Variable wall material passed: interpolation/integrals, domain, constant limit, independent ODE and all thermal routes.\n');

    function d=derivative(~,Y)
        dx=mat.thickness/5; cp=600+.5*Y;
        capacity=mat.density*dx*[.5;1;1;1;1;.5].*cp;
        flux=(Y(2:end)-Y(1:end-1)).*(.8+.0005*(Y(2:end)+Y(1:end-1)))/dx;
        q=300*(1200-Y(1))-.8*5.670374419e-8*(Y(1)^4-300^4);
        d=[q+flux(1);flux(2:end)-flux(1:end-1);-flux(end)]./capacity;
    end
end
