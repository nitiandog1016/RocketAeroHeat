function testTimedCoupledWall()
%TESTTIMEDCOUPLEDWALL 时间相关热边界的独立积分和自适应时刻核对。
% 使用人工边界以隔离时间耦合，不把该边界当作箭体流动物理模型。
addpath(fileparts(fileparts(mfilename('fullpath'))));
mat=struct('density',1000,'specificHeat',1000,'thickness',0.01, ...
    'thermalConductivity',0.5,'emissivity',0.8,'radiationTemperature',300, ...
    'thermalNodes',6);
heat=boundary(300,0);
coarse=solveCoupledWallConduction(@boundary,heat,0:.5:1.5,mat,300);
fine=solveCoupledWallConduction(@boundary,heat,0:.1:1.5,mat,300);
[~,ref]=ode45(@derivative,[0 1.5],300*ones(6,1), ...
    odeset('RelTol',1e-10,'AbsTol',1e-10));
assert(abs(fine.temperature(end)-ref(end,1))< ...
    0.4*abs(coarse.temperature(end)-ref(end,1)));
stiffInitial=stiff(300,0);
a=solveCoupledWallConduction(@stiff,stiffInitial,[0 .1],mat,300);
assert(a.coupling.timeSubdivisions>0 && numel(a.time)>2);
for j=1:numel(a.time)
    h=stiff(a.temperature(j),a.time(j));
    assert(abs(a.convectiveHeatFlux(j)-h.heatTransferCoefficient* ...
        (h.recoveryTemperature-a.temperature(j)))<1e-7);
    assert(a.coupling.boundaryMetadata{j}.time==a.time(j));
end
m=struct('F',[1 2 3],'localGeometry',struct('area',1),'wallThermal',a);
m=calcAerothermalExposure(m);
assert(m.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-4);
fprintf('Timed coupling passed: independent ODE, time refinement, adaptive boundary times and energy.\n');

    function h=boundary(T,t)
        h=struct('heatTransferCoefficient',100*(1+.4*t+.001*(T-300)), ...
            'recoveryTemperature',1000+100*t,'eligibleFace',true, ...
            'boundaryMetadata',struct('time',t));
    end
    function h=stiff(T,t)
        h=boundary(T,t);
        h.heatTransferCoefficient=(100+1000/(1+exp((T-320)/2)))*(1+.1*t);
    end
    function d=derivative(t,Y)
        dx=mat.thickness/5; g=mat.thermalConductivity/dx;
        cap=mat.density*mat.specificHeat*dx*[.5;1;1;1;1;.5];
        h=boundary(Y(1),t);
        q=h.heatTransferCoefficient*(h.recoveryTemperature-Y(1))- ...
            mat.emissivity*5.670374419e-8*(Y(1)^4-mat.radiationTemperature^4);
        d=[q+g*(Y(2)-Y(1));g*(Y(1:4)-2*Y(2:5)+Y(3:6));g*(Y(5)-Y(6))]./cap;
    end
end
