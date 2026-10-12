function testBackWallBoundary()
%TESTBACKWALLBOUNDARY 双面换热的矩阵指数、稳态热阻、非线性 ODE 与收支。
% 人工材料/内环境仅用于软件验证，不是实际材料或舱内状态建议。
addpath(fileparts(fileparts(mfilename('fullpath'))));
a=struct('density',1000,'specificHeat',800,'thickness',.004, ...
    'thermalConductivity',2,'thermalNodes',3);
b=struct('density',1500,'specificHeat',1000,'thickness',.006, ...
    'thermalConductivity',4,'thermalNodes',3);
mat=struct('layers',{{a,b}},'contactResistance',.001, ...
    'emissivity',0,'radiationTemperature',300);
mat.backBoundary=struct('type','convectionRadiation','heatTransferCoefficient',80, ...
    'fluidTemperature',350,'emissivity',0,'radiationTemperature',340);
boundary=wallBackBoundary(mat);
[q,d]=evaluateBackWallHeatFlux([300 350 400 NaN],boundary);
assert(isequal(q(1:3),[-4000 0 4000]) && all(d(1:3)==80) && isnan(q(4)));
% 独立写出这个六节点两材料系统，不调用生产网格或通量函数。
mass=[800;1600;800;2250;4500;2250];
g=[1000 1000 1000 4/.003 4/.003]; K=zeros(6);
for j=1:5, K(j:j+1,j:j+1)=K(j:j+1,j:j+1)+g(j)*[1 -1;-1 1]; end
K(1,1)=K(1,1)+150; K(6,6)=K(6,6)+80;
forcing=[150*900;0;0;0;0;80*350]; equilibrium=K\forcing;
reference=equilibrium+expm(-diag(1./mass)*K*2)*(500-equilibrium);
coarse=calcWallConduction1D(150,900,0:.25:2,mat,500,[]);
fine=calcWallConduction1D(150,900,0:.05:2,mat,500,[]);
assert(norm(fine.finalDepthTemperature'-reference)<.3*norm(coarse.finalDepthTemperature'-reference));
input=sum(.05*(150*(900-fine.surfaceTemperature(2:end))-80*(fine.backTemperature(2:end)-350)));
assert(abs(sum(mass'.*(fine.finalDepthTemperature-500))-input)<1e-4);
% 长时间后应等于独立串联热阻的稳态解。
steady=calcWallConduction1D(150,900,0:200:20000,mat,500,[]);
flux=(900-350)/(1/150+.004/2+.001+.006/4+1/80);
assert(abs(steady.surfaceTemperature(end)-(900-flux/150))<1e-6);
assert(abs(steady.backTemperature(end)-(350+flux/80))<1e-6);
assert(abs(steady.interfaceFinal.contactHeatFlux-flux)<1e-4);
% 两侧辐射和逐层线性 cp/k 用独立 ODE 检验。
T=[200;1000;2000];
a.propertyTable=struct('temperature',T,'specificHeat',600+.4*T,'thermalConductivity',1+.002*T);
b.propertyTable=struct('temperature',T,'specificHeat',700+.3*T,'thermalConductivity',3+.001*T);
mat.layers={a,b}; mat.emissivity=.8; mat.backBoundary.emissivity=.7;
[~,ode]=ode45(@derivative,[0 2],500*ones(6,1),odeset('RelTol',1e-10,'AbsTol',1e-10));
coarse=calcWallConduction1D(150,900,0:.25:2,mat,500,[]);
fine=calcWallConduction1D(150,900,0:.05:2,mat,500,[]);
assert(norm(fine.finalDepthTemperature-ode(end,:))<.3*norm(coarse.finalDepthTemperature-ode(end,:)));
checkRoutes(mat);
% 常物性、变物性单层都要进入双面完整导热，低 Bi 也不能忽略内侧。
single=rmfield(a,'propertyTable'); single.emissivity=.8; single.radiationTemperature=300;
single.backBoundary=mat.backBoundary; checkRoutes(single);
single.propertyTable=a.propertyTable; checkRoutes(single);
direct=solveVariableMaterialConduction(150,900,0:.05:2,single,500,[]);
wrapped=calcWallConduction1D(150,900,0:.05:2,single,500,[]);
assert(isequal(direct.finalDepthTemperature,wrapped.finalDepthTemperature));
single.backBoundary.fluidTemperature=1200; single.backBoundary.radiationTemperature=1200;
reverse=calcWallConduction1D(0,300,[0 .1 1],single,300,[]);
assert(reverse.backTemperature(end)>reverse.surfaceTemperature(end));
assert(evaluateBackWallHeatFlux(reverse.backTemperature(end),wallBackBoundary(single))<0);
single.backBoundary=struct('type','adiabatic');
old=rmfield(single,'backBoundary');
x=calcWallConduction1D(150,900,[0 .1 1],single,500,[]);
y=calcWallConduction1D(150,900,[0 .1 1],old,500,[]);
assert(isequal(x.finalDepthTemperature,y.finalDepthTemperature));
zero=old; zero.backBoundary=mat.backBoundary;
zero.backBoundary.heatTransferCoefficient=0; zero.backBoundary.emissivity=0;
z=calcWallConduction1D(150,900,[0 .1 1],zero,500,[]);
assert(max(abs(z.finalDepthTemperature-y.finalDepthTemperature))<1e-7);
% 内侧边界也必须在耦合自适应子步使用新背温，不能只在输出时扣热量。
adaptive=rmfield(single,'propertyTable');
adaptive.backBoundary=mat.backBoundary;
first=stiffBoundary(500,0);
adaptiveWall=solveCoupledWallConduction(@stiffBoundary,first,[0 1],adaptive,500);
auditMesh=struct('F',[1 2 3],'localGeometry',struct('area',1),'wallThermal',adaptiveWall);
auditMesh=calcAerothermalExposure(auditMesh);
assert(auditMesh.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-3);
assert(max(abs(adaptiveWall.backHeatFlux-evaluateBackWallHeatFlux( ...
    adaptiveWall.backTemperature,wallBackBoundary(adaptive))),[],'all')<1e-9);
bad=mat; bad.backBoundary.type='guess'; mustFail(@()wallBackBoundary(bad));
bad=mat; bad.backBoundary.heatTransferCoefficient=-1; mustFail(@()wallBackBoundary(bad));
bad=mat; bad.backBoundary.radiationTemperature=0; mustFail(@()wallBackBoundary(bad));
fprintf('Back wall boundary passed: matrix exponential, steady resistances, nonlinear ODE, all routes, signed heat and adiabatic compatibility.\n');

    function d=derivative(~,Y)
        capacity=[1000*.002*[.5;1;.5].*(600+.4*Y(1:3)); ...
            1500*.003*[.5;1;.5].*(700+.3*Y(4:6))];
        links=[diff(Y(1:3)).*(1+.001*(Y(1:2)+Y(2:3)))/.002; ...
            (Y(4)-Y(3))/.001;diff(Y(4:6)).*(3+.0005*(Y(4:5)+Y(5:6)))/.003];
        front=150*(900-Y(1))-.8*5.670374419e-8*(Y(1)^4-300^4);
        back=80*(Y(6)-350)+.7*5.670374419e-8*(Y(6)^4-340^4);
        d=[front+links(1);diff(links);-links(end)-back]./capacity;
    end
    function out=stiffBoundary(Tw,at)
        out=struct('heatTransferCoefficient',(100+1000/(1+exp((Tw-520)/2)))*(1+.1*at), ...
            'recoveryTemperature',1000+100*at,'eligibleFace',true);
    end
end

function checkRoutes(mat)
% 完整背温与内侧热量必须在三个入口中一致，并参与总储热核对。
cfg.flight.Twall=500; t=0:.05:2;
heat=struct('heatTransferCoefficient',[150;NaN],'recoveryTemperature',[900;NaN], ...
    'eligibleFace',[true;false]);
m=struct('F',ones(2,3),'heating',heat,'localGeometry',struct('area',[2;3]));
fixed=calcWallThermalResponse(m,t,mat,cfg);
w=solveCoupledWallConduction(@(~)heat,heat,t,mat,500);
history=struct('time',t,'initialWallTemperature',500, ...
    'heatTransferCoefficient',repmat(heat.heatTransferCoefficient,1,numel(t)), ...
    'recoveryTemperature',repmat(heat.recoveryTemperature,1,numel(t)), ...
    'eligibleFace',repmat(heat.eligibleFace,1,numel(t)));
flight=calcWallThermalHistory(m,history,mat,cfg);
for value={w,fixed.wallThermal,flight.wallThermal}
    local=value{1}; assert(~any(local.lumpedFace));
    assert(max(abs(local.temperature(1,:)-w.temperature(1,:)))<1e-7);
    assert(max(abs(local.backTemperature(1,:)-w.backTemperature(1,:)))<1e-7);
    assert(all(isnan(local.backHeatFlux(2,:))));
    assert(contains(local.method,'背面 Newton') && ~contains(local.scope,'绝热背面'));
    m.wallThermal=local; e=calcAerothermalExposure(m); dose=e.aerothermalExposure;
    assert(dose.summary.maximumEnergyBalanceResidual<1e-3);
    assert(max(abs(dose.cumulativeNetHeat(1,:)-dose.cumulativeOuterNetHeat(1,:)+dose.cumulativeBackHeat(1,:)))<1e-7);
    assert(abs(dose.summary.totalBackEnergyResolved-2*dose.cumulativeBackHeat(1,end))<1e-8);
end
% 不允许非绝热结果丢失内侧热流，或把外侧净输入冒充材料储热。
m.wallThermal=rmfield(w,'backHeatFlux'); mustFail(@()calcAerothermalExposure(m));
m.wallThermal=w; m.wallThermal.netHeatFlux=w.outerNetHeatFlux;
mustFail(@()calcAerothermalExposure(m));
end

function mustFail(action)
failed=false; try, action(); catch, failed=true; end
assert(failed);
end
