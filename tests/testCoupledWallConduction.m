function testCoupledWallConduction()
%TESTCOUPLEDWALLCONDUCTION 验证耦合求解的独立基准、非局部响应和储热。
% 先与整段一维导热对照，确认每步接续完整剖面且不重复推进时间。
% 再加入温度依赖与跨面元热边界，检查真正耦合残差及失效后留空。
addpath(fileparts(fileparts(mfilename('fullpath'))));
mat=struct('density',1000,'specificHeat',1000,'thickness',0.01, ...
    'thermalConductivity',0.5,'emissivity',0.8,'radiationTemperature',300);
t=[0 0.2 0.7 1.5];
heat=struct('heatTransferCoefficient',[100;180], ...
    'recoveryTemperature',[1000;1100],'eligibleFace',[true;true]);
fixed=solveCoupledWallConduction(@(~)heat,heat,t,mat,300);
reference=calcWallConduction1D(heat.heatTransferCoefficient, ...
    heat.recoveryTemperature,t,mat,300,[]);
assert(max(abs(fixed.temperature-reference.surfaceTemperature),[],'all')<1e-6);
assert(max(abs(fixed.oneDimensional.finalDepthTemperature- ...
    reference.finalDepthTemperature),[],'all')<1e-6);
coupled=solveCoupledWallConduction(@feedback,heat,t,mat,300);
assert(coupled.temperature(2,end)>fixed.temperature(2,end));
assert(any(coupled.coupling.iterationsByStep(2:end)>2));
for j=2:numel(t)
    actual=feedback(coupled.temperature(:,j));
    assert(max(abs(coupled.convectiveHeatFlux(:,j)- ...
        actual.heatTransferCoefficient.* ...
        (actual.recoveryTemperature-coupled.temperature(:,j))))<1e-8);
end
mesh.F=repmat([1 2 3],2,1); mesh.localGeometry.area=[1;2];
mesh.wallThermal=coupled;
mesh=calcAerothermalExposure(mesh);
assert(mesh.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-4);
assert(all(coupled.coupling.maximumFluxResidualByStep(2:end)<1e-4));
% 独立 ode45 解同一厚度半离散能量方程；缩小隐式时间步应趋近它。
% 此参考不调用壁温步进器，避免仅用同一个实现与自身互相比较。
mat.thermalNodes=6;
[~,referenceODE]=ode45(@thermalDerivative,[0 1.5],300*ones(12,1), ...
    odeset('RelTol',1e-10,'AbsTol',1e-10));
target=reshape(referenceODE(end,:),2,6);
coarse=solveCoupledWallConduction(@feedback,heat,0:0.5:1.5,mat,300);
fine=solveCoupledWallConduction(@feedback,heat,0:0.1:1.5,mat,300);
coarseError=max(abs(coarse.temperature(:,end)-target(:,1)));
fineError=max(abs(fine.temperature(:,end)-target(:,1)));
assert(fineError<0.4*coarseError);
mat=rmfield(mat,'thermalNodes');
% 回调返回行向量时也不得生成隐式广播的 N×N 数组。
row=heat; row.heatTransferCoefficient=row.heatTransferCoefficient';
row.recoveryTemperature=row.recoveryTemperature'; row.eligibleFace=row.eligibleFace';
rowResult=solveCoupledWallConduction(@(~)row,row,t,mat,300);
assert(isequaln(rowResult.temperature,fixed.temperature));
% 第一个面元在升温后失效；第二面元完全独立，必须继续求解。
lost=solveCoupledWallConduction(@failWhenHot,heat,t,mat,300);
assert(~lost.validFace(1) && lost.validFace(2));
assert(all(isnan(lost.temperature(1,:))));
assert(lost.coupling.firstFailureTime(1)>0);
assert(all(isfinite(lost.temperature(2,:))));
% 构造连续但陡峭的换热响应，强制时间细分；保留请求端点并用
% 实际内部站积分热量，不能把多个子步当成一个右端热流矩形。
stiffInitial=stiffBoundary([300;300]);
adaptive=solveCoupledWallConduction(@stiffBoundary,stiffInitial,[0 0.1],mat,300);
assert(adaptive.coupling.timeSubdivisions>0 && numel(adaptive.time)>2);
assert(adaptive.time(1)==0 && adaptive.time(end)==0.1);
assert(all(diff(adaptive.time)>0));
mesh.wallThermal=adaptive; mesh=calcAerothermalExposure(mesh);
assert(mesh.aerothermalExposure.summary.maximumEnergyBalanceResidual<1e-4);
fprintf('Coupled wall passed: profile continuation, nonlinear/nonlocal boundary, energy, row inputs and loss masking.\n');

    function h=feedback(T)
        h=heat;
        h.heatTransferCoefficient(1)=100*(1+0.001*(T(1)-300));
        h.heatTransferCoefficient(2)=180*(1+0.003*(T(1)-300)+0.001*(T(2)-300));
    end
    function h=failWhenHot(T)
        h=heat;
        h.eligibleFace(1)=isfinite(T(1)) && T(1)<301;
    end
    function h=stiffBoundary(T)
        h=heat;
        h.heatTransferCoefficient=100+1000./(1+exp((T-320)/2));
    end
    function dY=thermalDerivative(~,Y)
        P=reshape(Y,2,6); dx=mat.thickness/5;
        g=mat.thermalConductivity/dx;
        cap=mat.density*mat.specificHeat*dx*[0.5 1 1 1 1 0.5];
        h=feedback(P(:,1));
        q=h.heatTransferCoefficient.*(h.recoveryTemperature-P(:,1))- ...
            mat.emissivity*5.670374419e-8*(P(:,1).^4-mat.radiationTemperature^4);
        change=zeros(2,6);
        change(:,1)=q+g*(P(:,2)-P(:,1));
        change(:,2:5)=g*(P(:,1:4)-2*P(:,2:5)+P(:,3:6));
        change(:,6)=g*(P(:,5)-P(:,6));
        dY=reshape(change./cap,[],1);
    end
end
