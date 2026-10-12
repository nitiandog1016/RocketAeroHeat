function testLayeredWallConduction()
%TESTLAYEREDWALLCONDUCTION 独立验证分层热容、接触热阻及变物性导热。
% 人工材料只用于验证；矩阵指数/ODE 基准不调用生产通量或积分函数。
addpath(fileparts(fileparts(mfilename('fullpath'))));
a=struct('thickness',.004,'density',1200,'specificHeat',800, ...
    'thermalConductivity',2,'thermalNodes',3);
b=a; b.density=2700; b.specificHeat=900; b.thermalConductivity=8;
mat=struct('layers',{{a,a}},'contactResistance',0, ...
    'emissivity',0,'radiationTemperature',300);
single=a; single.thickness=.008; single.thermalNodes=5;
single.emissivity=0; single.radiationTemperature=300;
t=0:.1:2;
split=calcWallConduction1D([200;400],[1000;1200],t,mat,300,[]);
whole=calcWallConduction1D([200;400],[1000;1200],t,single,300,[]);
assert(max(abs(split.finalDepthTemperature-whole.finalDepthTemperature),[],'all')<1e-7);
mat.layers={a,b};
for rc=[0 .002]
    mat.contactResistance=rc;
    [mass,g]=independentGrid(a,b,rc);
    n=numel(mass); K=zeros(n);
    for j=1:n-1, K(j:j+1,j:j+1)=K(j:j+1,j:j+1)+g(j)*[1 -1;-1 1]; end
    K(1,1)=K(1,1)+400;
    ref=1200+expm(-diag(1./mass)*K*2)*(300-1200)*ones(n,1);
    coarse=calcWallConduction1D(400,1200,0:.2:2,mat,300,[]);
    fine=calcWallConduction1D(400,1200,0:.05:2,mat,300,[]);
    assert(norm(fine.finalDepthTemperature'-ref)<.4*norm(coarse.finalDepthTemperature'-ref));
    energy=sum(mass'.*(fine.finalDepthTemperature-300));
    input=sum(400*(1200-fine.surfaceTemperature(2:end))*.05);
    assert(abs(energy-input)<1e-4);
    audit=fine.interfaceFinal;
    assert(max(abs(audit.energyResidual),[],'all')<1e-5);
    if rc>0
        assert(audit.leftTemperature>audit.rightTemperature);
        assert(abs(audit.contactHeatFlux-audit.heatFluxFromLeft)<1e-5);
        assert(nnz(diff(fine.depth)==0)==1);
    else
        assert(audit.leftTemperature==audit.rightTemperature && all(diff(fine.depth)>0));
    end
end
% 每层独立线性 cp/k，基准 ODE 显式计算两种材料及接触热阻。
T=[200;800;2000];
a.propertyTable=struct('temperature',T,'specificHeat',500+.5*T,'thermalConductivity',1+.001*T);
b.propertyTable=struct('temperature',T,'specificHeat',700+.2*T,'thermalConductivity',4+.002*T);
mat.layers={a,b}; mat.contactResistance=.002;
[~,ref]=ode45(@derivative,[0 2],300*ones(6,1),odeset('RelTol',1e-10,'AbsTol',1e-10));
coarse=calcWallConduction1D(400,1200,0:.2:2,mat,300,[]);
fine=calcWallConduction1D(400,1200,0:.05:2,mat,300,[]);
assert(norm(fine.finalDepthTemperature-ref(end,:))<.4*norm(coarse.finalDepthTemperature-ref(end,:)));
assert(max(abs(fine.interfaceFinal.energyResidual),[],'all')<1e-5);
bad=mat; bad.contactResistance=[0 0]; mustFail(@()buildLayeredWallGrid(bad));
bad=mat; bad.contactResistance=-1; mustFail(@()buildLayeredWallGrid(bad));
bad=mat; bad.contactResistance='a'; mustFail(@()buildLayeredWallGrid(bad));
bad=mat; bad.contactResistance=1e-320;
mustFail(@()calcWallConduction1D(400,1200,[0 1],bad,300,[]));
bad=mat; bad.propertyTable=a.propertyTable; mustFail(@()normalizeLayeredWallMaterial(bad,300));
mustFail(@()calcWallConduction1D(400,1200,[0 1],mat,199,[]));
% 单层 layers 写法和只有初始站也要有定义良好的空界面输出。
mat.layers={a}; mat.contactResistance=[];
r=calcWallConduction1D(400,1200,[0 .1],mat,300,[]);
assert(isempty(r.interfaceFinal.energyResidual));
r=calcWallConduction1D(400,1200,0,mat,300,[]);
assert(~r.interfaceFinal.available && r.surfaceTemperature==300);
% 三层混合界面与不同间距，验证不把多于两层的质量/通量接错。
c=rmfield(a,'propertyTable'); c.thickness=.009; c.thermalNodes=4;
mat.layers={rmfield(a,'propertyTable'),rmfield(b,'propertyTable'),c};
mat.contactResistance=[0 .001];
r=calcWallConduction1D(300,1000,0:.1:1,mat,400,[]);
assert(r.thermalNodes==9 && nnz(diff(r.depth)==0)==1);
assert(max(abs(r.interfaceFinal.energyResidual),[],'all')<1e-5);
mass1=1200*800*.002*[.5 1 .5]; mass2=2700*900*.002*[.5 1 .5];
mass3=1200*800*.003*[.5 1 1 .5];
expectedMass=[mass1(1:2),mass1(3)+mass2(1),mass2(2:3),mass3];
assert(abs(sum(expectedMass.*(r.finalDepthTemperature-400))- ...
    sum(300*(1000-r.surfaceTemperature(2:end))*.1))<1e-4);
fprintf('Layered conduction passed: homogeneous limit, matrix exponential, nonlinear ODE, interface/storage balances and input domains.\n');

    function d=derivative(~,Y)
        cap=[1200*.002*[.5;1;.5].*(500+.5*Y(1:3)); ...
            2700*.002*[.5;1;.5].*(700+.2*Y(4:6))];
        flux=[diff(Y(1:3)).*(1+.0005*(Y(1:2)+Y(2:3)))/.002; ...
            (Y(4)-Y(3))/.002; ...
            diff(Y(4:6)).*(4+.001*(Y(4:5)+Y(5:6)))/.002];
        d=[400*(1200-Y(1))+flux(1);diff(flux);-flux(end)]./cap;
    end
end

function [mass,g]=independentGrid(a,b,rc)
% 仅为这个两层三节点算例独立写出热容与导热连接。
ma=a.density*a.specificHeat*.002; mb=b.density*b.specificHeat*.002;
if rc==0
    mass=[ma/2;ma;(ma+mb)/2;mb;mb/2];
    g=[a.thermalConductivity a.thermalConductivity b.thermalConductivity b.thermalConductivity]/.002;
else
    mass=[ma/2;ma;ma/2;mb/2;mb;mb/2];
    g=[a.thermalConductivity/.002 a.thermalConductivity/.002 1/rc b.thermalConductivity/.002 b.thermalConductivity/.002];
end
end

function mustFail(action)
failed=false; try, action(); catch, failed=true; end
assert(failed,'非法输入没有被拒绝。');
end
