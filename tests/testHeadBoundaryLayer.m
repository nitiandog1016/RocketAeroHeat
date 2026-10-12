function testHeadBoundaryLayer()
%TESTHEADBOUNDARYLAYER 验证 Head 守恒、梯度趋势、收敛与分离传递。
% 模块一：独立积算动量和卷吸通量，检查程序满足两个积分方程。
% 模块二：同一线性外缘多次细分，检查可变 H 与分离位置收敛。
% 模块三：分离后含热流输出空白和重新加速，不能本地重启附着解。
addpath(fileparts(fileparts(mfilename('fullpath'))));
nu=1.5e-5; th0=5e-4;
flat=marchHeadBoundaryLayerSegment(th0,1.4,100,100,nu,nu,0.3);
fast=marchHeadBoundaryLayerSegment(th0,1.4,100,110,nu,nu,0.3);
slow=marchHeadBoundaryLayerSegment(th0,1.4,100,70,nu,nu,0.3);
assert(all([flat.status fast.status slow.status]==0));
assert(fast.theta<flat.theta && slow.theta>flat.theta);
assert(fast.shapeFactor<flat.shapeFactor && slow.shapeFactor>flat.shapeFactor);
assert(slow.shapeFactor>1.6 && slow.shapeFactor<2.4);
for terminal=[100 70]
    [theta,H,~,~]=profile(400,terminal);
    x=linspace(0,0.3,401)'; U=linspace(100,terminal,401)';
    cf=0.246*10.^(-0.678*H).*(U.*theta/nu).^(-0.268);
    H1=3.3+0.8234./(H-1.1).^1.287;
    coefficient=0.8234/0.5^1.287*(1.6-0.6778)^3.064;
    ids=H>1.6; H1(ids)=3.3+coefficient./(H(ids)-0.6778).^3.064;
    E=0.0306./(H1-3).^0.6169;
    momentum=trapz(x,cf/2-(H+2).*theta./U*(terminal-100)/0.3);
    entrainment=trapz(x,U.*E);
    assert(abs(theta(end)-theta(1)-momentum)<2e-7);
    % 跨 H=1.6 的两段闭合必须保持卷吸通量连续。
    flux=U.*theta.*H1;
    assert(abs(flux(end)-flux(1)-entrainment)<5e-4*abs(entrainment));
    coarse=marchHeadBoundaryLayerSegment(th0,1.4,100,terminal,nu,nu,0.3);
    assert(abs(coarse.theta/theta(end)-1)<1e-4);
    assert(abs(coarse.shapeFactor-H(end))<1e-4);
end
[~,~,s1,x1]=profile(10,40); [~,~,s2,x2]=profile(100,40);
assert(s1==1 && s2==1 && abs(x1-x2)<2e-5);
domain=marchHeadBoundaryLayerSegment(1e-6,1.4,100,100,nu,nu,0.1);
assert(domain.status==2 && isnan(domain.Cf));
bad=marchHeadBoundaryLayerSegment(th0,1.05,100,100,nu,nu,0.1);
assert(bad.status==3 && isnan(bad.Cf));

%% 构造有分支、空白和强减速的准流线，核对无自动再附与分支隔离
n=103;
m.F=repmat([1 2 3],n,1);
m.boundaryLayerEdge.method='解析线性外缘';
m.boundaryLayerEdge.speed=[100;linspace(100,30,60)';100*ones(42,1)];
p=struct('velocitySource',m.boundaryLayerEdge.method, ...
    'distance',(1:n)'*0.005,'pathCost',(1:n)'*0.005, ...
    'predecessorFace',[0;(1:n-1)'],'stepLength',0.005*ones(n,1), ...
    'connectedFace',true(n,1));
% 最后一个面元是第一面元的独立分支，不应受到另一分支分离影响。
p.predecessorFace(end)=1;
m.surfaceStreamlines=p;
m.boundaryLayer.momentumThickness=nan(n,1);
m.boundaryLayer.momentumThickness(1)=th0;
t=struct('laminarFace',[true;false(n-1,1)], ...
    'turbulentFace',[false;true(n-1,1)]);
t.turbulentFace(70)=false;
state=calcTurbulentBoundaryLayerState(m,nu*ones(n,1),t);
first=find(state.possibleSeparationFace,1);
assert(~isempty(first) && first<60);
assert(all(state.separationUnresolvedFace(first:n-1)));
assert(~state.separationUnresolvedFace(n) && state.validFace(n));
assert(all(isnan(state.skinFrictionCoefficient(first:n-1))));
assert(~any(state.numericalUnresolvedFace));
assert(state.sourceID(70)==5 && state.sourceID(71)==5);
fprintf('Head BL passed: balances, refinement, variable H, separation, gaps and branches.\n');

    function [theta,H,status,position]=profile(count,lastU)
        theta=nan(count+1,1); H=theta; theta(1)=th0; H(1)=1.4;
        velocities=linspace(100,lastU,count+1); dx=0.3/count;
        status=0; position=0.3;
        for j=1:count
            r=marchHeadBoundaryLayerSegment(theta(j),H(j), ...
                velocities(j),velocities(j+1),nu,nu,dx);
            if r.status~=0
                status=r.status; position=(j-1+r.fraction)*dx; return;
            end
            theta(j+1)=r.theta; H(j+1)=r.shapeFactor;
        end
    end
end
