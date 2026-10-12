function testHeadDomainBlending()
%TESTHEADDOMAINBLENDING 验证 LT 适用边界附近不再硬切换整条下游状态。
% 模块作用：让层流交接 Re_theta 穿过下边界 300 与上边界 1e6，
% 检查 Cf、theta、H 及下游状态连续；衔接必须有独立来源标记。
addpath(fileparts(fileparts(mfilename('fullpath'))));
nu=1.5e-5; n=3;
m.F=repmat([1 2 3],n,1);
m.boundaryLayerEdge=struct('method','受控外缘','speed',100*ones(n,1));
m.surfaceStreamlines=struct('velocitySource','受控外缘', ...
    'distance',[0.1;0.11;0.12],'pathCost',[0;0.01;0.02], ...
    'predecessorFace',[0;1;2],'stepLength',[0;0.01;0.01], ...
    'connectedFace',true(n,1));
t=struct('laminarFace',[true;false;false],'turbulentFace',[false;true;true]);
for edge=[300 1e6]
    values=cell(1,2);
    for k=1:2
        re=edge*(1+(2*k-3)*1e-7);
        m.boundaryLayer.momentumThickness=[re*nu/100;NaN;NaN];
        values{k}=calcTurbulentBoundaryLayerState(m,nu*ones(n,1),t);
    end
    a=values{1}; b=values{2};
    assert(all(a.validFace(2:3)) && all(b.validFace(2:3)));
    assert(max(abs(a.skinFrictionCoefficient(2:3)-b.skinFrictionCoefficient(2:3)))<1e-7);
    assert(max(abs(a.momentumThickness(2:3)-b.momentumThickness(2:3)))<1e-7);
    assert(max(abs(a.shapeFactor(2:3)-b.shapeFactor(2:3)))<1e-5);
end
m.boundaryLayer.momentumThickness=[400*nu/100;NaN;NaN];
blend=calcTurbulentBoundaryLayerState(m,nu*ones(n,1),t);
assert(blend.sourceID(2)==7 && blend.integralWeight(2)>0 && blend.integralWeight(2)<1);
assert(blend.summary.blendedFaces>0);
fprintf('Head domain blending passed: both correlation limits and downstream continuity.\n');
end
