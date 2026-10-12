function testEulerEdgeConvergence()
%TESTEULEREDGECONVERGENCE 以独立 Taylor-Maccoll 解检验近壁状态及梯度。
% 直圆锥外缘状态沿母线恒定，因此真值 dU/ds=0；不由 Euler 拟合参考值。
% 此测试只证明该锥流算例的收敛趋势，不认证实际箭体或钝鼻驻点。
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);
reference=solveConicalShock(fs,10);assert(reference.attached);
counts=[32 64 128];errors=zeros(3,4);gradientError=zeros(3,1);solutions=cell(3,1);
for k=1:3
    g=buildAxisymmetricEulerGrid([-1 tand(10);0 0], ...
        struct('axialCells',counts(k),'radialCells',round(.75*counts(k)),'radialStretch',1));
    sol=solveAxisymmetricEuler(g,fs);assert(sol.available);solutions{k}=sol;
    edge=extractAxisymmetricEulerEdge(sol);selected=edge.x<-.35 & edge.x>-.8;
    names={'pressure','temperature','speed','totalPressure'};
    for j=1:numel(names)
        errors(k,j)=abs(mean(edge.(names{j})(selected))/reference.(names{j})-1);
    end
    gradientError(k)=max(abs(edge.velocityGradient(selected)))*g.length/fs.speed;
    fprintf('Cone edge %d: p %.6g, T %.6g, U %.6g, p0 %.6g, gradient %.6g\n', ...
        counts(k),errors(k,:),gradientError(k));
    assert(~edge.available && ~edge.gridAccuracyVerified);
end
% 独立精度要求：最细网格压力/总压 5%，温度/速度 1%，
% 速度梯度绝对误差不超过来流速度/模型长度的 2%；各项应随加密下降。
assert(all(all(diff(errors,1,1)<0)) && all(diff(gradientError)<0));
assert(all(errors(end,:)<[.05 .01 .01 .05]) && gradientError(end)<.02);
comparison=compareAxisymmetricEulerEdges(solutions{2},solutions{3});
assert(comparison.available && ~comparison.edgeStateAvailable);
assert(comparison.axialNodesChanged);
fprintf('Euler edge cone convergence passed; blunt-nose stagnation velocity remains unverified.\n');
end
