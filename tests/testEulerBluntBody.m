function testEulerBluntBody(refine)
%TESTEULERBLUNTBODY 半球头圆柱的驻点压力及总焓独立检查。
% 默认验证较粗网格；传 true 额外执行细网格，可能需要数分钟。
% 驻点参考值由法向激波跃迁和激波后等熵减速获得，非 Euler 结果拟合。
if nargin<1,refine=false;end
evalc('cfg=createCaseConfig();');fs=calcFreestream(cfg.flight);
theta=linspace(pi/2,0,65)';contour=[-1,.1;.1*(cos(theta)-1),.1*sin(theta)];
reference=solveNormalShock(fs);
if refine,counts=[48 96];else,counts=48;end
errors=zeros(numel(counts),2);
for k=1:numel(counts)
    g=buildAxisymmetricEulerGrid(contour,struct('axialCells',120,'radialCells',counts(k)));
    result=solveAxisymmetricEuler(g,fs,struct('maxIterations',150000));
    tip=find(g.bodyWall,1,'last');
    errors(k,:)=[result.wall.pressure(tip),max(result.wall.pressure)]/reference.stagnationPressure-1;
    assert(result.available && result.conservation.identityError<1e-10);
    assert(result.diagnostics.maximumTotalEnthalpyDeviation<1e-6);
    assert(abs(errors(k,1))<.02);
    assert(result.conservation.wallMassLeak<1e-12 && result.conservation.wallEnergyLeak<1e-12);
    fprintf('Sphere %d radial cells: stagnation error=%.6g, peak error=%.6g, H error=%.6g\n', ...
        counts(k),errors(k,1),errors(k,2),result.diagnostics.maximumTotalEnthalpyDeviation);
end
if refine
    assert(all(abs(errors(2,:))<abs(errors(1,:))));
    assert(abs(errors(2,2))<.03);
end
fprintf('Euler blunt-body benchmark passed; field-wide grid independence remains unverified.\n');
end
