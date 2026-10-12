function testAerothermalRegion()
% 区域分类回归测试：验证方向符号、阈值边界、数据保持及实际模型接口。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg = createCaseConfig();'); cfg.output.verbose = false;

%% 1. 构造封闭方柱，其六个表面法向沿坐标轴，便于解析验证
m.V = [0 -1 -1;0 1 -1;0 1 1;0 -1 1; ...
       1 -1 -1;1 1 -1;1 1 1;1 -1 1];
m.F = [1 2 6;1 6 5;2 3 7;2 7 6;3 4 8;3 8 7; ...
       4 1 5;4 5 8;1 4 3;1 3 2;5 6 7;5 7 8];
evalc('m = calcFaceGeometry(m); m = fixNormalDirection(m,cfg);');
evalc('m = buildSurfaceCoordinates(m,cfg);');
m = analyzeLocalGeometry(m,cfg); before = m;
flight = cfg.flight; flight.alpha = 0; flight.beta = 0;
fs = calcFreestream(flight);
m = classifyAerothermalRegion(m,fs,cfg); r = m.aerothermalRegion;
assert(isequaln(rmfield(m,'aerothermalRegion'),before));
assert(nnz(r.isWindward) == 2 && nnz(r.isLeeward) == 2 && nnz(r.isGrazing) == 8);
assert(all(r.isWindward(r.isNoseCap)) && all(r.isLeeward(r.isAftCap)));
assert(all(r.incidenceAngleDeg(r.isNoseCap) == 90));
assert(all(r.incidenceAngleDeg(r.isAftCap) == -90));
assert(abs(r.summary.windwardArea+r.summary.leewardArea+ ...
    r.summary.grazingArea-r.summary.totalArea) < 1e-12);

%% 2. 正负攻角和侧滑角：使用正式来流模块验证符号而非手填方向
N = m.localGeometry.normal;
for angle = [-10 10]
    flight.alpha = angle; flight.beta = 0;
    tilted = classifyAerothermalRegion(m,calcFreestream(flight),cfg);
    rt = tilted.aerothermalRegion;
    assert(all(rt.isWindward(N(:,3) == -sign(angle))));
    assert(all(rt.isLeeward(N(:,3) == sign(angle))));
    flight.alpha = 0; flight.beta = angle;
    tilted = classifyAerothermalRegion(m,calcFreestream(flight),cfg);
    rt = tilted.aerothermalRegion;
    assert(all(rt.isWindward(N(:,2) == sign(angle))));
    assert(all(rt.isLeeward(N(:,2) == -sign(angle))));
end

%% 3. 阈值边界、反向流、向量归一化与非法输入
flight.alpha = 0.5; flight.beta = 0;
edge = classifyAerothermalRegion(m,calcFreestream(flight),cfg);
assert(all(edge.aerothermalRegion.isGrazing(abs(N(:,3)) == 1)));
flight.alpha = 0.501;
edge = classifyAerothermalRegion(m,calcFreestream(flight),cfg);
assert(all(edge.aerothermalRegion.isWindward(N(:,3) == -1)));
reverse = fs; reverse.unitFlowDirection = -fs.unitFlowDirection;
rev = classifyAerothermalRegion(m,reverse,cfg);
assert(isequal(rev.aerothermalRegion.isWindward,r.isLeeward));
scaled = fs; scaled.unitFlowDirection = 5*fs.unitFlowDirection;
scaledMesh = classifyAerothermalRegion(m,scaled,cfg);
assert(isequal(scaledMesh.aerothermalRegion.mu,r.mu));
bad = fs; bad.unitFlowDirection = [0 0 0];
mustFail(@()classifyAerothermalRegion(m,bad,cfg),'RocketAeroHeat:Region:FlowDirection');
badMesh = m; badMesh.F(1,:) = badMesh.F(1,[1 3 2]);
mustFail(@()classifyAerothermalRegion(badMesh,fs,cfg),'RocketAeroHeat:Region:StaleGeometry');
badMesh = rmfield(m,'localGeometry');
mustFail(@()classifyAerothermalRegion(badMesh,fs,cfg),'RocketAeroHeat:Region:Geometry');

%% 4. 默认 STL 全链路：检查分类完备性和重复调用不修改几何
evalc('raw = loadSTLMesh(cfg); raw = calcFaceGeometry(raw);');
[raw,~] = alignGeometryAxis(raw,'pca');
evalc('m = autoRemeshSurface(raw,cfg); m = fixNormalDirection(m,cfg);');
evalc('m = buildSurfaceCoordinates(m,cfg);');
m = analyzeLocalGeometry(m,cfg); before = m;
tic; m = classifyAerothermalRegion(m,calcFreestream(cfg.flight),cfg); elapsed = toc;
r = m.aerothermalRegion;
assert(all(double(r.isWindward)+double(r.isLeeward)+double(r.isGrazing) == 1));
assert(all(isfinite(r.mu)) && all(abs(r.mu) <= 1));
assert(isequaln(rmfield(m,'aerothermalRegion'),before));
fprintf('区域分类测试通过：%d 面元；迎风 %d，背风 %d，近切向 %d；耗时 %.3f 秒。\n', ...
    size(m.F,1),nnz(r.isWindward),nnz(r.isLeeward),nnz(r.isGrazing),elapsed);
end

function mustFail(callback,identifier)
% 用途：确认非法输入被明确拒绝，而不是返回看似有效的分类。
try
    callback();
catch ME
    assert(strcmp(ME.identifier,identifier),'错误类型不符：%s',ME.message); return;
end
error('预期错误未发生：%s',identifier);
end
