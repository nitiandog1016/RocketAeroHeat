function testLocalGeometry()
% Analytic cylinder/cone checks plus the actual project STL pipeline.
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg = createCaseConfig();'); cfg.output.verbose = false;
for cone = [false true]
    m = body(32,cone);
    evalc('m = calcFaceGeometry(m); m = fixNormalDirection(m,cfg);');
    evalc('m = buildSurfaceCoordinates(m,cfg);');
    before = m;
    m = analyzeLocalGeometry(m,cfg); g = m.localGeometry;
    assert(isequaln(rmfield(m,'localGeometry'),before));
    valid = g.frameValid;
    assert(max(abs(vecnorm(g.tangentFromNose(valid,:),2,2)-1)) < 1e-12);
    assert(max(abs(sum(g.tangentFromNose(valid,:).*g.normal(valid,:),2))) < 1e-12);
    assert(max(vecnorm(cross(g.tangentFromNose(valid,:), ...
        g.tangentCircumferential(valid,:),2)-g.normal(valid,:),2,2)) < 1e-12);
    assert(all(g.tangentFromNose(valid,1) < 0));
    assert(all(isnan(g.tangentFromNose(~valid,:)),'all'));
    assert(nnz(g.isAftCap) == 32);
    assert(g.summary.boundaryEdges == 0 && g.summary.nonManifoldEdges == 0);
    if cone
        assert(nnz(g.isNoseTipAdjacent) == 32);
        expected = atand(cos(pi/32));
        assert(max(abs(g.inclinationDeg(g.isSidewall)-expected)) < 1e-10);
    else
        assert(nnz(g.isNoseCap) == 32);
        assert(max(abs(g.inclinationDeg(g.isSidewall))) < 1e-12);
        assert(all(g.inclinationDeg(g.isNoseCap) == 90));
        assert(all(g.inclinationDeg(g.isAftCap) == -90));
    end
end
bad = before; bad.F(1,2) = bad.F(1,1);
mustFail(@()analyzeLocalGeometry(bad,cfg),'RocketAeroHeat:LocalGeometry:DegenerateFace');
bad = rmfield(before,'surfaceCoord');
mustFail(@()analyzeLocalGeometry(bad,cfg),'RocketAeroHeat:LocalGeometry:SurfaceCoordinates');
evalc('raw = loadSTLMesh(cfg); raw = calcFaceGeometry(raw);');
[raw,~] = alignGeometryAxis(raw,'pca');
evalc('m = autoRemeshSurface(raw,cfg); m = fixNormalDirection(m,cfg);');
evalc('m = buildSurfaceCoordinates(m,cfg);');
old = m; tic; m = analyzeLocalGeometry(m,cfg); elapsed = toc;
assert(isequaln(rmfield(m,'localGeometry'),old));
g = m.localGeometry;
assert(all(isfinite(g.inclinationDeg)));
assert(all(isfinite(g.minimumAngleDeg)));
assert(g.summary.boundaryEdges == 0 && g.summary.nonManifoldEdges == 0);
evalc('q = checkMeshQuality(m,cfg);');
fprintf('Local geometry tests passed: %d faces, %.3f s; minimum angle %.4f deg, face %d.\n', ...
    size(m.F,1),elapsed,g.summary.minimumAngleDeg,g.summary.worstFaceID);
fprintf('Caps: aft %d, nose %d; nose-tip adjacent %d; small angles %d.\n', ...
    nnz(g.isAftCap),nnz(g.isNoseCap),nnz(g.isNoseTipAdjacent),nnz(g.poorAngle));
end

function m = body(n,cone)
theta = (0:n-1)'*2*pi/n; ring = [cos(theta) sin(theta)];
if cone
    m.V = [zeros(n,1) ring;0.5*ones(n,1) 0.5*ring;1 0 0;0 0 0];
    m.F = zeros(4*n,3);
    for i=1:n
        j=mod(i,n)+1;
        m.F(2*i-1,:)=[i j n+i]; m.F(2*i,:)=[j n+j n+i];
        m.F(2*n+i,:)=[n+i n+j 2*n+1]; m.F(3*n+i,:)=[i 2*n+2 j];
    end
else
    m.V = [zeros(n,1) ring;ones(n,1) ring;0 0 0;1 0 0];
    m.F = zeros(4*n,3);
    for i=1:n
        j=mod(i,n)+1;
        m.F(2*i-1,:)=[i j n+i]; m.F(2*i,:)=[j n+j n+i];
        m.F(2*n+i,:)=[i 2*n+1 j]; m.F(3*n+i,:)=[n+i n+j 2*n+2];
    end
end
end

function mustFail(f,id)
try
    f();
catch ME
    assert(strcmp(ME.identifier,id),'Unexpected error: %s',ME.message); return;
end
error('Expected error %s.',id);
end
