function testNumericalRegression()
% Deterministic regressions for the second complete code review.
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg = createCaseConfig();');

% Cached atmospheric data must not survive an altitude change.
flight = cfg.flight;
flight.atmosphere = standardAtmosphere(0);
fs = calcFreestream(flight);
expected = standardAtmosphere(flight.altitude);
assert(abs(fs.density-expected.density) < 1e-12);
flight.alpha = [0 1];
mustFail(@()calcFreestream(flight));
flight = cfg.flight;
flight.atmosphere = expected;
flight.atmosphere.gamma = 1;
mustFail(@()calcFreestream(flight));

% Correctly retain a sharp fold, even if one input face has reversed winding.
fold.V = [0 0 0;1 0 0;0 1 0;0 cosd(30) sind(30)];
fold.F = [1 2 3;2 1 4];
evalc('fold = calcFaceGeometry(fold); c = analyzeSurfaceCurvature(fold,cfg);');
assert(abs(c.edgeInfo.angleDeg-150) < 1e-10);
fold.F(2,:) = fold.F(2,[1 3 2]);
evalc('fold = calcFaceGeometry(fold); c = analyzeSurfaceCurvature(fold,cfg);');
assert(abs(c.edgeInfo.angleDeg-150) < 1e-10);
mustFail(@()alignGeometryAxis(fold,'matrix',nan(3)));
mustFail(@()alignGeometryAxis(fold,'matrix',eye(3),[0 Inf 0]));

% Translation and initially inward winding must not affect signed volume.
tetra.V = [0 0 0;1 0 0;0 1 0;0 0 1] + 1e8;
tetra.F = [1 3 2;1 2 4;2 3 4;3 1 4];
for reverse = [false true]
    t = tetra;
    if reverse, t.F = t.F(:,[1 3 2]); end
    evalc('t = fixNormalDirection(t,cfg);');
    assert(abs(t.normalInfo.componentSignedVolume-1/6) < 1e-12);
    assert(isequal(t.F,tetra.F));
end

% Synthetic circular body: curvature sensitivity should refine, not coarsen.
body = cylinderMesh(32,0.1);
small = cfg;
small.remesh.hMin = 0.001;
small.remesh.hMax = 0.1;
small.remesh.maxStations = 3;
small.remesh.maxCircumPoints = 256;
evalc('coarse = autoRemeshSurface(body,small);');
fineCfg = small;
fineCfg.remesh.curvatureSensitivity = 2;
evalc('fine = autoRemeshSurface(body,fineCfg);');
assert(all(fine.param.sectionPointCount > coarse.param.sectionPointCount));
assert(max(abs(fine.param.sectionTargetSpacing*2-coarse.param.sectionTargetSpacing)) < 1e-12);
evalc('q = checkMeshQuality(fine,small);');
assert(abs(q.meanAspectRatio-fine.remeshInfo.meanAspectRatio) < 1e-12);
assert(abs(q.maxAspectRatio-fine.remeshInfo.maxAspectRatio) < 1e-12);

limited = small;
limited.remesh.maxFaces = 10;
mustFail(@()autoRemeshSurface(body,limited),'RocketAeroHeat:Remesh:FaceLimit');
limited = small;
limited.remesh.curvatureSensitivity = NaN;
mustFail(@()autoRemeshSurface(body,limited));

% Never drop a second body or fabricate missing cross-section edges.
multiple.V = [body.V;body.V+[0 1 0]];
multiple.F = [body.F;body.F+size(body.V,1)];
mustFail(@()autoRemeshSurface(multiple,small), ...
    'RocketAeroHeat:Remesh:MultipleSectionLoops');
open = body;
open.F(1:2,:) = [];
mustFail(@()autoRemeshSurface(open,small),'RocketAeroHeat:Remesh:OpenSection');
gap.V = [body.V;body.V+[2 0 0]];
gap.F = multiple.F;
mustFail(@()autoRemeshSurface(gap,small),'RocketAeroHeat:Remesh:InvalidSection');

% Fallback inference must separate flat-cap centers from perimeter nodes.
square = cylinderMesh(4,1);
evalc('inferred = buildSurfaceCoordinates(square,cfg);');
assert(all(inferred.surfaceCoord.sectionPointCount == 4));
assert(inferred.surfaceCoord.nodeIsCap(9) && inferred.surfaceCoord.nodeIsCap(10));
assert(inferred.surfaceCoord.nodeSFromAft(9) == 0);
assert(inferred.surfaceCoord.nodeSFromNose(10) == 0);

% A single nose cap must not also be assigned as an aft cap.
oneCap.V = [square.V(1:8,:);2 0 0];
oneCap.F = [square.F(1:8,:);5 6 9;6 7 9;7 8 9;8 5 9];
evalc('oneCap = buildSurfaceCoordinates(oneCap,cfg);');
assert(max(abs(oneCap.surfaceCoord.nodeSFromAft(1:4))) < 1e-12);
assert(oneCap.surfaceCoord.nodeSFromNose(9) == 0);

% Periodic interpolation also works when the first ring fraction is not zero.
square.param.sectionX = [0;1];
square.param.sectionNodeIDs = {(1:4)';(5:8)'};
square.param.sectionCenterYZ = zeros(2,2);
square.param.nodeSectionID = [ones(4,1);2*ones(4,1);0;0];
square.param.nodeRingFraction = [repmat((0.1:0.25:0.85)',2,1);NaN;NaN];
square.param.aftCapNodeID = 9;
square.param.noseCapNodeID = 10;
evalc('square = buildSurfaceCoordinates(square,cfg);');
seamPoint = squeeze(square.surfaceCoord.sectionPointGrid(1,8,:));
assert(norm(seamPoint-[0;0.1;-0.9]) < 1e-12);

% Full current project pipeline, including surface coordinates and flow.
evalc('raw = loadSTLMesh(cfg); raw = calcFaceGeometry(raw);');
[raw,~] = alignGeometryAxis(raw,'pca');
evalc('m = autoRemeshSurface(raw,cfg); m = fixNormalDirection(m,cfg);');
evalc('m = buildSurfaceCoordinates(m,cfg); q = checkMeshQuality(m,cfg);');
assert(m.normalInfo.isClosedSurface && m.normalInfo.statusOK);
assert(all(isfinite(m.surfaceCoord.nodeS)));
yc = interp1(m.param.sectionX,m.param.sectionCenterYZ(:,1),m.center(:,1),'linear','extrap');
zc = interp1(m.param.sectionX,m.param.sectionCenterYZ(:,2),m.center(:,1),'linear','extrap');
expectedTheta = atan2(m.center(:,3)-zc,m.center(:,2)-yc);
angleError = angle(exp(1i*(expectedTheta-m.surfaceCoord.faceTheta)));
assert(max(abs(angleError)) < 1e-12);
assert(abs(q.meanAspectRatio-m.remeshInfo.meanAspectRatio) < 1e-12);
fprintf('Numerical regressions passed. Default faces=%d, poor faces=%.3f%%.\n', ...
    size(m.F,1),q.badFacePercent);
end

function mesh = cylinderMesh(n,radius)
theta = (0:n-1)'*2*pi/n;
ring = radius*[cos(theta) sin(theta)];
mesh.V = [zeros(n,1) ring;ones(n,1) ring;0 0 0;1 0 0];
mesh.F = zeros(4*n,3);
for i = 1:n
    j = mod(i,n)+1;
    mesh.F(2*i-1,:) = [i j n+i];
    mesh.F(2*i,:) = [j n+j n+i];
    mesh.F(2*n+i,:) = [i 2*n+1 j];
    mesh.F(3*n+i,:) = [n+i n+j 2*n+2];
end
end

function mustFail(callback,identifier)
caught = false;
try
    evalc('callback();');
catch ME
    caught = true;
    if nargin > 1
        assert(strcmp(ME.identifier,identifier), ...
            'Expected %s, received %s: %s',identifier,ME.identifier,ME.message);
    end
end
assert(caught,'Expected this invalid input to be rejected.');
end
