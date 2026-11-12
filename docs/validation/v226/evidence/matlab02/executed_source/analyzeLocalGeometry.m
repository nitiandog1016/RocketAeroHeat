function mesh = analyzeLocalGeometry(mesh,cfg)
%ANALYZELOCALGEOMETRY Final-mesh geometric diagnostics (coordinates in m).
% Call after fixNormalDirection and buildSurfaceCoordinates:
%   mesh = analyzeLocalGeometry(mesh,cfg);
% Does not change V, F, surfaceCoord or remeshing parameters.
%
% mesh.localGeometry fields (one row per face unless stated otherwise):
%   normal                  unit normal from current face winding
%   tangentFromNose          projection of -X onto the face
%   tangentToNose            opposite tangent
%   tangentCircumferential   cross(normal,tangentFromNose)
% The frame [tangentFromNose,tangentCircumferential,normal] is right handed.
% For an outward axisymmetric sidewall, circumferential points toward
% increasing theta = atan2(z-zc,y-yc). These are geometric directions,
% NOT streamlines or principal curvature directions on arbitrary bodies.
%   frameValid              false where X is normal to the face (caps)
%   inclinationDeg          asin(normal_x), signed, [-90,90] degrees;
%                           positive for a nose-facing surface
%   isAftCap/isNoseCap       triangles on an axial extreme plane
%   isAftTipAdjacent/isNoseTipAdjacent : incident on an isolated end vertex
%   region                  sidewall/aftCap/noseCap/aftTip/noseTip
%                           (topological end labels, not flow regions)
%   minimumAngleDeg/poorAngle : diagnostic threshold 20 degrees by default
%   discreteCurvatureMean/Max, discreteRadiusProxy : reused edge metric
% IMPORTANT: discreteRadiusProxy is NOT a principal, meridional or physical
% nose radius. It must not be used directly in stagnation-heating formulas.
% NaN denotes an unavailable frame; Inf radius proxy denotes zero metric.
% cfg.localGeometry.minimumAngleDeg optionally changes the diagnostic flag.
% Outward winding is a prerequisite; this routine never silently flips it.

if nargin < 2, cfg = struct; end
if ~isstruct(mesh) || ~isfield(mesh,'V') || ~isfield(mesh,'F')
    error('RocketAeroHeat:LocalGeometry:Input','mesh.V and mesh.F are required.');
end
V = mesh.V; F = mesh.F;
validateattributes(V,{'numeric'},{'real','finite','2d','ncols',3,'nonempty'});
validateattributes(F,{'numeric'},{'real','finite','integer','positive','2d','ncols',3,'nonempty','<=',size(V,1)});
V = double(V); F = double(F); nf = size(F,1);
if ~isfield(mesh,'surfaceCoord')
    error('RocketAeroHeat:LocalGeometry:SurfaceCoordinates', ...
        'Run buildSurfaceCoordinates on the final mesh first.');
end
sc = mesh.surfaceCoord;
for name = {'faceX','faceR','faceTheta','faceS','faceSFromNose','faceSFromAft'}
    key = name{1};
    if ~isfield(sc,key) || ~isequal(size(sc.(key)),[nf 1]) || any(~isfinite(sc.(key)))
        error('RocketAeroHeat:LocalGeometry:SurfaceCoordinates', ...
            'Missing or invalid surfaceCoord.%s; rebuild surface coordinates.',key);
    end
end
threshold = 20;
if isfield(cfg,'localGeometry') && isfield(cfg.localGeometry,'minimumAngleDeg')
    threshold = cfg.localGeometry.minimumAngleDeg;
end
validateattributes(threshold,{'numeric'},{'real','finite','scalar','>=',0,'<=',60});

% Recompute from current connectivity; do not trust stale face normals.
p1 = V(F(:,1),:); p2 = V(F(:,2),:); p3 = V(F(:,3),:);
e12 = p2-p1; e13 = p3-p1; e23 = p3-p2;
raw = cross(e12,e13,2); twiceArea = vecnorm(raw,2,2);
edgeLengths = [vecnorm(e12,2,2),vecnorm(e23,2,2),vecnorm(e13,2,2)];
if any(twiceArea <= 64*eps*max(edgeLengths,[],2).^2)
    error('RocketAeroHeat:LocalGeometry:DegenerateFace', ...
        'Degenerate triangles must be repaired before local analysis.');
end
N = raw./twiceArea;
C = p1+(e12+e13)/3;
tol = max(1e-10*max(max(V)-min(V)),64*eps(max(abs(V(:)))));
if any(abs(sc.faceX-C(:,1)) > tol)
    error('RocketAeroHeat:LocalGeometry:SurfaceCoordinates', ...
        'Surface coordinates do not match current mesh; rebuild them.');
end
t = repmat([-1 0 0],nf,1)+N(:,1).*N;
tnorm = vecnorm(t,2,2); frameValid = tnorm > 1e-10;
t(frameValid,:) = t(frameValid,:)./tnorm(frameValid);
t(~frameValid,:) = NaN;
c = cross(N,t,2);

% End-plane detection includes caps triangulated without a center vertex.
x = V(:,1); xmin = min(x); xmax = max(x);
if xmax-xmin <= tol
    error('RocketAeroHeat:LocalGeometry:AxialExtent','Body must have nonzero X extent.');
end
aftNodes = abs(x-xmin) <= tol; noseNodes = abs(x-xmax) <= tol;
aftCap = all(aftNodes(F),2); noseCap = all(noseNodes(F),2);
aftTip = false(nf,1); noseTip = false(nf,1);
if nnz(aftNodes) == 1, aftTip = any(aftNodes(F),2) & ~aftCap; end
if nnz(noseNodes) == 1, noseTip = any(noseNodes(F),2) & ~noseCap; end
region = repmat({'sidewall'},nf,1);
region(aftCap) = {'aftCap'}; region(noseCap) = {'noseCap'};
region(aftTip) = {'aftTip'}; region(noseTip) = {'noseTip'};

% atan2 avoids loss of precision for small angles.
angles = [atan2d(twiceArea,sum(e12.*e13,2)), ...
          atan2d(twiceArea,sum(-e12.*e23,2)), ...
          atan2d(twiceArea,sum(e13.*e23,2))];
minimumAngle = min(angles,[],2);

% Reuse the existing curvature algorithm on the FINAL mesh. Keep its
% target-size suggestions local so no remesh settings/results are changed.
curvatureCfg = cfg;
defaults = struct('featureAngle',30,'curvatureAngleTarget',10, ...
    'curvatureSensitivity',1,'hMin',0.001,'hMax',0.1,'featureSizeFactor',1);
if ~isfield(curvatureCfg,'remesh'), curvatureCfg.remesh = struct; end
for name = fieldnames(defaults)'
    key = name{1};
    if ~isfield(curvatureCfg.remesh,key)
        curvatureCfg.remesh.(key) = defaults.(key);
    end
end
tmp = struct('V',V,'F',F,'normal',N);
evalc('tmp = analyzeSurfaceCurvature(tmp,curvatureCfg);');
radius = inf(nf,1); curved = tmp.curvatureMax > 0;
radius(curved) = 1./tmp.curvatureMax(curved);
edges = sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2);
[~,~,groups] = unique(edges,'rows'); counts = accumarray(groups,1);
badEdge = counts(groups) ~= 2;
topologyIssue = any(reshape(badEdge,nf,3),2);

g = struct;
g.version = '1.0';
g.normal = N; g.center = C; g.area = twiceArea/2;
g.tangentFromNose = t; g.tangentToNose = -t;
g.tangentCircumferential = c; g.frameValid = frameValid;
g.inclinationDeg = asind(max(-1,min(1,N(:,1))));
g.isAftCap = aftCap; g.isNoseCap = noseCap;
g.isAftTipAdjacent = aftTip; g.isNoseTipAdjacent = noseTip;
g.isSidewall = ~(aftCap|noseCap); g.region = region;
g.x = sc.faceX; g.r = sc.faceR; g.theta = sc.faceTheta;
g.sFromNose = sc.faceSFromNose; g.sFromAft = sc.faceSFromAft;
g.minimumAngleDeg = minimumAngle;
g.minimumAngleThresholdDeg = threshold; g.poorAngle = minimumAngle < threshold;
g.discreteCurvatureMean = tmp.curvatureMean;
g.discreteCurvatureMax = tmp.curvatureMax;
g.discreteRadiusProxy = radius;
g.curvatureDefinition = 'Adjacent-normal angle / shared-edge length; diagnostic only';
g.featureFace = tmp.featureFace; g.topologyIssue = topologyIssue;
g.needsReview = g.poorAngle | g.featureFace | topologyIssue;
[g.summary.minimumAngleDeg,g.summary.worstFaceID] = min(minimumAngle);
g.summary.numberFaces = nf; g.summary.poorAngleFaces = nnz(g.poorAngle);
g.summary.undefinedFrames = nnz(~frameValid);
g.summary.boundaryEdges = nnz(counts == 1);
g.summary.nonManifoldEdges = nnz(counts > 2);
mesh.localGeometry = g;
verbose = true;
if isfield(cfg,'output') && isfield(cfg.output,'verbose'), verbose = cfg.output.verbose; end
if verbose
    fprintf('\nLocal geometry completed: %d faces, %d small-angle faces.\n',nf,nnz(g.poorAngle));
    fprintf('Minimum angle: %.4f deg (face %d); undefined axial frames: %d.\n', ...
        g.summary.minimumAngleDeg,g.summary.worstFaceID,nnz(~frameValid));
    fprintf('Curvature radius proxy is diagnostic only, not a heating-model radius.\n');
end
end
