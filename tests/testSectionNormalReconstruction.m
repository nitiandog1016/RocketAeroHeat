function testSectionNormalReconstruction()
%TESTSECTIONNORMALRECONSTRUCTION 解析锥面法向、折角/椭圆保护和真实判据回归。
root=fileparts(fileparts(mfilename('fullpath')));addpath(root);
evalc('cfg=createCaseConfig();');cfg.output.verbose=false;
[m,meta]=buildCanonicalBody('coneCylinder',struct('resolution',48));
m=prepare(m,cfg);old=m;frame=reconstructSectionNormals(m,cfg);
assert(isequaln(old,m),'法向重构不得修改计算网格');
cone=strcmp(meta.faceRegion,'cone');cylinder=strcmp(meta.faceRegion,'cylinder');
use=cone & frame.appliedFace;assert(nnz(use)>500);
radial=m.center(use,2:3);radial=radial./vecnorm(radial,2,2);
expected=[repmat(sind(15),nnz(use),1),cosd(15)*radial];
assert(max(vecnorm(frame.normal(use,:)-expected,2,2))<1e-10);
assert(max(abs(frame.normal(cylinder & frame.appliedFace,1)))<1e-10);
assert(max(frame.correctionAngleDeg(use))>3,'该网格应触发可测的拉链法向修正');
% 圆柱/锥面接缝两侧保留离散法向，不能用跨折角拟合平滑实际几何。
join=any(abs(reshape(m.V(m.F(:),1),size(m.F))-0.4)<1e-10,2);
assert(~any(frame.appliedFace(join)));
cfgFacet=cfg;cfgFacet.externalFlow.normalModel='facet';f=reconstructSectionNormals(m,cfgFacet);
assert(~any(f.appliedFace) && isequal(f.normal,m.localGeometry.normal));
bad=m;bad.V(:,2)=1.25*bad.V(:,2);bad=prepare(rmfield(bad,'param'),cfg);
elliptic=reconstructSectionNormals(bad,cfg);assert(~any(elliptic.appliedFace));

%% 2. 单独的平移/单位缩放不能改变几何判据或重构方向
shifted=old;shifted.V=2*old.V+[7 -3 5];
shifted.param.sectionX=2*old.param.sectionX+7;
shifted.param.sectionCenterYZ=repmat([-3 5],numel(shifted.param.sectionX),1);
shifted=prepare(shifted,cfg);f2=reconstructSectionNormals(shifted,cfg);
assert(isequal(frame.appliedFace,f2.appliedFace));
assert(max(vecnorm(frame.normal-f2.normal,2,2))<1e-9);

%% 3. 有攻角端到端：解除网格法向导致的伪分离，不改变转捩/分离阈值
cfg.flight.Mach=4;cfg.flight.altitude=20000;cfg.flight.alpha=20;
new=calcFlightStationState(m,cfg.flight,cfg);
cfgFacet.flight=cfg.flight;before=calcFlightStationState(m,cfgFacet.flight,cfgFacet);
assert(~any(before.heating.eligibleFace),'旧面法向路径应复现原缺陷');
assert(~any(new.heating.separationTrace.onsetFace));
fraction=sum(new.area(new.heating.eligibleFace))/sum(new.area);
assert(fraction>0.65);
assert(isequaln(new.pressure,before.pressure),'不得改变修正牛顿压力结果');
assert(isequal(new.heating.transition.criticalReynolds,before.heating.transition.criticalReynolds));
assert(all(isnan(new.heating.heatFlux(strcmp(meta.faceRegion,'base')))));
valid=new.boundaryLayerEdge.algebraicallyValid;
assert(max(abs(sum(new.boundaryLayerEdge.velocityVector(valid,:).*m.normal(valid,:),2)))<1e-9);
data=buildAerothermalViewerData(new);field=data.fields(strcmp({data.fields.key},'edgePressure'));
assert(any(contains(string(field.methodNames),'重构入射角')));
fprintf('Section normal reconstruction passed: analytic cone, corner/noncircular protection, invariant geometry and oblique heat area %.2f%%.\n',100*fraction);
end
function mesh=prepare(mesh,cfg)
% 复用实际生产几何流程，不用预置法向替代计算。
evalc('mesh=calcFaceGeometry(mesh);mesh=fixNormalDirection(mesh,cfg);mesh=buildSurfaceCoordinates(mesh,cfg);');
mesh=analyzeLocalGeometry(mesh,cfg);
end
