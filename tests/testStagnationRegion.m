function testStagnationRegion(out)
% 解析球体、尖锥、平端方柱和当前 STL：测试候选定位与半径保护。
addpath(fileparts(fileparts(mfilename('fullpath'))));
evalc('cfg=createCaseConfig();'); cfg.output.verbose=false;
flight=cfg.flight; flight.alpha=0; flight.beta=0;
fs=calcFreestream(flight);
% 规则纬线球体，X 为极轴，半径 0.2 m。
n=32; nt=16; R=0.2; V=[-R 0 0];
for k=1:nt-1
    angle=pi*k/nt; theta=(0:n-1)'*2*pi/n;
    V=[V;repmat(-R*cos(angle),n,1),R*sin(angle)*cos(theta),R*sin(angle)*sin(theta)]; %#ok<AGROW>
end
V=[V;R 0 0]; F=convhull(V); sphere=prepare(V,F,cfg);
sphere=findStagnationRegion(sphere,fs,cfg); s=sphere.stagnation;
assert(s.representativePosition(1)>0.17);
assert(s.numberComponents==1 && ~s.flags.weakNormalAlignment);
assert(abs(s.sphereFit.radius-R)<1e-10);
assert(s.hasPhysicalRadius && abs(s.physicalRadius-R)<1e-9);
assert(~s.readyForHeating && ~s.sphereFit.approvedForHeating);
known=cfg; known.stagnation.physicalNoseRadius=R;
verified=findStagnationRegion(sphere,fs,known);
assert(verified.stagnation.physicalRadius==R && ~verified.stagnation.readyForHeating);
% 攻角反号时，候选位置应向相反 Z 侧移动。
positions=zeros(2,3);
for k=1:2
    flight.alpha=(2*k-3)*30;
    tilted=findStagnationRegion(sphere,calcFreestream(flight),cfg);
    positions(k,:)=tilted.stagnation.representativePosition;
end
assert(positions(1,3)>0 && positions(2,3)<0);
% 平移不能改变候选面集合；物理半径也不随坐标平移变化。
translated=prepare(V+[20 -10 4],F,cfg);
translated=findStagnationRegion(translated,fs,cfg);
assert(isequal(translated.stagnation.candidateFace,s.candidateFace));
assert(abs(translated.stagnation.sphereFit.radius-R)<1e-8);
assert(abs(translated.stagnation.physicalRadius-R)<1e-9);
% 尖锥：即使迎流法向未接近 1，仍保留最前端邻接面。
theta=(0:n-1)'*2*pi/n;
coneV=[zeros(n,1),cos(theta),sin(theta);1 0 0];
cone=prepare(coneV,convhull(coneV),cfg);
cone=findStagnationRegion(cone,fs,cfg);
assert(cone.stagnation.flags.weakNormalAlignment);
assert(any(cone.stagnation.candidateFace & cone.stagnation.leadingFace));
assert(isnan(cone.stagnation.physicalRadius));
% 密集尖锥即使采样充分，也不能被局部球拟合误认成球形钝鼻。
coneDenseV=[0.05 0 0];
for j=1:24
    depth=j*0.001; theta=(0:n-1)'*2*pi/n;
    coneDenseV=[coneDenseV;repmat(0.05-depth,n,1), ...
        0.5*depth*cos(theta),0.5*depth*sin(theta)]; %#ok<AGROW>
end
coneEstimate=estimateNoseRadius(struct('V',coneDenseV));
assert(~coneEstimate.available && isnan(coneEstimate.radius));
% 方柱平端面不可凭压力峰值认定为圆钝头。
boxV=[0 -1 -1;0 1 -1;0 1 1;0 -1 1;1 -1 -1;1 1 -1;1 1 1;1 -1 1];
box=prepare(boxV,convhull(boxV),cfg);
box=findStagnationRegion(box,fs,cfg);
assert(box.stagnation.flags.flatCapCandidate && ~box.stagnation.readyForHeating);
% 实际项目完整数值链路与显示字段。
evalc('raw=loadSTLMesh(cfg); raw=calcFaceGeometry(raw);');
[raw,~]=alignGeometryAxis(raw,'pca');
evalc('m=autoRemeshSurface(raw,cfg); m=fixNormalDirection(m,cfg); m=buildSurfaceCoordinates(m,cfg);');
m=analyzeLocalGeometry(m,cfg); m=modifiedNewtonianPressure(m,fs,cfg);
m=integrateAerodynamicLoads(m,cfg); before=m;
m=findStagnationRegion(m,fs,cfg);
assert(isequal(m.V,before.V) && isequal(m.F,before.F) && isequal(m.pressure,before.pressure));
assert(m.stagnation.hasPhysicalRadius);
assert(abs(m.stagnation.physicalRadius-0.05)<0.005);
assert(m.stagnation.numberComponents==1);
m=calcStagnationHeating(m,fs,cfg);
assert(m.stagnationHeating.available);
expected=cfg.heating.stagnationCoefficient* ...
    sqrt(fs.density/m.stagnation.physicalRadius)*fs.speed^3;
assert(abs(m.stagnationHeating.heatFlux-expected)<1e-9*expected);
cfg.pressureGUI.visible='off'; fig=pressureGUI(m,cfg);
cleanup=onCleanup(@()delete(fig)); %#ok<NASGU>
ui=fig.UserData; ui.selector.Value=4; ui.selector.ValueChangedFcn(ui.selector,[]);
assert(isequal(ui.surface.FaceVertexCData,double(m.stagnation.candidateFace)));
if nargin>0, drawnow; exportapp(fig,fullfile(out,'stagnation_gui.png')); end
fprintf('Stagnation tests passed: sphere, cone, flat cap, incidence, translation and GUI.\n');
fprintf('Project: candidates=%d, regions=%d, mu=%.6f, radiusSource=%s\n', ...
    nnz(m.stagnation.candidateFace),m.stagnation.numberComponents,m.stagnation.peakMu,m.stagnation.radiusSource);
end

function mesh=prepare(V,F,cfg)
% 解析测试只需面元几何；为分类构造与生产模块相同的最小字段。
mesh=struct('V',V,'F',F);
evalc('mesh=calcFaceGeometry(mesh); mesh=fixNormalDirection(mesh,cfg);');
nf=size(mesh.F,1); g=struct('normal',mesh.normal,'center',mesh.center,'area',mesh.area);
g.region=repmat({'sidewall'},nf,1);
% 显式索引后重排为 Nface x 3，避免 MATLAB 单下标形状差异。
x=reshape(mesh.V(mesh.F(:),1),size(mesh.F));
g.isAftCap=all(abs(x-min(V(:,1)))<1e-10,2);
g.isNoseCap=all(abs(x-max(V(:,1)))<1e-10,2);
for name={'isAftTipAdjacent','isNoseTipAdjacent','poorAngle','featureFace','topologyIssue'}
    g.(name{1})=false(nf,1);
end
mesh.localGeometry=g;
end
