function frame=reconstructSectionNormals(mesh,cfg)
%RECONSTRUCTSECTIONNORMALS 在已验证的圆形直母线区域重构外流转角法向。
% 模块用途：抑制不同周向点数拉链连接造成的三角面法向振荡。
% 不修改网格、三角面面积/法向、牛顿压力积分或准流线的实际切平面。
% 只接受共轴圆截面、三个相邻环带斜率一致的区域；端盖、折角、
% 非圆截面、弯曲母线和过大修正保留面法向。可显式选择 facet 模式。
if nargin<2,cfg=struct;end
g=mesh.localGeometry;nf=size(mesh.F,1);mode='validatedSections';
if isfield(cfg,'externalFlow') && isfield(cfg.externalFlow,'normalModel')
    mode=validatestring(cfg.externalFlow.normalModel,{'validatedSections','facet'});
end
frame=struct('normal',g.normal,'appliedFace',false(nf,1), ...
    'correctionAngleDeg',zeros(nf,1),'mode',mode, ...
    'method','共轴圆截面检查 + 三环带直母线法向重构；不满足条件时用三角面法向', ...
    'scope','仅用于局部工程外流入射角；网格面积、受力法向和实际面内流向不变', ...
    'status','未应用截面重构','sectionValid',false(0,1));
if strcmp(mode,'facet'),frame.status='显式采用三角面法向';return;end
if ~isfield(mesh,'surfaceCoord'),return;end
sc=mesh.surfaceCoord;
if ~all(isfield(sc,{'sectionX','sectionNodeIDs','faceBandID'})),return;end
x=sc.sectionX(:);ns=numel(x);nb=ns-1;
if ns<4 || any(diff(x)<=0),return;end

%% 1. 独立拟合每圈的圆，检查最大径向误差与周向覆盖
radius=nan(ns,1);center=nan(ns,2);valid=false(ns,1);
for k=1:ns
    ids=sc.sectionNodeIDs{k};yz=double(mesh.V(ids,2:3));
    if numel(ids)<8,continue;end
    if max(abs(mesh.V(ids,1)-x(k)))>max(1e-10,1e-8*(max(x)-min(x))),continue;end
    origin=mean(yz,1);q=yz-origin;
    scale=max(vecnorm(q,2,2));if scale<=0,continue;end
    q=q/scale;A=[2*q,ones(numel(ids),1)];
    if rcond(A'*A)<1e-10,continue;end
    fit=A\sum(q.^2,2);rr=sqrt(sum(fit(1:2).^2)+fit(3));
    c=fit(1:2)';distance=vecnorm(q-c,2,2);
    if ~isreal(rr) || rr<=0,continue;end
    theta=sort(mod(atan2(q(:,2)-c(2),q(:,1)-c(1)),2*pi));
    gaps=diff([theta;theta(1)+2*pi]);
    valid(k)=max(abs(distance/rr-1))<=0.003 && max(gaps)<=pi/4+1e-6;
    radius(k)=rr*scale;center(k,:)=origin+c*scale;
end
frame.sectionValid=valid;
inclination=atan2d(-diff(radius),diff(x));
band=double(sc.faceBandID(:));

%% 2. 三个环带检验直母线；跨越肩部/尾锥折角时回退原始法向
for k=1:nb
    start=max(1,min(k-1,nb-2));bands=start:start+2;sections=start:start+3;
    if ~all(valid(sections)),continue;end
    angles=inclination(bands);
    if max(angles)-min(angles)>0.5,continue;end
    c=mean(center(sections,:),1);
    if max(vecnorm(center(sections,:)-c,2,2))>0.003*min(radius(sections)),continue;end
    % 最小二乘直线只在通过检查的局部母线上估计斜率。
    localX=x(sections)-mean(x(sections));
    slope=localX\(radius(sections)-mean(radius(sections)));
    ids=find(band==k & g.isSidewall & ~g.featureFace & ~g.topologyIssue & ...
        ~g.isNoseTipAdjacent & ~g.isAftTipAdjacent);
    if isempty(ids),continue;end
    radial=g.center(ids,2:3)-c;radial=radial./vecnorm(radial,2,2);
    normals=[repmat(-slope,numel(ids),1),radial];normals=normals./vecnorm(normals,2,2);
    correction=acosd(max(-1,min(1,sum(normals.*g.normal(ids,:),2))));
    keep=isfinite(correction) & correction<=15;ids=ids(keep);
    frame.normal(ids,:)=normals(keep,:);frame.appliedFace(ids)=true;
    frame.correctionAngleDeg(ids)=correction(keep);
end
frame.status=sprintf('已验证直母线重构 %d/%d 面元；其余保留三角面法向',nnz(frame.appliedFace),nf);
end
