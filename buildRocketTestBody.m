function [mesh,meta] = buildRocketTestBody(cfg)
%BUILDROCKETTESTBODY 生成带光滑突起的完整箭体几何夹具，单位 m。
% 模块用途：为几何、网格保真和区域适用性试验提供可重复的完整模型。
% 主轴 X、鼻尖 +X；箭体依次含底面、后体、收缩肩、柱体、锥段和球鼻。
% 突起是同一径向图上的光滑整流罩/纵脊，既不叠放闭合实体，也不代表任意尾翼。
% 本程序只生成几何真值，不求解流场，不增加已经验证的物理解面积。
%
% cfg 常用字段：resolution(默认48)、protrusionEnabled(默认true)、
% protrusionType('fairing'或'ridge')；protrusions可覆盖为结构体数组。
% 每个突起可设 type、height、length、centerFraction、azimuthDeg、halfWidthDeg。
% centerFraction 是从收缩肩到鼻部方向的柱段长度比例，不是流向坐标。
% mesh.V/F 可直接供既有几何链路使用；meta 含逐面区域与突起精确支持域标签。
if nargin<1,cfg=struct;end
assert(isstruct(cfg)&&isscalar(cfg),'RocketAeroHeat:RocketFixture:Config','cfg 必须是标量结构体。');
allowed={'bodyRadius','baseRadius','afterbodyLength','shoulderLength', ...
    'cylinderLength','noseRadius','coneAngleDeg','resolution', ...
    'protrusionEnabled','protrusionType','protrusions'};
assert(all(ismember(fieldnames(cfg),allowed)),'RocketAeroHeat:RocketFixture:UnknownOption','存在不支持的 cfg 字段。');

%% 1. 校验尺寸及突起参数，拒绝会跨越柱段接缝的配置
R=option('bodyRadius',0.12);Rb=option('baseRadius',0.08);
La=option('afterbodyLength',0.24);Ls=option('shoulderLength',0.067);
Lc=option('cylinderLength',0.8);Rn=option('noseRadius',0.025);
angle=option('coneAngleDeg',18);n=option('resolution',48);
values={R,Rb,La,Ls,Lc,Rn};
for k=1:numel(values),validateattributes(values{k},{'numeric'},{'scalar','real','finite','positive'});end
validateattributes(angle,{'numeric'},{'scalar','real','finite','>',0,'<',80});
validateattributes(n,{'numeric'},{'scalar','real','finite','integer','>=',16,'<=',512});
assert(Rb<R && Rn*cosd(angle)<R,'RocketAeroHeat:RocketFixture:Profile','后体必须比柱体细，且球锥切点必须位于柱体半径内。');
enabled=option('protrusionEnabled',true);
assert(islogical(enabled)&&isscalar(enabled),'RocketAeroHeat:RocketFixture:Enabled','protrusionEnabled 必须为逻辑标量。');
defaultType=validatestring(option('protrusionType','fairing'),{'fairing','ridge','none'});
features=option('protrusions',struct('type',defaultType));
assert(isstruct(features),'RocketAeroHeat:RocketFixture:Features','protrusions 必须为结构体数组。');
if ~enabled || strcmp(defaultType,'none'),features=struct([]);end
nf=numel(features);spec=repmat(struct('type','','height',0,'length',0, ...
    'centerFraction',0,'azimuthDeg',0,'halfWidthDeg',0,'centerX',0),nf,1);
featureAllowed={'type','height','length','centerFraction','azimuthDeg','halfWidthDeg'};
xShoulder=La+Ls;xCylinder=xShoulder+Lc;
for k=1:nf
    assert(all(ismember(fieldnames(features(k)),featureAllowed)), ...
        'RocketAeroHeat:RocketFixture:FeatureOption','存在不支持的突起字段。');
    type=validatestring(featureOption(k,'type',defaultType),{'fairing','ridge'});
    if strcmp(type,'ridge'),lengthDefault=0.65*Lc;widthDefault=10;
    else,lengthDefault=0.3*Lc;widthDefault=22;end
    f=struct('type',type,'height',featureOption(k,'height',0.16*R), ...
        'length',featureOption(k,'length',lengthDefault), ...
        'centerFraction',featureOption(k,'centerFraction',0.55), ...
        'azimuthDeg',featureOption(k,'azimuthDeg',35), ...
        'halfWidthDeg',featureOption(k,'halfWidthDeg',widthDefault),'centerX',0);
    for key={'height','length'},validateattributes(f.(key{1}),{'numeric'},{'scalar','real','finite','positive'});end
    validateattributes(f.centerFraction,{'numeric'},{'scalar','real','finite','>',0,'<',1});
    validateattributes(f.azimuthDeg,{'numeric'},{'scalar','real','finite'});
    validateattributes(f.halfWidthDeg,{'numeric'},{'scalar','real','finite','>',0,'<',90});
    f.centerX=xShoulder+Lc*f.centerFraction;f.azimuthDeg=mod(f.azimuthDeg,360);
    assert(f.centerX-f.length/2>xShoulder && f.centerX+f.length/2<xCylinder, ...
        'RocketAeroHeat:RocketFixture:FeatureSupport','突起支持域必须完全位于柱段内部。');
    spec(k)=f;
end

%% 2. 解析母线与突起边界采样，所有支持域边界均为显式网格接缝
ds=2*pi*R/n;x=0;r=Rb;band={};
lineTo(La,Rb,'afterbody');lineTo(xShoulder,R,'contraction');
pins=[xShoulder,xCylinder];theta=(0:n-1)'*2*pi/n;
for k=1:nf
    f=spec(k);np=max(8,2*ceil(f.length/(2*ds)));
    pins=[pins,linspace(f.centerX-f.length/2,f.centerX+f.length/2,np+1)]; %#ok<AGROW>
    if strcmp(f.type,'ridge'),pins=[pins,f.centerX-f.length/4,f.centerX+f.length/4];end %#ok<AGROW>
    aw=deg2rad(f.halfWidthDeg);nt=max(8,2*ceil(2*aw/(4*pi/n)));
    theta=[theta;mod(deg2rad(f.azimuthDeg)+linspace(-aw,aw,nt+1)',2*pi)]; %#ok<AGROW>
end
pins=sort(pins);xTol=128*eps(max(1,xCylinder));
pins=pins([true,diff(pins)>xTol]);
for k=2:numel(pins),lineTo(pins(k),R,'cylinder');end
rt=Rn*cosd(angle);H=(R-rt)/tand(angle);xt=xCylinder+H;xc=xt-Rn*sind(angle);
lineTo(xt,rt,'cone');
na=max(4,ceil(Rn*(pi/2-deg2rad(angle))/ds));
t=linspace(pi/2-deg2rad(angle),0,na+1);
x=[x,xc+Rn*cos(t(2:end))];r=[r,Rn*sin(t(2:end))];
band=[band,repmat({'sphericalNose'},1,na)];
% 合并舍入到同一角度的节点，避免周期接缝附近产生零面积细片。
theta=sort(theta);tol=128*eps(2*pi);theta=theta([true;diff(theta)>tol]);
if 2*pi-theta(end)+theta(1)<tol,theta(end)=[];end
nc=numel(theta);ns=numel(x)-1;nv=ns*nc+2;
V=zeros(nv,3);baseRadius=zeros(nv,1);offset=zeros(nv,1);
nodeFeatureHeight=zeros(nv,nf);rings=cell(ns,1);
for k=1:ns
    ids=(k-1)*nc+(1:nc);rings{k}=ids;
    dh=featureHeights(repmat(x(k),nc,1),theta);rr=r(k)+sum(dh,2);
    V(ids,:)=[repmat(x(k),nc,1),rr.*cos(theta),rr.*sin(theta)];
    baseRadius(ids)=r(k);offset(ids)=sum(dh,2);nodeFeatureHeight(ids,:)=dh;
end
tip=nv-1;cap=nv;V(tip,:)=[x(end),0,0];V(cap,:)=[0,0,0];

%% 3. 一致有向的侧壁、鼻尖和底面三角化；逐面标签由参数支持域给出
nSide=2*nc*(ns-1)+nc;F=zeros(nSide+nc,3);
faceBaseRegion=cell(size(F,1),1);faceBand=zeros(size(F,1),1);
faceParamX=zeros(size(F,1),1);faceParamTheta=zeros(size(F,1),1);
next=[2:nc,1];thetaNext=theta(next);thetaNext(end)=theta(1)+2*pi;
cursor=0;
for k=1:ns-1
    a=rings{k};b=rings{k+1};ids=cursor+(1:2*nc);cursor=cursor+2*nc;
    F(ids(1:nc),:)=[a(:),a(next)',b(:)];
    F(ids(nc+1:end),:)=[a(next)',b(next)',b(:)];
    faceBaseRegion(ids)={band{k}};faceBand(ids)=k;
    faceParamX(ids(1:nc))=(2*x(k)+x(k+1))/3;
    faceParamX(ids(nc+1:end))=(x(k)+2*x(k+1))/3;
    faceParamTheta(ids(1:nc))=(2*theta+thetaNext)/3;
    faceParamTheta(ids(nc+1:end))=(theta+2*thetaNext)/3;
end
a=rings{end};ids=cursor+(1:nc);F(ids,:)=[a(:),a(next)',repmat(tip,nc,1)];
faceBaseRegion(ids)={band{end}};faceBand(ids)=ns;
faceParamX(ids)=(2*x(end-1)+x(end))/3;faceParamTheta(ids)=(theta+thetaNext)/2;
ids=nSide+(1:nc);a=rings{1};F(ids,:)=[repmat(cap,nc,1),a(next)',a(:)];
faceBaseRegion(ids)={'base'};faceBand(ids)=-1;faceParamTheta(ids)=(theta+thetaNext)/2;
faceContribution=featureHeights(faceParamX,faceParamTheta);
faceFeatureMask=faceContribution>0;faceRegion=faceBaseRegion;
for k=1:size(F,1)
    active=find(faceFeatureMask(k,:));
    if numel(active)>1,faceRegion{k}='protrusionOverlap';
    elseif numel(active)==1,faceRegion{k}=['protrusion_',spec(active).type];end
end
normalArea=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
area=vecnorm(normalArea,2,2)/2;normal=normalArea./(2*area);
center=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
assert(all(area>0)&&volume>0,'RocketAeroHeat:RocketFixture:Degenerate','几何出现退化或错误朝向。');
mesh=struct('V',V,'F',F,'center',center,'normal',normal,'area',area);

%% 4. 与既有程序兼容的原生截面信息，附带独立的逐节点/逐面几何真值
p=struct('sectionX',x(1:ns)','sectionCenterYZ',zeros(ns,2), ...
    'sectionNodeIDs',{rings},'nodeSectionID',zeros(nv,1),'nodeCircumID',zeros(nv,1), ...
    'nodeRingFraction',nan(nv,1),'sectionEquivalentRadius',zeros(ns,1), ...
    'sectionPerimeter',zeros(ns,1),'nodeIsCap',false(nv,1), ...
    'aftCapNodeID',cap,'noseCapNodeID',tip,'faceBandID',faceBand);
for k=1:ns
    ids=rings{k};p.nodeSectionID(ids)=k;p.nodeCircumID(ids)=(1:nc)';
    yz=V(ids,2:3);edge=vecnorm(yz(next,:)-yz,2,2);per=sum(edge);
    p.nodeRingFraction(ids)=[0;cumsum(edge(1:end-1))]/per;
    p.sectionEquivalentRadius(k)=sqrt(polyarea(yz(:,1),yz(:,2))/pi);p.sectionPerimeter(k)=per;
end
p.nodeIsCap([tip,cap])=true;mesh.param=p;
exactBaselineArea=pi*Rb^2+2*pi*Rb*La+pi*(Rb+R)*hypot(Ls,R-Rb)+ ...
    2*pi*R*Lc+pi*(R+rt)*hypot(H,R-rt)+2*pi*Rn^2*(1-sind(angle));
h=Rn*(1-sind(angle));exactBaselineVolume=pi*Rb^2*La+ ...
    pi*Ls*(Rb^2+Rb*R+R^2)/3+pi*R^2*Lc+ ...
    pi*H*(R^2+R*rt+rt^2)/3+pi*h^2*(Rn-h/3);
exactVolume=exactBaselineVolume;hasOverlap=false;
for k=1:nf
    f=spec(k);aw=deg2rad(f.halfWidthDeg);
    if strcmp(f.type,'fairing'),ix1=f.length/2;ix2=3*f.length/8;
    else,ix1=3*f.length/4;ix2=11*f.length/16;end
    exactVolume=exactVolume+R*f.height*ix1*aw+0.5*f.height^2*ix2*(3*aw/4);
    for j=1:k-1
        angularDistance=abs(atan2(sind(f.azimuthDeg-spec(j).azimuthDeg),cosd(f.azimuthDeg-spec(j).azimuthDeg)));
        if abs(f.centerX-spec(j).centerX)<(f.length+spec(j).length)/2 && ...
                angularDistance<deg2rad(f.halfWidthDeg+spec(j).halfWidthDeg),hasOverlap=true;end
    end
end
if hasOverlap,exactVolume=NaN;end % 不漏掉重叠突起的交叉体积项来伪称解析真值。
edges=[F(:,[1,2]);F(:,[2,3]);F(:,[3,1])];
[~,~,eid]=unique(sort(edges,2),'rows');counts=accumarray(eid,1);
directions=accumarray(eid,sign(edges(:,2)-edges(:,1)));
adj=sparse([F(:,1);F(:,2);F(:,3)],[F(:,2);F(:,3);F(:,1)],true,nv,nv);
component=conncomp(graph(adj+adj'));
audit=struct('boundaryEdges',sum(counts==1),'nonManifoldEdges',sum(counts>2), ...
    'inconsistentEdges',sum(directions~=0),'eulerCharacteristic',nv-max(eid)+size(F,1), ...
    'connectedComponents',max(component),'minimumTriangleArea',min(area), ...
    'surfaceArea',sum(area),'enclosedVolume',volume,'featureFaceCount',sum(faceFeatureMask,1), ...
    'featurePeakHeight',max(nodeFeatureHeight,[],1));
meta=struct('kind','parametricRocketRadialGraph','coordinateConvention','SI; nose +X; flow -X', ...
    'geometryOnly',true,'physicalValidation',false,'resolution',n,'bodyRadius',R, ...
    'baseRadius',Rb,'noseRadius',Rn,'coneAngleDeg',angle,'shoulderAngleDeg',atand((R-Rb)/Ls), ...
    'cylinderX',[xShoulder,xCylinder],'totalLength',x(end),'features',spec, ...
    'faceRegion',{faceRegion},'faceBaseRegion',{faceBaseRegion},'faceFeatureMask',faceFeatureMask, ...
    'faceParamX',faceParamX,'faceParamTheta',mod(faceParamTheta,2*pi), ...
    'nodeBaselineRadius',baseRadius,'nodeRadialOffset',offset,'nodeFeatureHeight',nodeFeatureHeight, ...
    'exactBaselineArea',exactBaselineArea,'exactBaselineVolume',exactBaselineVolume, ...
    'exactVolume',exactVolume,'featuresOverlap',hasOverlap,'exactArea',NaN, ...
    'geometryAudit',audit,'supportDescription','single watertight radial graph; smooth fairings/ridges on cylinder; no independent fins or undercuts');
if nf==0,meta.exactArea=exactBaselineArea;end

    function v=option(key,fallback)
        v=fallback;if isfield(cfg,key),v=cfg.(key);end
    end
    function v=featureOption(k,key,fallback)
        v=fallback;if isfield(features(k),key)&&~isempty(features(k).(key)),v=features(k).(key);end
    end
    function lineTo(xEnd,rEnd,label)
        count=max(1,ceil(hypot(xEnd-x(end),rEnd-r(end))/ds));
        xx=linspace(x(end),xEnd,count+1);rr=linspace(r(end),rEnd,count+1);
        x=[x,xx(2:end)];r=[r,rr(2:end)];band=[band,repmat({label},1,count)];
    end
    function dh=featureHeights(xx,tt)
        dh=zeros(numel(xx),nf);
        for j=1:nf
            f=spec(j);xi=2*(xx-f.centerX)/f.length;
            da=atan2(sin(tt-deg2rad(f.azimuthDeg)),cos(tt-deg2rad(f.azimuthDeg)));
            eta=da/deg2rad(f.halfWidthDeg);inside=abs(xi)<1 & abs(eta)<1;
            if strcmp(f.type,'fairing'),fx=0.5*(1+cos(pi*xi(inside)));
            else,q=abs(xi(inside));fx=ones(size(q));ramp=q>0.5;fx(ramp)=0.5*(1+cos(2*pi*(q(ramp)-0.5)));end
            dh(inside,j)=f.height*fx.*(0.5*(1+cos(pi*eta(inside))));
        end
    end
end
