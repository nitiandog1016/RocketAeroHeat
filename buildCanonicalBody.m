function [mesh,meta] = buildCanonicalBody(kind,options)
%BUILDCANONICALBODY 建立经典轴对称验证模型，单位 m，鼻尖位于 +X。
% 模块用途：直接生成闭合三角面，供既有几何/外流/热流链路校核。
% 不修改 STL 重构算法；区域真值由解析母线给出，不依赖自动分类器。
% kind: sphere、cone、coneCylinder、hemisphereCylinder、bluntedCone。
% options: radius、cylinderLength、noseRadius、coneAngleDeg、resolution。
if nargin<2,options=struct;end
R=option('radius',0.1); L=option('cylinderLength',0.4);
Rn=option('noseRadius',0.04); angle=option('coneAngleDeg',15);
n=option('resolution',32);
validateattributes(R,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(L,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(n,{'numeric'},{'scalar','integer','>=',12});
validateattributes(Rn,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(angle,{'numeric'},{'scalar','real','finite','>',0,'<',90});
kind=validatestring(kind,{'sphere','cone','coneCylinder','hemisphereCylinder','bluntedCone'});
ds=2*pi*R/n; x=[]; r=[]; band={}; refArea=0; refVolume=0;

%% 1. 解析母线及精确面积/体积，保留接缝节点
if strcmp(kind,'sphere')
    arc(0,R,pi,pi/2,'aftSphere'); arc(0,R,pi/2,0,'sphericalNose');
    refArea=4*pi*R^2; refVolume=4*pi*R^3/3;
else
    if strcmp(kind,'cone'), L=0; end
    x=0;r=R;
    if L>0,lineTo(L,R,'cylinder');end
    refArea=2*pi*R*L+pi*R^2; refVolume=pi*R^2*L;
    switch kind
        case {'cone','coneCylinder'}
            H=R/tand(angle);lineTo(L+H,0,'cone');
            refArea=refArea+pi*R*hypot(R,H);
            refVolume=refVolume+pi*R^2*H/3;
        case 'hemisphereCylinder'
            arc(L,R,pi/2,0,'sphericalNose');
            refArea=refArea+2*pi*R^2;refVolume=refVolume+2*pi*R^3/3;
        case 'bluntedCone'
            rt=Rn*cosd(angle);
            assert(rt<R,'RocketAeroHeat:Canonical:Radius','球头切点半径必须小于柱体半径。');
            H=(R-rt)/tand(angle);xt=L+H;xc=xt-Rn*sind(angle);
            lineTo(xt,rt,'cone');arc(xc,Rn,pi/2-deg2rad(angle),0,'sphericalNose');
            h=Rn*(1-sind(angle));
            refArea=refArea+pi*(R+rt)*hypot(H,R-rt)+2*pi*Rn*h;
            refVolume=refVolume+pi*H*(R^2+R*rt+rt^2)/3+pi*h^2*(Rn-h/3);
    end
end

%% 2. 周向点数随半径变化，避免极点附近仍保留整圈高密度节点
V=zeros(0,3);rings=cell(numel(x),1);
for k=1:numel(x)
    if abs(r(k))<1e-12*R
        ids=size(V,1)+1;V(ids,:)=[x(k),0,0];
    else
        count=max(8,round(n*r(k)/R));theta=(0:count-1)'*2*pi/count;
        ids=size(V,1)+(1:count);V(ids,:)=[repmat(x(k),count,1),r(k)*cos(theta),r(k)*sin(theta)];
    end
    rings{k}=ids;
end
F=zeros(0,3);regions=cell(0,1);
for k=1:numel(rings)-1
    first=size(F,1)+1;a=rings{k};b=rings{k+1};na=numel(a);nb=numel(b);
    if na==1
        F=[F;repmat(a,nb,1),b(:),b([2:end 1])']; %#ok<AGROW>
    elseif nb==1
        F=[F;a(:),repmat(b,na,1),a([2:end 1])']; %#ok<AGROW>
    else
        % 两环按归一化周向进度拉链连接，保持每条环边只使用一次。
        i=0;j=0;
        while i<na || j<nb
            ai=a(mod(i,na)+1);bj=b(mod(j,nb)+1);
            if i<na && (j==nb || (i+1)/na<=(j+1)/nb)
                F(end+1,:)=[ai,bj,a(mod(i+1,na)+1)];i=i+1; %#ok<AGROW>
            else
                F(end+1,:)=[ai,bj,b(mod(j+1,nb)+1)];j=j+1; %#ok<AGROW>
            end
        end
    end
    regions(first:size(F,1),1)={band{k}};
end
if ~strcmp(kind,'sphere')
    a=rings{1};center=size(V,1)+1;V(center,:)=[0 0 0];
    F=[F;a(:),a([2:end 1])',repmat(center,numel(a),1)];
    regions(end+1:end+numel(a),1)={'base'};
end

%% 3. 整体有向体积校正，输出独立几何真值
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
if volume<0,F=F(:,[1 3 2]);end
mesh=struct('V',V,'F',F);
% 原生截面信息用于轴向波系推进；仅有推断曲面坐标会使该模块跳过。
positive=find(cellfun(@numel,rings)>1);ns=numel(positive);
p=struct('sectionX',x(positive)', 'sectionCenterYZ',zeros(ns,2), ...
    'sectionNodeIDs',{rings(positive)},'nodeSectionID',zeros(size(V,1),1), ...
    'nodeCircumID',zeros(size(V,1),1),'nodeRingFraction',nan(size(V,1),1), ...
    'sectionEquivalentRadius',zeros(ns,1),'sectionPerimeter',zeros(ns,1));
for k=1:ns
    ids=rings{positive(k)};count=numel(ids);rr=r(positive(k));
    p.nodeSectionID(ids)=k;p.nodeCircumID(ids)=(1:count)';
    p.nodeRingFraction(ids)=(0:count-1)'/count;
    p.sectionEquivalentRadius(k)=rr*sqrt(count*sin(2*pi/count)/(2*pi));
    p.sectionPerimeter(k)=2*count*rr*sin(pi/count);
end
p.nodeIsCap=p.nodeSectionID==0;p.aftCapNodeID=find(p.nodeIsCap & V(:,1)==min(V(:,1)),1);
p.noseCapNodeID=find(p.nodeIsCap & V(:,1)==max(V(:,1)),1);
p.faceBandID=zeros(size(F,1),1);
for k=1:size(F,1)
    sections=p.nodeSectionID(F(k,:));nonzero=sections(sections>0);
    if any(sections==0)
        if max(nonzero)==1,p.faceBandID(k)=-1;else,p.faceBandID(k)=ns;end
    else
        p.faceBandID(k)=min(nonzero);
    end
end
mesh.param=p;
meta=struct('kind',kind,'radius',R,'cylinderLength',L,'noseRadius',NaN, ...
    'coneAngleDeg',angle,'resolution',n,'faceRegion',{regions}, ...
    'exactArea',refArea,'exactVolume',refVolume,'coordinateConvention','SI; nose +X; flow -X');
if ismember(kind,{'sphere','hemisphereCylinder'}),meta.noseRadius=R;
elseif strcmp(kind,'bluntedCone'),meta.noseRadius=Rn;end

    function value=option(key,fallback)
        value=fallback;if isfield(options,key),value=options.(key);end
    end
    function lineTo(xend,rend,label)
        count=max(1,ceil(hypot(xend-x(end),rend-r(end))/ds));
        xx=linspace(x(end),xend,count+1);rr=linspace(r(end),rend,count+1);
        x=[x xx(2:end)];r=[r rr(2:end)];band=[band repmat({label},1,count)];
    end
    function arc(xc,rad,start,finish,label)
        count=max(2,ceil(rad*abs(finish-start)/ds));
        theta=linspace(start,finish,count+1);xx=xc+rad*cos(theta);rr=rad*sin(theta);
        if isempty(x),x=xx(1);r=rr(1);end
        x=[x xx(2:end)];r=[r rr(2:end)];band=[band repmat({label},1,count)];
    end
end
