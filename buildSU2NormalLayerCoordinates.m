function [fv,quality]=buildSU2NormalLayerCoordinates(base,domain,options)
%BUILDSU2NORMALLAYERCOORDINATES 保留原母线，建立平行偏置近壁层并连接同一外域。
% 模块一：相邻真实墙段的偏置直线求交；角点使用斜接，不抹平肩部/底缘。
% 模块二：近壁首层为真实段法向距离，外部层从近壁壳连接原远场。
% 模块三：全部四边形必须严格凸，失败即拒绝，不裁节点或删除难区。
% 这是几何构造，不宣称y+、边界层厚度或外域独立性已合格。
if nargin<3,options=struct;end
defaults=struct('normalLayers',24,'normalGrowth',1.25);
assert(all(ismember(fieldnames(options),fieldnames(defaults))),'未知近壁层选项。');
for key=fieldnames(defaults)',k=key{1};if ~isfield(options,k),options.(k)=defaults.(k);end,end
layers=size(base.X,2)-1;count=options.normalLayers;growth=options.normalGrowth;
validateattributes(count,{'numeric'},{'scalar','integer','>=',1,'<',layers});
validateattributes(growth,{'numeric'},{'scalar','finite','>=',1,'<=',1.4});
first=domain.options.firstRadialSpacing;wall=domain.wallPolylineXR;
outer=[base.X(:,end),base.R(:,end)];edge=diff(wall);normal=[edge(:,2),-edge(:,1)]./vecnorm(edge,2,2);
direction=zeros(size(wall));
direction(1,:)=[1/normal(1,1),0];direction(end,:)=[1/normal(end,1),0];
den=1+dot(normal(1:end-1,:),normal(2:end,:),2);
assert(all(den>.1),'RocketAeroHeat:SU2:NormalCorner','墙面折角过强，偏置层不能在此构造。');
direction(2:end-1,:)=(normal(1:end-1,:)+normal(2:end,:))./den;
distances=[0 cumsum(first*growth.^(0:count-1))];
X=zeros(size(base.X));R=X;
for j=1:count+1
    p=wall+distances(j)*direction;X(:,j)=p(:,1);R(:,j)=p(:,2);
end
assert(all(R(:)>=0),'RocketAeroHeat:SU2:NormalAxis','偏置层跨过轴。');
remaining=layers-count;lastStep=first*growth^(count-1);
shell=[X(:,count+1),R(:,count+1)];vector=outer-shell;span=vecnorm(vector,2,2);outerGrowth=zeros(size(span));
for i=1:size(wall,1)
    step=lastStep*growth;assert(step*remaining<span(i),'近壁壳与外域之间空间不足。');
    lo=1;hi=2;while step*sum(hi.^(0:remaining-1))<span(i),hi=hi*2;end
    for k=1:80
        mid=(lo+hi)/2;
        if step*sum(mid.^(0:remaining-1))<span(i),lo=mid;else,hi=mid;end
    end
    q=(lo+hi)/2;outerGrowth(i)=q;fraction=cumsum(step*q.^(0:remaining-1))/span(i);fraction(end)=1;
    X(i,count+2:end)=shell(i,1)+vector(i,1)*fraction;
    R(i,count+2:end)=shell(i,2)+vector(i,2)*fraction;
end
X(:,1)=wall(:,1);R(:,1)=wall(:,2);X(:,end)=outer(:,1);R(:,end)=outer(:,2);
R([1 end],:)=0;
fv=buildAxisymmetricFVGeometry(X,R);
displace=[X(:,2)-X(:,1),R(:,2)-R(:,1)];
left=dot(displace(1:end-1,:),normal,2);right=dot(displace(2:end,:),normal,2);
error=max(abs([left;right]/first-1));assert(error<2e-9,'近壁法向距离未命中指定首层。');
midDisplace=(displace(1:end-1,:)+displace(2:end,:))/2;
alignment=asind(min(1,dot(midDisplace,normal,2)./vecnorm(midDisplace,2,2)));
old=[base.X(:,2)-base.X(:,1),base.R(:,2)-base.R(:,1)];oldMid=(old(1:end-1,:)+old(2:end,:))/2;
oldAlignment=asind(min(1,dot(oldMid,normal,2)./vecnorm(oldMid,2,2)));
quality=struct('method','piecewise-parallel-offset-normal-layers-with-miter-corners', ...
    'normalLayers',count,'normalGrowth',growth,'normalShellThickness',distances(end), ...
    'firstNormalSpacing',first,'maximumFirstNormalSpacingRelativeError',error, ...
    'wallAlignmentDeg',alignment,'minimumWallAlignmentDeg',min(alignment), ...
    'baselineWallAlignmentDeg',oldAlignment,'baselineMinimumWallAlignmentDeg',min(oldAlignment), ...
    'outerGrowthRatio',outerGrowth,'farfieldUnchanged',true,'originalWallUnchanged',true, ...
    'physicsQualified',false,'nearWallResolutionValidated',false);
end
