function [mesh,meta]=buildNasaBoattailBody(angleDeg,resolution)
%BUILDNASABOATTAILBODY 生成NASA TN D-2761名义15度锥鼻、柱体和直尾锥闭合模型。
% 模块一：严格区分14英寸鼻底至尾底、10英寸基本柱体、4英寸可换尾段和1.18英寸收缩。
% 模块二：独立几何夹具按圆环三角化，不修改既有STL重网格器或生产模型。
% 模块三：输出解析面积/体积和逐面区域；未生成支杆、砂粒带，不能冒充实验安装全几何。
if nargin<2,resolution=48;end
validateattributes(angleDeg,{'numeric'},{'real','finite','scalar'});
assert(ismember(angleDeg,[6 12 18]),'NasaBoattailBody:Angle','当前只生成三组选定直尾锥。');
validateattributes(resolution,{'numeric'},{'scalar','integer','>=',16,'<=',192});
D=3.375*.0254;R=D/2;Lbt=1.18*.0254;Lc=(14-1.18)*.0254;
Ln=R/tand(15);Rb=R-Lbt*tand(angleDeg);ds=2*pi*R/resolution;
x=0;r=Rb;bands={};
segment(Lbt,R,'boattail');segment(Lbt+Lc,R,'cylinder');segment(Lbt+Lc+Ln,0,'cone');
nc=resolution;nr=numel(x)-1;theta=(0:nc-1)'*2*pi/nc;V=zeros(nr*nc+2,3);
for k=1:nr,ids=(k-1)*nc+(1:nc);V(ids,:)=[repmat(x(k),nc,1),r(k)*cos(theta),r(k)*sin(theta)];end
tip=nr*nc+1;base=tip+1;V(tip,:)=[x(end),0,0];
F=zeros(2*nc*nr,3);region=cell(size(F,1),1);j=0;
for k=1:nr-1
    a=(k-1)*nc+(1:nc);b=k*nc+(1:nc);next=[2:nc 1];
    idx=j+(1:2*nc);F(idx,:)=[a(:),b(:),a(next)';a(next)',b(:),b(next)'];
    region(idx)={bands{k}};j=j+2*nc;
end
a=(nr-1)*nc+(1:nc);idx=j+(1:nc);F(idx,:)=[a(:),repmat(tip,nc,1),a([2:nc 1])'];region(idx)={bands{nr}};j=j+nc;
a=1:nc;idx=j+(1:nc);F(idx,:)=[a(:),a([2:nc 1])',repmat(base,nc,1)];region(idx)={'base'};
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
if volume<0,F=F(:,[1 3 2]);volume=-volume;end
normal=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);area=vecnorm(normal,2,2)/2;
assert(all(area>0));normal=normal./(2*area);
mesh=struct('V',V,'F',F,'normal',normal,'area',area, ...
    'center',(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3, ...
    'localGeometry',struct('region',{region}));
exactVolume=pi*Lbt*(R^2+R*Rb+Rb^2)/3+pi*R^2*Lc+pi*R^2*Ln/3;
exactArea=pi*Rb^2+pi*(R+Rb)*hypot(Lbt,R-Rb)+2*pi*R*Lc+pi*R*hypot(R,Ln);
meta=struct('sourceReport','NASA TN D-2761 Fig.1 nominal nose section 1', ...
    'boattailAngleDeg',angleDeg,'diameterM',D,'boattailLengthM',Lbt,'cylinderLengthBeforeBoattailM',Lc, ...
    'noseBaseToBaseM',14*.0254,'noseConeLengthM',Ln,'baseRadiusM',Rb,'resolution',nc, ...
    'exactAreaM2',exactArea,'exactVolumeM3',exactVolume,'facetAreaM2',sum(area),'facetVolumeM3',volume, ...
    'faceRegion',{region},'nosePositiveX',true,'baseX_m',0,'stingIncluded',false,'gritBandGeometryIncluded',false, ...
    'geometryOnly',true,'physicalValidationComplete',false,'newProductionSolvedAreaM2',0);
    function segment(endpoint,radius,label)
        n=max(1,ceil(hypot(endpoint-x(end),radius-r(end))/ds));
        xx=linspace(x(end),endpoint,n+1);rr=linspace(r(end),radius,n+1);
        x=[x xx(2:end)];r=[r rr(2:end)];bands=[bands repmat({label},1,n)];
    end
end
