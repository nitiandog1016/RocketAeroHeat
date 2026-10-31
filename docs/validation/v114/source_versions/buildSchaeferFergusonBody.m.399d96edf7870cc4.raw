function [mesh,meta]=buildSchaeferFergusonBody(options)
%BUILDSCHAEFERFERGUSONBODY 生成公开Mach5试验的完整锥—柱—10°扩张体。
% 模块一：按原图英寸尺寸建模；鼻锥角由直径/长度决定，不强行兼用舍入后的15°。
% 模块二：接缝共节点、鼻尖单节点、底面封闭；保留全部区域并给出解析几何账。
% 模块三：项目坐标鼻尖+X；报告流向距离s=L-X。底盖不代表实验支杆装置。
if nargin<1,options=struct;end
defaults=struct('circumferentialPoints',32,'axialStep',.008);
assert(all(ismember(fieldnames(options),fieldnames(defaults))),'不支持的经典模型参数。');
for k=fieldnames(defaults)',key=k{1};if ~isfield(options,key),options.(key)=defaults.(key);end,end
n=options.circumferentialPoints;ds=options.axialStep;
validateattributes(n,{'numeric'},{'scalar','integer','>=',16,'<=',256});
validateattributes(ds,{'numeric'},{'scalar','finite','positive','<=',.03});
inch=.0254;R=.9*inch;Ln=3.33*inch;Lc=8*inch;Lf=6.24*inch;L=Ln+Lc+Lf;
Rb=R+Lf*tand(10);x=0;r=Rb;bands={};
append(Lf,R,'flare');append(Lf+Lc,R,'cylinder');append(L,0,'cone');
nr=numel(x)-1;theta=(0:n-1)'*2*pi/n;V=zeros(nr*n+2,3);rings=cell(nr,1);
for j=1:nr
    ids=(j-1)*n+(1:n);rings{j}=ids;
    V(ids,:)=[repmat(x(j),n,1),r(j)*cos(theta),r(j)*sin(theta)];
end
tip=nr*n+1;cap=tip+1;V(tip,:)=[L 0 0];V(cap,:)=[0 0 0];
F=zeros(2*n*nr,3);region=cell(size(F,1),1);bandID=zeros(size(F,1),1);cursor=0;next=[2:n 1];
for j=1:nr-1
    a=rings{j};b=rings{j+1};ii=cursor+(1:2*n);cursor=cursor+2*n;
    F(ii,:)=[a(:) a(next)' b(:);a(next)' b(next)' b(:)];
    region(ii)={bands{j}};bandID(ii)=j;
end
a=rings{end};ii=cursor+(1:n);F(ii,:)=[a(:) a(next)' repmat(tip,n,1)];
region(ii)={'cone'};bandID(ii)=nr;
a=rings{1};ii=cursor+n+(1:n);F(ii,:)=[repmat(cap,n,1) a(next)' a(:)];region(ii)={'base'};bandID(ii)=-1;
av=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2)/2;area=vecnorm(av,2,2);
center=(V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))/3;
volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
e=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];[~,~,g]=unique(sort(e,2),'rows');
counts=accumarray(g,1);balance=accumarray(g,sign(e(:,2)-e(:,1)));
assert(all(area>0) && volume>0 && all(counts==2) && all(balance==0),'经典模型非闭合或方向错误。');
assert(size(V,1)-max(g)+size(F,1)==2);
exactVolume=pi*R^2*(Lc+Ln/3)+pi*Lf*(R^2+R*Rb+Rb^2)/3;
exactArea=pi*R*hypot(Ln,R)+2*pi*R*Lc+pi*(R+Rb)*hypot(Lf,Rb-R)+pi*Rb^2;
mesh=struct('V',V,'F',F,'center',center,'normal',av./area,'area',area, ...
    'param',struct('faceBandID',bandID));
meta=struct('kind','Schaefer-Ferguson-1961-cone-cylinder-10deg-flare', ...
    'source','https://ntrs.nasa.gov/citations/20150019680','figure',1, ...
    'coordinateConvention','SI; nose +X; zero-incidence velocity -X; source s=totalLength-X', ...
    'diameter',2*R,'noseLength',Ln,'cylinderLength',Lc,'flareLength',Lf,'totalLength',L, ...
    'flareAngleDeg',10,'nominalNoseAngleDeg',15,'actualNoseAngleDeg',atand(R/Ln), ...
    'dimensionRule','1.80in diameter, 3.33in cone, 8.00in cylinder, 6.24in flare; exact 10deg flare', ...
    'baseRadius',Rb,'faceRegion',{region},'faceDownstreamDistance',L-center(:,1), ...
    'exactSmoothArea',exactArea,'exactSmoothVolume',exactVolume, ...
    'facetedArea',sum(area),'facetedVolume',volume,'circumferentialPoints',n, ...
    'facetedVolumeRatio',sin(2*pi/n)/(2*pi/n),'options',options, ...
    'baseBoundaryAssumption','closed flat base; tunnel sting/support not represented; no base validation', ...
    'geometryOnly',true,'physicsQualified',false,'experimentalComparisonPerformed',false);
    function append(xEnd,rEnd,label)
        count=max(1,ceil((xEnd-x(end))/ds));xx=linspace(x(end),xEnd,count+1);rr=linspace(r(end),rEnd,count+1);
        x=[x xx(2:end)];r=[r rr(2:end)];bands=[bands repmat({label},1,count)];
    end
end
