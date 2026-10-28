function domain=buildRocketTriangularShellDomain(V,F,options)
%BUILDROCKETTRIANGULARSHELLDOMAIN 从原闭合三角面构造完整三维近壁层/外域控制体。
% 模块一：保留原V/F，寻找轴线上严格位于全部有向面内侧的核点，拒绝无核几何。
% 模块二：近壁按角度加权顶点法向外伸，外域围绕最后一层整体相似放大。
% 模块三：每个原三角形跨层成八面三角化控制体，共享面唯一，不采用二维流场复制。
% 模块四：由真实有向三角边界积分体积/质心，并检查每个CV的向量及九矩闭合。
% 这里只构造数值域；网格可用不等于流场已收敛、近壁充分或新增物理解面积。
if nargin<3,options=struct;end
defaults=struct('normalLayerHeightsM',[0 20e-6 80e-6 320e-6], ...
    'outerScaleFactors',[1.0002 1.001 1.005 1.02 1.08 1.32 2.28 6.12], ...
    'axisCenterYZ',[0 0],'observer',[]);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))), ...
    'RocketShell:Options','不支持未知网格选项。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
V=double(V);F=double(F);validateattributes(V,{'double'},{'finite','real','ncols',3,'nonempty'});
validateattributes(F,{'double'},{'finite','integer','ncols',3,'>=',1,'<=',size(V,1)});
h=double(options.normalLayerHeightsM(:));assert(numel(h)>=2&&h(1)==0&&all(isfinite(h))&&all(diff(h)>0), ...
    'RocketShell:Layers','层高度必须从0严格递增。');
validateattributes(options.axisCenterYZ,{'numeric'},{'real','finite','numel',2});
scales=double(options.outerScaleFactors(:));assert(~isempty(scales)&&all(isfinite(scales))&&scales(1)>1&&all(diff(scales)>0), ...
    'RocketShell:OuterScales','外域相似尺度须从大于1开始严格递增。');
assert(isempty(options.observer)||isa(options.observer,'function_handle'),'RocketShell:Observer','观察入口须为空或函数句柄。');
nf=size(F,1);nv=size(V,1);nl=numel(h)-1+numel(scales);nc=nf*nl;
a=V(F(:,1),:);b=V(F(:,2),:);c=V(F(:,3),:);av=cross(b-a,c-a,2)/2;area=vecnorm(av,2,2);
assert(all(area>0),'RocketShell:Degenerate','原面不能退化。');n=av./area;wallCenter=(a+b+c)/3;
edges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];source=repmat((1:nf)',3,1);
[edgeNodes,~,id]=unique(sort(edges,2),'rows');ne=size(edgeNodes,1);
assert(all(accumarray(id,1)==2)&&all(accumarray(id,sign(edges(:,2)-edges(:,1)))==0), ...
    'RocketShell:Topology','原表面必须闭合且方向一致。');
first=accumarray(id,(1:3*nf)',[],@min);last=accumarray(id,(1:3*nf)',[],@max);
left=source(first);right=source(last);edgeSign=sign(edges(first,2)-edges(first,1));

%% 模块一：核点与原面正四面体；不移动/光顺/旋转源面
yz=double(options.axisCenterYZ(:)');rhs=sum(n.*wallCenter,2)-n(:,2:3)*yz';
fun=@(x)-min(rhs-n(:,1)*x);
cx=fminbnd(fun,min(V(:,1)),max(V(:,1)),optimset('TolX',1e-12,'Display','off'));
kernel=[cx yz];height=rhs-n(:,1)*cx;
assert(all(height>0),'RocketShell:NoKernel','当前原几何在指定轴线上没有严格核点；不能强行外伸。');
bodyVolume=sum(area.*height)/3;assert(bodyVolume>0);
vertexNormal=zeros(nv,3);
for k=1:3
    p=V(F(:,k),:);u=V(F(:,mod(k,3)+1),:)-p;w=V(F(:,mod(k+1,3)+1),:)-p;
    angle=atan2(vecnorm(cross(u,w,2),2,2),sum(u.*w,2));
    for j=1:3,vertexNormal(:,j)=vertexNormal(:,j)+accumarray(F(:,k),angle.*n(:,j),[nv 1]);end
end
vertexNormal=vertexNormal./vecnorm(vertexNormal,2,2);
assert(all(isfinite(vertexNormal),'all'));
alignment=zeros(nf,3);for k=1:3,alignment(:,k)=sum(vertexNormal(F(:,k),:).*n,2);end
assert(all(alignment>0,'all'),'RocketShell:Normal','原网格角点法向不能向任一相邻面内侧外伸。');
nodes=zeros(nv,3,nl+1);nodes(:,:,1)=V;
for j=2:numel(h),nodes(:,:,j)=V+h(j)*vertexNormal;end
collar=nodes(:,:,numel(h));collarAV=cross(collar(F(:,2),:)-collar(F(:,1),:),collar(F(:,3),:)-collar(F(:,1),:),2)/2;
collarC=(collar(F(:,1),:)+collar(F(:,2),:)+collar(F(:,3),:))/3;
assert(all(sum((collarC-kernel).*collarAV,2)>0),'RocketShell:CollarKernel', ...
    '近壁最后一层不再相对原核点严格有向，不能构造外域。');
for j=1:numel(scales),nodes(:,:,numel(h)+j)=kernel+scales(j)*(collar-kernel);end

%% 模块二：每个共享四边形沿同一对角线分两三角面，内/外帽均保留三角身份
count=nl*(nf+2*ne);faceAV=zeros(count,3);faceC=zeros(count,3);minus=zeros(count,1);plus=zeros(count,1);
vertexIDs=zeros(count,3);kind=zeros(count,1);references=zeros(nc,3);cursor=0;
for j=1:nl
    inner=nodes(:,:,j);outer=nodes(:,:,j+1);base=(j-1)*nf;
    references(base+(1:nf),:)=(inner(F(:,1),:)+inner(F(:,2),:)+inner(F(:,3),:)+ ...
        outer(F(:,1),:)+outer(F(:,2),:)+outer(F(:,3),:))/6;
    ii=cursor+(1:nf);cursor=cursor+nf;
    store(ii,outer(F(:,1),:),outer(F(:,2),:),outer(F(:,3),:),ones(nf,1));
    minus(ii)=base+(1:nf);if j<nl,plus(ii)=j*nf+(1:nf);end
    vertexIDs(ii,:)=F+j*nv;kind(ii)=1+(j==nl);
    e0=inner(edgeNodes(:,1),:);e1=inner(edgeNodes(:,2),:);e2=outer(edgeNodes(:,2),:);e3=outer(edgeNodes(:,1),:);
    ii=cursor+(1:ne);cursor=cursor+ne;store(ii,e0,e1,e2,edgeSign);
    minus(ii)=base+left;plus(ii)=base+right;vertexIDs(ii,:)=[edgeNodes(:,1)+(j-1)*nv,edgeNodes(:,2)+(j-1)*nv,edgeNodes(:,2)+j*nv];kind(ii)=3;
    reverse=ii(edgeSign<0);vertexIDs(reverse,:)=vertexIDs(reverse,[1 3 2]);
    ii=cursor+(1:ne);cursor=cursor+ne;store(ii,e0,e2,e3,edgeSign);
    minus(ii)=base+left;plus(ii)=base+right;vertexIDs(ii,:)=[edgeNodes(:,1)+(j-1)*nv,edgeNodes(:,2)+j*nv,edgeNodes(:,1)+j*nv];kind(ii)=3;
    reverse=ii(edgeSign<0);vertexIDs(reverse,:)=vertexIDs(reverse,[1 3 2]);
end
faceArea=vecnorm(faceAV,2,2);assert(all(faceArea>0),'RocketShell:Face','层间共享面退化。');
faces=struct('area',faceArea,'normal',faceAV./faceArea,'centroid',faceC,'cellMinus',minus,'cellPlus',plus, ...
    'axis',2*ones(count,1),'areaVector',faceAV,'vertexIDs',vertexIDs,'kind',kind);
walls=struct('area',area,'normal',-n,'areaVector',-av,'centroid',wallCenter,'ownerCell',(1:nf)', ...
    'originalFaceID',(1:nf)','vertexIDs',F(:,[1 3 2]));

%% 模块三：局部参考点积分，避免薄层体积由大数相减获得
active=plus>0;owner=[minus;plus(active);walls.ownerCell];
boundaryC=[faceC;faceC(active,:);wallCenter];boundaryAV=[faceAV;-faceAV(active,:);-av];
relative=boundaryC-references(owner,:);tetraVolume=sum(relative.*boundaryAV,2)/3;
volume=accumarray(owner,tetraVolume,[nc 1]);
if ~isempty(options.observer)
    options.observer(struct('kind','volume-candidate','volume',volume,'references',references, ...
        'openFaces',faces,'walls',walls,'nodeLayers',nodes,'sourceGeometry',struct('V',V,'F',F), ...
        'kernelPoint',kernel,'normalLayerHeightsM',h,'outerScaleFactors',scales,'sourceFaces',nf));
end
assert(all(isfinite(volume)&volume>0),'RocketShell:Volume','出现非正真实CV体积；不能裁成正值。');
moment=zeros(nc,3);for j=1:3,moment(:,j)=accumarray(owner,tetraVolume.*(.75*relative(:,j)),[nc 1]);end
center=references+moment./volume;delta=boundaryC-center(owner,:);
signedHeight=sum(delta.*boundaryAV,2)./vecnorm(boundaryAV,2,2);
assert(all(signedHeight>0),'RocketShell:Centroid','某CV质心越过其边界面，不接受当前层构造。');
vector=zeros(nc,3);for j=1:3,vector(:,j)=accumarray(owner,boundaryAV(:,j),[nc 1]);end
tensorError=zeros(nc,1);
for j=1:3
    for k=1:3
        value=accumarray(owner,delta(:,j).*boundaryAV(:,k),[nc 1])-double(j==k)*volume;
        tensorError=max(tensorError,abs(value)./volume);
    end
end
areaSum=accumarray(owner,vecnorm(boundaryAV,2,2),[nc 1]);vectorError=vecnorm(vector,2,2)./areaSum;
assert(max(vectorError)<2e-12&&max(tensorError)<1e-8,'RocketShell:Closure','三维CV向量/九矩闭合未通过。');
wallDistance=sum((wallCenter-center(1:nf,:)).*(-n),2);
domain=struct('kind','original-triangle-body-fitted-three-dimensional-shell', ...
    'inputGeometry',struct('V',V,'F',F),'cells',struct('volume',volume,'centroid',center), ...
    'openFaces',faces,'walls',walls,'nodeLayers',nodes,'options',options, ...
    'kernelPoint',kernel,'kernelFaceHeightM',height,'vertexNormal',vertexNormal, ...
    'quality',struct('cellCount',nc,'layerCount',nl,'sourceFaces',nf, ...
        'minimumVolumeM3',min(volume),'maximumVectorClosureRelative',max(vectorError), ...
        'maximumTensorClosureOverVolume',max(tensorError),'minimumCVBoundaryDistanceM',min(signedHeight), ...
        'wallNormalDistanceRangeM',[min(wallDistance),max(wallDistance)],'bodyVolumeM3',bodyVolume, ...
        'physicalSolution',false,'eligibleForSurfaceCoverage',false,'newProductionSolvedAreaM2',0));
    function store(ids,p,q,r,orientation)
        faceAV(ids,:)=cross(q-p,r-p,2)/2.*orientation;
        faceC(ids,:)=(p+q+r)/3;
    end
end
