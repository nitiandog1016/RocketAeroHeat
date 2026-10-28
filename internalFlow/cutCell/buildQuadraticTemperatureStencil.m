function stencil=buildQuadraticTemperatureStencil(cut,moments,boundary,options)
%BUILDQUADRATICTEMPERATURESTENCIL 真实FV温度均值的二次多项式重构模板草稿。
% 模块用途：p_i(r)=Tbar_i+g_i*r'+.5*H_i:(r*r'-Cov_i)。
% 未知量为3个梯度及6个对称Hessian分量；邻居均值约束保留各自Cov，
% 点值Dirichlet约束减去本格Cov。只沿正开放面扩展邻居，不跨实体补点。
% 局部长度及列范数缩放后用SVD，秩不足明确拒绝，不悄悄改回线性模型。
% 这是标量温度空间模板；一般Q经EOS所得primitive不自动等于温度FV均值。
if nargin<4,options=struct;end
if ~isfield(options,'maximumRing'),options.maximumRing=3;end
if ~isfield(options,'minimumRelativeSingularValue'),options.minimumRelativeSingularValue=1e-10;end
validateattributes(options.maximumRing,{'numeric'},{'integer','scalar','>=',1,'<=',8});
validateattributes(options.minimumRelativeSingularValue,{'numeric'},{'real','finite','scalar','positive','<',1});
volume=cut.cells.volume(:);centroid=cut.cells.centroid;n=numel(volume);fluid=volume>0;
assert(isfield(moments,'geometryIdentity')&& ...
    isequaln(moments.geometryIdentity.volume,volume)&&isequaln(moments.geometryIdentity.centroid,centroid)&& ...
    isequaln(moments.geometryIdentity.openFaces,cut.openFaces)&&isequaln(moments.geometryIdentity.walls,cut.walls), ...
    'CutQuadratic:MomentGeometryIdentity','二阶矩须由当前完整CV几何生成；不能借用同体积/同质心的其他几何。');
assert(isreal(volume)&&isreal(centroid)&&isreal(moments.volume)&&isreal(moments.centroid)&& ...
    all(isfinite(volume)&volume>=0)&&all(isfinite(moments.volume(:)))&& ...
    all(isfinite(moments.centroid(fluid,:)),'all'),'CutQuadratic:Moments','体积或流体质心非有限。');
Cov=moments.covariance;
assert(isreal(Cov)&&isequal(size(Cov),[n 3 3])&&all(isfinite(Cov(fluid,:,:)),'all'), ...
    'CutQuadratic:Moments','缺少真实实数CV二阶矩。');
assert(max(abs(moments.volume(:)-volume)./max(volume,realmin))<1e-9,'CutQuadratic:Moments','CV体积与矩身份不一致。');
assert(max(abs(moments.centroid(fluid,:)-centroid(fluid,:)),[],'all')<1e-10,'CutQuadratic:Moments','CV质心与矩身份不一致。');
for id=find(fluid)'
    cv=reshape(Cov(id,:,:),3,3);assert(norm(cv-cv','fro')<1e-10*max(norm(cv,'fro'),realmin)&& ...
        min(eig((cv+cv')/2))>=-1e-10*max(norm(cv,'fro'),realmin),'CutQuadratic:Moments','非物理Cov，未裁剪本征值。');
end
f=cut.openFaces;positive=f.area>0 & f.cellMinus>0 & f.cellPlus>0;
m=f.cellMinus(positive);p=f.cellPlus(positive);a=f.area(positive);
adj=sparse([m;p],[p;m],[a;a],n,n);samples=boundarySamples(cut,boundary);
sampleByCell=cell(n,1);for k=1:numel(samples.owner),id=samples.owner(k);sampleByCell{id}(end+1)=k;end
local=cell(n,1);rings=zeros(n,1);relativeSV=nan(n,1);rows=zeros(n,1);rankValue=zeros(n,1);
for id=find(fluid)'
    L=volume(id)^(1/3);ci=centroid(id,:);covi=reshape(Cov(id,:,:),3,3);bids=sampleByCell{id}(:);
    visited=false(n,1);visited(id)=true;front=id;neighbors=[];accepted=false;
    for ring=1:options.maximumRing
        candidates=find(any(adj(front,:)>0,1));candidates=candidates(:);new=candidates(~visited(candidates));
        visited(new)=true;neighbors=[neighbors;new];front=new; %#ok<AGROW>
        nb=numel(neighbors);ns=numel(bids);design=zeros(nb+ns,9);weight=zeros(nb+ns,1);
        for k=1:nb
            j=neighbors(k);delta=centroid(j,:)-ci;M=delta'*delta+reshape(Cov(j,:,:),3,3)-covi;
            design(k,:)=basis(delta/L,M/L^2);
            if adj(id,j)>0,area=full(adj(id,j));else,area=min(volume(id),volume(j))^(2/3);end
            % 权重尺度下限只避免权重奇异，不改变真实坐标/体积/墙距或任何输入值。
            weight(k)=sqrt(area)/max(norm(delta),.25*L);
        end
        for k=1:ns
            s=bids(k);delta=samples.point(s,:)-ci;
            if samples.kind(s)==1
                design(nb+k,:)=basis(delta/L,(delta'*delta-covi)/L^2);
            else
                normal=samples.normal(s,:);d=delta/L;
                design(nb+k,:)=[normal,normal.*d,normal(1)*d(2)+normal(2)*d(1), ...
                    normal(1)*d(3)+normal(3)*d(1),normal(2)*d(3)+normal(3)*d(2)];
            end
            weight(nb+k)=sqrt(samples.area(s))/max(norm(delta),.25*L);
        end
        weighted=design.*weight;columnScale=sqrt(sum(weighted.^2,1));columnScale(columnScale==0)=1;
        scaled=weighted./columnScale;[U,S,V]=svd(scaled,0);s=diag(S);
        if numel(s)>=9&&s(1)>0,rankValue(id)=nnz(s>options.minimumRelativeSingularValue*s(1));relativeSV(id)=s(9)/s(1);
        elseif ~isempty(s)&&s(1)>0,rankValue(id)=nnz(s>options.minimumRelativeSingularValue*s(1));relativeSV(id)=0;
        else,rankValue(id)=0;relativeSV(id)=0;end
        if rankValue(id)==9&&relativeSV(id)>=options.minimumRelativeSingularValue
            inverse=diag(1./columnScale)*V*diag(1./s)*U';inverse=inverse.*weight';
            local{id}=struct('neighbors',neighbors,'boundarySampleIDs',bids,'length',L, ...
                'inverse',inverse,'design',design,'weight',weight,'singularValues',s,'columnScale',columnScale);
            rings(id)=ring;rows(id)=size(design,1);accepted=true;break;
        end
        if isempty(front),break;end
    end
    assert(accepted,'CutQuadratic:InsufficientStencil', ...
        'CV%d在%d圈正孔隙邻居后仍不能识别9个二次系数：rank=%d, relativeSV=%.9g；没有降为线性或补零。', ...
        id,ring,rankValue(id),relativeSV(id));
end
stencil=struct('local',{local},'samples',samples,'covariance',Cov,'volume',volume,'centroid',centroid, ...
    'geometryOpenFaces',cut.openFaces,'geometryWalls',cut.walls, ...
    'ringCount',rings,'constraintCount',rows,'scaledRelativeSingularValue',relativeSV,'rank',rankValue, ...
    'options',options,'method','FV平均/Cov修正二次重构；局部长度/列缩放SVD；正孔隙邻接', ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false,'fullNavierStokesStabilityProven',false);
end
function row=basis(d,M)
row=[d,.5*M(1,1),.5*M(2,2),.5*M(3,3),M(1,2),M(1,3),M(2,3)];
end
function s=boundarySamples(c,b)
% 模块用途：明确point温度或零法向梯度，不能把面均值当点值或推断入射剖面。
s=struct('owner',[],'point',zeros(0,3),'normal',zeros(0,3),'area',[],'kind',[], ...
    'sourceID',[],'isCurvedWall',false(0,1),'boundaryName',{{}});
f=c.openFaces;w=c.walls;names={'xmin','xmax';'ymin','ymax';'zmin','zmax'};
for k=find(f.area>0 & (f.cellMinus==0|f.cellPlus==0))'
    if f.cellMinus(k)>0,id=f.cellMinus(k);normal=f.normal(k,:);side=2;else,id=f.cellPlus(k);normal=-f.normal(k,:);side=1;end
    name=names{f.axis(k),side};assert(isfield(b,name),'CutQuadratic:Boundary','缺少明确外边界。');
    add(id,f.centroid(k,:),normal,f.area(k),b.(name),k,false,name);
end
if any(w.ownerCell>0),assert(isfield(b,'solidWall'),'CutQuadratic:Boundary','缺少真实曲面墙热条件。');end
for k=find(w.ownerCell>0)',add(w.ownerCell(k),w.centroid(k,:),w.normal(k,:),w.area(k),b.solidWall,k,true,'solidWall');end
    function add(id,point,normal,area,bc,source,isWall,name)
        if strcmp(bc.type,'temperaturePoint'),kind=1;assert(isfield(bc,'value'),'CutQuadratic:Boundary','须给定点温度。');
        elseif strcmp(bc.type,'zeroNormalGradient'),kind=2;
        else,error('CutQuadratic:Boundary','仅支持明确temperaturePoint或zeroNormalGradient，不自动解释完整NS特征边界。');end
        s.owner(end+1,1)=id;s.point(end+1,:)=point;s.normal(end+1,:)=normal;s.area(end+1,1)=area;
        s.kind(end+1,1)=kind;s.sourceID(end+1,1)=source;s.isCurvedWall(end+1,1)=isWall;s.boundaryName{end+1,1}=name;
    end
end
