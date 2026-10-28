function [moments,details]=buildCutCellVolumeMoments(source,group)
%BUILDCUTCELLVOLUMEMOMENTS 原曲面控制体/保守聚合组的真实二阶体积矩。
% 模块用途：在MATLAB内部由原有向曲面的六项prefix帽积分补全开放面
% 二阶面积矩，再用平面r·n常数恒等式求V、∫r和∫rr'；无需Python运行输入。
% 输入source须为原Cartesian cut-cell，group可选且须保留原格/共享面映射。
% 输出Cov不裁负特征值、不移动几何、不以质心温度冒充二次场体积平均。
% 整面/闭面只继承已验证的source拓扑；真实cut面保留原始积分。
started=tic;if nargin<2 || isempty(group),group=source;aggregated=false;else,aggregated=true;end
[o,h,n,V,F]=validateSource(source);target=validateTarget(source,group,aggregated);
q=(V-o)./h;faces=source.openFaces;nf=numel(faces.area);
area=zeros(nf,1);anchorLocal=zeros(nf,3);first=zeros(nf,3);second=zeros(nf,3,3);
rawUnit=zeros(nf,6);chosenUnit=rawUnit;normalizationCount=0;maximumTopologyCorrection=0;

%% 1. 六个有向投影单项式：[1,u,v,u²,uv,v²]，坐标相对当前tile左下角
for a=1:3
    uv=setdiff(1:3,a);nu=n(uv(1));nv=n(uv(2));na=n(a);
    delta=zeros(na+2,nu,nv,6);direct=zeros(na+1,nu,nv,6);coplanarPositive=zeros(na+1,nu,nv,6);
    for t=1:size(F,1)
        tri=q(F(t,:),:);lo=min(tri,[],1);hi=max(tri,[],1);
        u0=max(1,min(nu,floor(lo(uv(1)))+1));u1=max(1,min(nu,ceil(hi(uv(1)))));
        v0=max(1,min(nv,floor(lo(uv(2)))+1));v1=max(1,min(nv,ceil(hi(uv(2)))));
        for v=v0:v1
            for u=u0:u1
                poly=clip(tri,uv(1),u-1,true);poly=clip(poly,uv(1),u,false);
                poly=clip(poly,uv(2),v-1,true);poly=clip(poly,uv(2),v,false);
                if size(poly,1)<3,continue;end
                ref=zeros(1,3);ref(uv)=[u-1 v-1];full=projectedSix(poly-ref,a,uv);
                lower=min(poly(:,a));upper=max(poly(:,a));start=max(0,ceil(upper));
                if start<=na,delta(start+1,u,v,:)=delta(start+1,u,v,:)-reshape(full,1,1,1,6);end
                for s=max(0,floor(lower)+1):min(na,ceil(upper)-1)
                    p=clip(poly,a,s,false);m=projectedSix(p-ref,a,uv);
                    direct(s+1,u,v,:)=direct(s+1,u,v,:)-reshape(m,1,1,1,6);
                end
                nearest=round(tri(1,a));constant=all(V(F(t,:),a)==V(F(t,1),a));
                if constant && abs(tri(1,a)-nearest)<=64*eps(max([1 n])) && nearest>=0 && nearest<=na && full(1)>0
                    coplanarPositive(nearest+1,u,v,:)=coplanarPositive(nearest+1,u,v,:)+reshape(full,1,1,1,6);
                end
            end
        end
    end
    solid=cumsum(delta(1:na+1,:,:,:),1)+direct+coplanarPositive;
    ids=find(faces.axis==a);lookup=zeros(na+1,nu,nv);
    assert(numel(ids)==(na+1)*nu*nv,'CutMom:FaceTable','原共享面表数量错误。');
    for id=ids'
        gi=double(faces.gridIndex(id,:));sub=[gi(a)+1 gi(uv)];
        assert(all(sub>=1 & sub<=[na+1 nu nv]) && lookup(sub(1),sub(2),sub(3))==0,'CutMom:FaceTable','共享面索引非法/重复。');
        lookup(sub(1),sub(2),sub(3))=id;
    end
    scale=h(uv(1))*h(uv(2));fullSix=[1 .5 .5 1/3 .25 1/3];
    for v=1:nv
        for u=1:nu
            for s=0:na
                id=lookup(s+1,u,v);sm=reshape(solid(s+1,u,v,:),1,6);raw=fullSix-sm;chosen=raw;
                state=faces.topologicalState(id);
                if state==0,chosen=zeros(1,6);elseif state==2,chosen=fullSix;
                else,assert(state==1 && raw(1)>=0 && raw(1)<=1,'CutMom:Aperture','真实cut面面积积分非法，未裁剪。');end
                change=max(abs(chosen-raw));
                assert(change<3e-11,'CutMom:TopologyMismatch','完整/闭面高阶矩与已验证拓扑不一致。');
                normalizationCount=normalizationCount+(change~=0);maximumTopologyCorrection=max(maximumTopologyCorrection,change);
                rawUnit(id,:)=raw;chosenUnit(id,:)=chosen;area(id)=chosen(1)*scale;
                r=zeros(1,3);r(a)=s*h(a);r(uv)=([u v]-.5).*h(uv);anchorLocal(id,:)=r;
                first(id,uv)=scale*h(uv).*(chosen(2:3)-.5*chosen(1));
                S=zeros(3);S(uv(1),uv(1))=scale*h(uv(1))^2*(chosen(4)-chosen(2)+.25*chosen(1));
                S(uv(2),uv(2))=scale*h(uv(2))^2*(chosen(6)-chosen(3)+.25*chosen(1));
                S(uv(1),uv(2))=scale*prod(h(uv))*(chosen(5)-.5*(chosen(2)+chosen(3))+.25*chosen(1));
                S(uv(2),uv(1))=S(uv(1),uv(2));second(id,:,:)=reshape(S,1,3,3);
            end
        end
    end
end
areaDifference=max(abs(area-faces.area));
assert(areaDifference<3e-11*max(h)^2,'CutMom:ApertureMismatch','重建开放面积与source不一致。');
active=area>0;independentCentroid=nan(nf,3);independentCentroid(active,:)=o+anchorLocal(active,:)+first(active,:)./area(active);
centroidDifference=max(abs(independentCentroid(active,:)-faces.centroid(active,:)),[],'all');
assert(centroidDifference<2e-8*max(h),'CutMom:ApertureMismatch','开放面一阶矩/质心与source不一致。');

%% 2. 原墙片在自身首顶点参考下精确积分；单位法向继承原有序三角面
w=source.walls;nw=numel(w.area);wallAnchor=zeros(nw,3);wallFirst=zeros(nw,3);wallSecond=zeros(nw,3,3);wallArea=zeros(nw,1);
maximumWallRelativeAreaError=0;negativeProjectedFanCount=0;
for k=1:nw
    poly=w.polygonVerticesGridLocal{k};anchor=poly(1,:);wallAnchor(k,:)=anchor;S=zeros(3);M=zeros(1,3);A=0;
    for j=2:size(poly,1)-1
        local=poly([1 j j+1],:)-anchor;
        projected=dot(cross(local(2,:),local(3,:))/2,-w.normal(k,:));
        % 原凸片按原平面正投影定义面积，沿用source对浮点dust扇片的定义。
        % 保留每个正片，记录负投影诊断；不使用dust叉积范数来建立法向。
        negativeProjectedFanCount=negativeProjectedFanCount+(projected<0);ar=max(0,projected);total=sum(local,1);
        A=A+ar;M=M+ar*total/3;S=S+ar*(total'*total+local'*local)/12;
    end
    relative=abs(A-w.area(k))/w.area(k);maximumWallRelativeAreaError=max(maximumWallRelativeAreaError,relative);
    assert(relative<128*eps,'CutMom:WallArea','保留原片面积与可靠平面积分不一致。');
    wallArea(k)=A;wallFirst(k,:)=M;wallSecond(k,:,:)=reshape(S,1,3,3);
end

%% 3. 在每个真实体积质心的局部参考下组装V、∫r、∫rr'，不构造全局大坐标矩
volume=target.cells.volume(:);center=target.cells.centroid;N=numel(volume);positive=volume>0;
reconstructedVolume=zeros(N,1);M1=zeros(N,3);M2=zeros(N,3,3);boundaryArea=zeros(N,1);
tf=target.openFaces;if aggregated,original=double(tf.originalCartesianFaceID(:));else,original=(1:nf)';end
for k=1:numel(tf.area)
    idFace=original(k);if area(idFace)==0,continue;end
    for side=1:2
        if side==1,id=tf.cellMinus(k);sgn=1;else,id=tf.cellPlus(k);sgn=-1;end
        if id==0,continue;end
        offset=anchorLocal(idFace,:)-(center(id,:)-o);normal=sgn*faces.normal(idFace,:);
        addBoundary(id,area(idFace),first(idFace,:),reshape(second(idFace,:,:),3,3),offset,normal);
    end
end
for k=1:nw
    id=target.walls.ownerCell(k);if id==0,continue;end
    offset=wallAnchor(k,:)-(center(id,:)-o);
    addBoundary(id,wallArea(k),wallFirst(k,:),reshape(wallSecond(k,:,:),3,3),offset,w.normal(k,:));
end
assert(all(reconstructedVolume(positive)>0),'CutMom:Volume','真实流体体积重建非正，未裁剪。');
relativeVolume=nan(N,1);relativeVolume(positive)=abs(reconstructedVolume(positive)-volume(positive))./volume(positive);
assert(max(relativeVolume,[],'omitnan')<3e-8,'CutMom:VolumeMismatch','局部边界体积与source/group不一致。');
meanOffset=nan(N,3);meanOffset(positive,:)=M1(positive,:)./reconstructedVolume(positive);
assert(max(abs(meanOffset(positive,:)),[],'all')<2e-8*max(h),'CutMom:CentroidMismatch','局部一阶矩与source/group质心不一致。');
covariance=nan(N,3,3);centralSecond=nan(N,3,3);minEigen=inf;maxSymmetry=0;
for id=find(positive)'
    C=reshape(M2(id,:,:),3,3)/reconstructedVolume(id)-meanOffset(id,:)'*meanOffset(id,:);
    symmetry=max(abs(C-C'),[],'all');maxSymmetry=max(maxSymmetry,symmetry);
    assert(symmetry<512*eps*max(h)^2,'CutMom:Symmetry','真实Cov矩阵不对称。');
    eigenvalues=eig((C+C')/2);minEigen=min(minEigen,min(eigenvalues));
    assert(all(eigenvalues>0),'CutMom:Covariance','真实Cov非正定；未裁特征值或补造二阶矩。');
    covariance(id,:,:)=reshape(C,1,3,3);centralSecond(id,:,:)=reconstructedVolume(id)*reshape(C,1,3,3);
end
moments=struct('volume',reconstructedVolume,'sourceVolume',volume,'referenceCenter',center, ...
    'centroid',center+meanOffset,'centralSecondMoment',centralSecond,'covariance',covariance, ...
    'secondMomentAboutReferenceCenter',M2,'firstMomentAboutReferenceCenter',M1, ...
    'physicalSolution',false,'eligibleForSurfaceCoverage',false, ...
    'method','MATLAB oriented prefix six-monomial apertures + source-plane polygon moments + local homogeneous divergence', ...
    'units',struct('volume','m^3','firstVolumeMoment','m^4','secondVolumeMoment','m^5','covariance','m^2'), ...
    'stateInterpretation','真实FV几何矩；可用于二次保守量平均，不能把质心primitive样本冒充平均');
% 模块用途：矩与实际CV边界绑定；只校验体积/质心不足以辨别不同形状的Cov。
% 保存完整共享面/墙几何以便严格比较，不以近似数值摘要冒充几何指纹。
moments.geometryIdentity=struct('volume',target.cells.volume(:),'centroid',target.cells.centroid, ...
    'openFaces',target.openFaces,'walls',target.walls);
% 面/墙Cov关于各自面积质心；给定point壁温不等于非仿射温度的面积平均。
faceCovariance=nan(nf,3,3);wallCovariance=nan(nw,3,3);wallCentroid=zeros(nw,3);
for k=find(area>0)'
    mu=first(k,:)/area(k);faceCovariance(k,:,:)=reshape(reshape(second(k,:,:),3,3)/area(k)-mu'*mu,1,3,3);
end
for k=1:nw
    mu=wallFirst(k,:)/wallArea(k);wallCovariance(k,:,:)=reshape(reshape(wallSecond(k,:,:),3,3)/wallArea(k)-mu'*mu,1,3,3);
    wallCentroid(k,:)=o+wallAnchor(k,:)+mu;
end
details=struct('openFaces',struct('area',area,'referencePointGridLocal',anchorLocal,'referencePoint',o+anchorLocal, ...
    'firstMoment',first,'secondMoment',second,'covariance',faceCovariance,'centroidFromMoments',independentCentroid, ...
    'rawUnitMoments',rawUnit,'chosenUnitMoments',chosenUnit), ...
    'walls',struct('area',wallArea,'referencePointGridLocal',wallAnchor,'referencePoint',o+wallAnchor, ...
    'firstMoment',wallFirst,'secondMoment',wallSecond,'covariance',wallCovariance,'centroidFromMoments',wallCentroid), ...
    'targetFaceSourceID',original,'targetWallOwner',target.walls.ownerCell, ...
    'maximumApertureAreaDifference',areaDifference,'maximumApertureCentroidDifference',centroidDifference, ...
    'maximumWallRelativeAreaError',maximumWallRelativeAreaError,'negativeProjectedFanCount',negativeProjectedFanCount, ...
    'topologicalFaceNormalizationCount',normalizationCount,'maximumTopologicalUnitMomentCorrection',maximumTopologyCorrection, ...
    'maximumRelativeVolumeDifference',max(relativeVolume,[],'omitnan'), ...
    'maximumCentroidOffset',max(abs(meanOffset(positive,:)),[],'all'), ...
    'minimumCovarianceEigenvalue',minEigen,'maximumCovarianceSymmetryError',maxSymmetry, ...
    'boundaryArea',boundaryArea,'aggregated',aggregated,'elapsedSeconds',toc(started), ...
    'covarianceClippedOrSymmetrized',false,'physicalSolution',false,'eligibleForSurfaceCoverage',false);

    function addBoundary(id,A,F1,F2,offset,normal)
        d=dot(offset,normal);firstLocal=F1+A*offset;
        secondLocal=F2+offset'*F1+F1'*offset+A*(offset'*offset);
        reconstructedVolume(id)=reconstructedVolume(id)+d*A/3;
        M1(id,:)=M1(id,:)+d*firstLocal/4;M2(id,:,:)=M2(id,:,:)+reshape(d*secondLocal/5,1,3,3);
        boundaryArea(id)=boundaryArea(id)+A;
    end
end

function [o,h,n,V,F]=validateSource(s)
% 模块用途：显式验证原格/共享面/有向曲面合同；离域或损坏输入不能伪造矩。
required={'grid','cells','openFaces','walls','inputGeometry'};
assert(all(isfield(s,required)),'CutMom:Schema','缺少原始cut-cell结构。');
o=double(s.grid.origin(:)');h=double(s.grid.spacing(:)');n=double(s.grid.cellCounts(:)');
if isscalar(h),h=repmat(h,1,3);end
assert(numel(o)==3 && numel(h)==3 && numel(n)==3 && all(isfinite([o h n])) && all(h>0) && all(n>=1 & n==fix(n)), ...
    'CutMom:Grid','原Cartesian网格参数非法。');
V=double(s.inputGeometry.V);F=double(s.inputGeometry.F);
assert(size(V,2)==3 && all(isfinite(V),'all') && size(F,2)==3 && all(F>=1 & F<=size(V,1) & F==fix(F),'all'), ...
    'CutMom:Geometry','原三角曲面非法。');
local=(V-o)./h;tol=128*eps(max([1 n]));
assert(all(local>=-tol & local<=n+tol,'all'),'CutMom:OutsideDomain','实体顶点离开原计算域；不能截掉实体伪造体积矩。');
[~,~,map]=unique(V,'rows');f=reshape(map(F),size(F));edges=[f(:,[1 2]);f(:,[2 3]);f(:,[3 1])];
[~,~,j]=unique(sort(edges,2),'rows');count=accumarray(j,1);balance=accumarray(j,2*(edges(:,1)<edges(:,2))-1);
assert(all(count==2 & balance==0),'CutMom:Orientation','曲面须闭合且共享边反向。');
r=mean(V,1);a=V(F(:,1),:)-r;b=V(F(:,2),:)-r;c=V(F(:,3),:)-r;
assert(sum(sum(a.*cross(b,c,2),2))/6>0,'CutMom:Orientation','原实体有向体积必须为正。');
vol=s.cells.volume(:);assert(numel(vol)==prod(n) && all(isfinite(vol)&vol>=0),'CutMom:Cells','source须为全原格，体积不能非法。');
assert(any(vol>0),'CutMom:NoFluid','没有流体控制体。');
assert(isequal(size(s.cells.centroid),[numel(vol) 3]) && all(isfinite(s.cells.centroid(vol>0,:)),'all'),'CutMom:Cells','流体质心非法。');
assert(all(isfield(s.openFaces,{'axis','gridIndex','normal','area','centroid','cellMinus','cellPlus','topologicalState'})), ...
    'CutMom:FaceTable','必须保存原全局共享面及拓扑状态。');
assert(all(isfield(s.walls,{'normal','area','areaVector','ownerCell','originalFaceID','polygonVerticesGridLocal'})), ...
    'CutMom:WallSchema','墙须保留可靠原法向、原面身份与grid-local多边形。');
crossOriginal=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);length=sqrt(sum(crossOriginal.^2,2));
assert(all(length>0),'CutMom:Geometry','源三角形退化。');normal=-crossOriginal./length;
w=s.walls;assert(all(isfinite(w.area)&w.area>0),'CutMom:WallArea','原正墙片不能删除或补零。');
assert(max(abs(w.normal-normal(w.originalFaceID,:)),[],'all')<64*eps,'CutMom:WallNormal','墙法向须继承可靠原三角面。');
assert(max(abs(w.areaVector-w.area.*w.normal),[],'all')<128*eps*max(w.area),'CutMom:WallNormal','面积向量与可靠法向不一致。');
end

function target=validateTarget(source,target,aggregated)
if ~aggregated,return;end
assert(isequaln(source.inputGeometry,target.inputGeometry) && isequaln(source.grid,target.grid),'CutMom:GroupGeometry','group与source必须对应同一几何/原网格。');
assert(isfield(target,'originalCellToGroup') && isfield(target.openFaces,'originalCartesianFaceID'),'CutMom:GroupMap','缺少原格/原共享面到组的映射。');
v=source.cells.volume(:);g=target.cells.volume(:);map=double(target.originalCellToGroup(:));positive=v>0;
assert(numel(map)==numel(v) && all(map(~positive)==0) && all(map(positive)>=1 & map(positive)<=numel(g) & map(positive)==fix(map(positive))), ...
    'CutMom:GroupMap','原流体格到组的映射非法。');
summed=accumarray(map(positive),v(positive),[numel(g) 1]);
assert(all(abs(summed-g)<=512*eps*max(summed,g)),'CutMom:GroupVolume','聚合组体积不是原格守恒和。');
reference=source.grid.origin;first=source.cells.volume(positive).*(source.cells.centroid(positive,:)-reference);
center=zeros(numel(g),3);for d=1:3,center(:,d)=reference(d)+accumarray(map(positive),first(:,d),[numel(g) 1])./summed;end
assert(max(abs(center-target.cells.centroid),[],'all')<512*eps(max(1,max(abs(center),[],'all'))),'CutMom:GroupCentroid','组质心不是守恒体积矩。');
ids=double(target.openFaces.originalCartesianFaceID(:));f=source.openFaces;tf=target.openFaces;
assert(all(ids>=1 & ids<=numel(f.area)) && numel(unique(ids))==numel(ids),'CutMom:GroupMap','原共享面映射重复/非法。');
assert(isequaln(tf.area,f.area(ids)) && isequaln(tf.normal,f.normal(ids,:)) && isequaln(tf.centroid,f.centroid(ids,:)), ...
    'CutMom:GroupMap','组共享面几何必须保持原样。');
fm=zeros(size(f.cellMinus));fp=fm;use=f.cellMinus>0;fm(use)=map(f.cellMinus(use));use=f.cellPlus>0;fp(use)=map(f.cellPlus(use));
keep=f.area>0 & fm~=fp;assert(isequal(ids,find(keep)) && isequaln(double(tf.cellMinus),double(fm(ids))) && isequaln(double(tf.cellPlus),double(fp(ids))), ...
    'CutMom:GroupMap','组内面须消失、外露原共享面须唯一且左右归属正确。');
w=source.walls;tw=target.walls;owner=zeros(size(w.ownerCell));use=w.ownerCell>0;owner(use)=map(w.ownerCell(use));
assert(isequaln(tw.ownerCell,owner) && isequaln(tw.area,w.area) && isequaln(tw.normal,w.normal) && isequaln(tw.originalFaceID,w.originalFaceID) && ...
    isequaln(tw.centroid,w.centroid) && isequaln(tw.polygonVertices,w.polygonVertices) && ...
    isequaln(tw.polygonVerticesGridLocal,w.polygonVerticesGridLocal), ...
    'CutMom:GroupMap','全部原墙片及真正fluid-owner必须保持。');
end

function m=projectedSix(poly,a,uv)
m=zeros(1,6);if size(poly,1)<3,return;end
anchor=poly(1,:);
for k=2:size(poly,1)-1
    p=poly([1 k k+1],uv);av=cross(poly(k,:)-anchor,poly(k+1,:)-anchor)/2;A=av(a);total=sum(p,1);
    m=m+A*[1 total/3 (total(1)^2+sum(p(:,1).^2))/12 ...
        (prod(total)+sum(p(:,1).*p(:,2)))/12 (total(2)^2+sum(p(:,2).^2))/12];
end
end

function out=clip(poly,a,s,lower)
if isempty(poly),out=zeros(0,3);return;end
out=zeros(size(poly,1)+2,3);used=0;previous=poly(end,:);d0=previous(a)-s;
if lower,in0=d0>=0;else,in0=d0<=0;end
for k=1:size(poly,1)
    current=poly(k,:);d1=current(a)-s;if lower,in1=d1>=0;else,in1=d1<=0;end
    if in0~=in1
        p=previous+d0/(d0-d1)*(current-previous);p(a)=s;
        if used==0 || any(out(used,:)~=p),used=used+1;out(used,:)=p;end
    end
    if in1 && (used==0 || any(out(used,:)~=current)),used=used+1;out(used,:)=current;end
    previous=current;d0=d1;in0=in1;
end
out=out(1:used,:);if size(out,1)>1 && all(out(1,:)==out(end,:)),out(end,:)=[];end
end
