function cut = buildOriginalSurfaceCutCells(V,F,grid)
%BUILDORIGINALSURFACECUTCELLS 原始闭合定向曲面与笛卡尔网格的流体控制体。
% 中文模块说明：利用有向曲面的平面帽积分构造共享开放面，随后在每个
% 单元的局部坐标内积分真实边界，得到流体体积和体积质心。非凸曲面
% 保留原始三角面；不体素化、不用等效圆、不把负体积裁成零。
% 共面墙的几何半开归属须改成实体外法向一侧的真实流体所有者。
% 几何整格/整面只按无内部曲面及独立绕数分类规范化；保留原始积分值。

% 同目录完整裁片模块由setupCutCellFlow加入路径；不依赖outputs草稿目录。
parts=partitionOrientedTriangles(V,F,grid);
grid=parts.grid; o=grid.origin; h=grid.spacing; n=grid.cellCounts;
q=(double(V)-o)./h; F=double(F); N=prod(n); cellVol=prod(h);
validateClosedOrientation(q,F);
[ii,jj,kk]=ndgrid(1:n(1),1:n(2),1:n(3)); idx=[ii(:) jj(:) kk(:)];
centers=o+(idx-.5).*h;

%% 原曲面墙面：跨格片归所在格，共面片归真实流体侧
walls=parts;
sourceCross=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
sourceNormal=sourceCross./sqrt(sum(sourceCross.^2,2));
walls.normal=-sourceNormal(parts.originalFaceID,:);
walls.rawPolygonAreaVector=-parts.areaVector;
walls.areaVector=parts.area.*walls.normal;
walls.areaVectorCorrection=walls.areaVector-walls.rawPolygonAreaVector;
correctionNorm=sqrt(sum(walls.areaVectorCorrection.^2,2));
walls.normalCorrection=struct('source','negative normalized cross product of original ordered source triangle edges', ...
    'changedPatchCount',sum(any(walls.areaVectorCorrection~=0,2)), ...
    'maxAbsoluteAreaVectorChange',max(correctionNorm), ...
    'maxRelativeAreaVectorChange',max(correctionNorm./walls.area), ...
    'maxUnitNormalError',max(abs(sqrt(sum(walls.normal.^2,2))-1)));
assert(walls.normalCorrection.maxUnitNormalError<=16*eps);
walls.polygonVerticesGridLocal=cell(size(parts.polygonVertices));
walls.polygonVerticesCellLocal=cell(size(parts.polygonVertices));
walls.geometricCellID=parts.cellID; walls.geometricCellSubscripts=parts.cellSubscripts;
walls.ownerCell=parts.cellID; ownerSubscripts=parts.cellSubscripts; interiorWall=false(N,1);
for p=1:numel(parts.area)
    sub=parts.cellSubscripts(p,:);
    original=F(parts.originalFaceID(p),:); moved=false;
    % 网格局部顶点重放原六平面裁切，避免极小正片在加上大原点后丢失位数。
    poly=recoverPatch(q(original,:),V(original,:),sub,n);
    walls.polygonVerticesGridLocal{p}=poly.*h;
    for a=1:3
        if all(V(original,a)==V(original(1),a))
            s=round(q(original(1),a));
            if abs(q(original(1),a)-s)<=64*eps(max([1 n]))
                outward=parts.areaVector(p,a);
                if outward>0, sub(a)=s+1; elseif outward<0, sub(a)=s; end
                moved=true;
            end
        end
    end
    if any(sub<1 | sub>n), walls.ownerCell(p)=0; ownerSubscripts(p,:)=0;
    else
        id=sub(1)+(sub(2)-1)*n(1)+(sub(3)-1)*n(1)*n(2);
        walls.ownerCell(p)=id;
        ownerSubscripts(p,:)=sub;
        walls.polygonVerticesCellLocal{p}=(poly-(sub-.5)).*h;
        if ~moved && any(all(poly>sub-1 & poly<sub,2))
            interiorWall(id)=true;
        elseif ~moved
            % 即使没有内部顶点，面片内部仍可跨越单元内部。
            c=(parts.centroid(p,:)-o)./h;
            if all(c>sub-1 & c<sub), interiorWall(id)=true; end
        end
    end
end
walls.fluidOwnerCellID=walls.ownerCell;
walls.cellID=walls.ownerCell; walls.cellSubscripts=ownerSubscripts;
walls.fluidOwnerDetermined=true;
walls.orientationConvention='fluid outward normal = negative original body outward normal';
walls.ownershipPolicy='original body outward side fluid cell; geometricCellID separately retains half-open clipping ownership; zero means exterior domain';

%% 三轴共享面：曲面前缀帽积分，矩在每个面格的切向局部坐标中计算
face=struct('axis',[],'gridIndex',[],'normal',[],'area',[],'centroid',[], ...
    'cellMinus',[],'cellPlus',[],'rawArea',[],'rawCentroid',[],'topologicalState',[]);
faceNormalizationCount=0; maxFaceAreaChange=0; maxFaceMomentChange=0;
for a=1:3
    uv=setdiff(1:3,a); nu=n(uv(1)); nv=n(uv(2)); na=n(a);
    delta=zeros(na+2,nu,nv,3); direct=zeros(na+1,nu,nv,3);
    coplanarPositive=zeros(na+1,nu,nv,3);
    boundary=false(na+1,nu,nv); coplanar=false(na+1,nu,nv);
    for t=1:size(F,1)
        tri=q(F(t,:),:); lo=min(tri,[],1); hi=max(tri,[],1);
        u0=max(1,min(nu,floor(lo(uv(1)))+1)); u1=max(1,min(nu,ceil(hi(uv(1)))));
        v0=max(1,min(nv,floor(lo(uv(2)))+1)); v1=max(1,min(nv,ceil(hi(uv(2)))));
        if u1<u0 || v1<v0, continue; end
        for v=v0:v1
            for u=u0:u1
                poly=tri;
                poly=clip(poly,uv(1),u-1,true); poly=clip(poly,uv(1),u,false);
                poly=clip(poly,uv(2),v-1,true); poly=clip(poly,uv(2),v,false);
                if size(poly,1)<3, continue; end
                ref=zeros(1,3); ref(uv)=[u-1 v-1];
                full=projectedMoments(poly-ref,a,uv);
                lower=min(poly(:,a)); upper=max(poly(:,a));
                start=max(0,ceil(upper));
                if start<=na
                    delta(start+1,u,v,:)=delta(start+1,u,v,:)-reshape(full,1,1,1,3);
                end
                for s=max(0,floor(lower)+1):min(na,ceil(upper)-1)
                    p=clip(poly,a,s,false);
                    m=projectedMoments(p-ref,a,uv);
                    direct(s+1,u,v,:)=direct(s+1,u,v,:)-reshape(m,1,1,1,3);
                end
                constant=all(V(F(t,:),a)==V(F(t,1),a)); nearest=round(tri(1,a));
                if constant && abs(tri(1,a)-nearest)<=64*eps(max([1 n])) && nearest>=0 && nearest<=na
                    if full(1)~=0, coplanar(nearest+1,u,v)=true; end
                    if full(1)>0
                        coplanarPositive(nearest+1,u,v,:)=coplanarPositive(nearest+1,u,v,:)+reshape(full,1,1,1,3);
                    end
                elseif lower<upper
                    for s=max(0,ceil(lower)):min(na,floor(upper))
                        points=planeLine(poly,a,s);
                        if size(points,1)>=2
                            mid=mean(points,1);
                            if any(points(1,:)~=points(end,:)) && all(mid(uv)>[u-1 v-1] & mid(uv)<[u v])
                                boundary(s+1,u,v)=true;
                            end
                        end
                    end
                end
            end
        end
    end
    solid=cumsum(delta(1:na+1,:,:,:),1)+direct+coplanarPositive;
    nf=(na+1)*nu*nv; ax=repmat(a,nf,1); gi=zeros(nf,3); area=zeros(nf,1); rawArea=area;
    centroid=nan(nf,3); rawCentroid=centroid; minus=zeros(nf,1); plus=minus; normal=zeros(nf,3); normal(:,a)=1;
    state=zeros(nf,1); used=0; faceScale=h(uv(1))*h(uv(2));
    for v=1:nv
        for u=1:nu
            for s=0:na
                used=used+1; g=zeros(1,3); g(a)=s; g(uv)=[u v]; gi(used,:)=g;
                gminus=g; gplus=g; gminus(a)=s; gplus(a)=s+1;
                if s>=1, minus(used)=linearID(gminus,n); end
                if s<na, plus(used)=linearID(gplus,n); end
                sm=reshape(solid(s+1,u,v,:),1,3); raw=([1 .5 .5]-sm); chosen=raw;
                if ~boundary(s+1,u,v)
                    p=zeros(1,3); p(a)=s; p(uv)=[u-.5 v-.5];
                    inside=coplanar(s+1,u,v) || pointInside(p,q,F);
                    if inside, chosen=[0 0 0]; state(used)=0;
                    else, chosen=[1 .5 .5]; state(used)=2; end
                    change=abs(chosen-raw);
                    if any(change>2e-12)
                        error('CutCell:TopologicalFaceMismatch','无内部截线面的积分与独立拓扑不一致 axis=%d, grid=%s.',a,mat2str(g));
                    end
                    if any(change~=0), faceNormalizationCount=faceNormalizationCount+1; end
                    maxFaceAreaChange=max(maxFaceAreaChange,change(1)*faceScale);
                    maxFaceMomentChange=max(maxFaceMomentChange,max(change(2:3))*faceScale*max(h(uv)));
                else
                    state(used)=1;
                    if chosen(1)<0 || chosen(1)>1
                        error('CutCell:InvalidAperture','真实截切面的原始开放面积超出[0,1]: %.17g axis=%d grid=%s.',chosen(1),a,mat2str(g));
                    end
                end
                rawArea(used)=raw(1)*faceScale; area(used)=chosen(1)*faceScale;
                if raw(1)~=0
                    cq=zeros(1,3); cq(a)=s; cq(uv)=[u-1 v-1]+raw(2:3)/raw(1); rawCentroid(used,:)=o+cq.*h;
                end
                if chosen(1)>0
                    cq=zeros(1,3); cq(a)=s; cq(uv)=[u-1 v-1]+chosen(2:3)/chosen(1); centroid(used,:)=o+cq.*h;
                end
            end
        end
    end
    names=fieldnames(face); values={ax,gi,normal,area,centroid,minus,plus,rawArea,rawCentroid,state};
    for k=1:numel(names), face.(names{k})=[face.(names{k}); values{k}]; end
    fprintf('  shared Cartesian axis %d complete: %d faces\n',a,nf);
end
face.lowerCellID=face.cellMinus; face.upperCellID=face.cellPlus;

%% 每个单元局部边界矩：体积张量与精确二阶墙面积矩
av=zeros(N,3); tensor=zeros(N,3,3); second=zeros(N,3); boundaryArea=zeros(N,1);
for f=1:numel(face.area)
    if face.area(f)==0, continue; end
    for side=1:2
        if side==1, id=face.cellMinus(f); sgn=1; else, id=face.cellPlus(f); sgn=-1; end
        if id==0, continue; end
        an=sgn*face.area(f)*face.normal(f,:); r=face.centroid(f,:)-centers(id,:);
        av(id,:)=av(id,:)+an; tensor(id,:,:)=tensor(id,:,:)+reshape(r'*an,1,3,3);
        boundaryArea(id)=boundaryArea(id)+face.area(f);
        second(id,:)=second(id,:)+.5*(r.^2).*an;
    end
end
for p=1:numel(walls.area)
    id=walls.ownerCell(p); if id==0, continue; end
    an=walls.areaVector(p,:); r=walls.centroid(p,:)-centers(id,:);
    av(id,:)=av(id,:)+an; tensor(id,:,:)=tensor(id,:,:)+reshape(r'*an,1,3,3);
    boundaryArea(id)=boundaryArea(id)+walls.area(p);
    local=walls.polygonVerticesGridLocal{p}; poly=walls.polygonVerticesCellLocal{p};
    nwall=walls.normal(p,:); anchor=poly(1,:); anchorGrid=local(1,:); quadratic=zeros(1,3); integratedArea=0;
    for k=2:size(poly,1)-1
        b=poly(k,:); c=poly(k+1,:);
        % 沿原始可靠法向投影，与原裁切模块的面积定义完全相同；不取dust叉积范数。
        projected=dot(cross(local(k,:)-anchorGrid,local(k+1,:)-anchorGrid)/2,-nwall);
        ar=max(0,projected); integratedArea=integratedArea+ar;
        quadratic=quadratic+ar*(anchor.^2+b.^2+c.^2+anchor.*b+b.*c+c.*anchor)/6;
    end
    if abs(integratedArea-walls.area(p))>64*eps*walls.area(p)
        error('CutCell:WallMomentAreaMismatch','原面法向二阶面积积分与保留patch面积不一致 patch=%d.',p);
    end
    second(id,:)=second(id,:)+.5*quadratic.*nwall;
end
rawVolume=(tensor(:,1,1)+tensor(:,2,2)+tensor(:,3,3))/3;
volume=rawVolume; rawCellCentroid=centers+second./rawVolume; cellCentroid=rawCellCentroid;
cellState=ones(N,1); cellNormalizationCount=0; maxVolumeChange=0;
for id=1:N
    if ~interiorWall(id)
        inside=pointInside((centers(id,:)-o)./h,q,F);
        if inside, exact=0; cellState(id)=0; cellCentroid(id,:)=NaN;
        else, exact=cellVol; cellState(id)=2; cellCentroid(id,:)=centers(id,:); end
        change=abs(rawVolume(id)-exact);
        if change>2e-12*cellVol
            error('CutCell:TopologicalCellMismatch','无内部曲面单元积分与独立绕数不一致 cell=%d raw=%.17g exact=%.17g.',id,rawVolume(id),exact);
        end
        volume(id)=exact;
        if change~=0, cellNormalizationCount=cellNormalizationCount+1; end
        maxVolumeChange=max(maxVolumeChange,change);
    elseif ~(rawVolume(id)>0 && rawVolume(id)<cellVol)
        error('CutCell:InvalidFluidVolume','真实切割单元体积必须严格处于(0,Vcell)，未裁剪 cell=%d raw=%.17g.',id,rawVolume(id));
    end
end
tensorResidual=tensor;
for a=1:3, tensorResidual(:,a,a)=tensorResidual(:,a,a)-volume; end
areaScale=sum(abs(av),2); %#ok<NASGU>
areaError=max(abs(av),[],'all')/max(h)^2;
tensorError=max(abs(tensorResidual),[],'all')/cellVol;
areaRelative=nan(N,1); momentRelative=nan(N,1); positive=volume>0;
areaRelative(positive)=max(abs(av(positive,:)),[],2)./boundaryArea(positive);
momentRelative(positive)=max(reshape(abs(tensorResidual(positive,:,:)),sum(positive),9),[],2)./volume(positive);
if areaError>2e-10 || tensorError>2e-10
    error('CutCell:LocalClosureFailure','逐单元闭合失败 area=%.17g tensor=%.17g.',areaError,tensorError);
end
distance=nan(numel(walls.area),1);
for p=1:numel(walls.area)
    id=walls.ownerCell(p);
    if id>0, distance(p)=dot(cellCentroid(id,:)-walls.centroid(p,:),-walls.areaVector(p,:)/walls.area(p)); end
end
cut=struct('grid',grid,'cells',struct('volume',volume,'centroid',cellCentroid,'center',centers, ...
    'index',idx,'state',cellState,'rawVolume',rawVolume,'rawCentroid',rawCellCentroid), ...
    'openFaces',face,'walls',walls,'closure',struct('areaVectorResidual',av,'firstMomentResidual',tensorResidual, ...
    'rawFirstMomentTensor',tensor,'boundaryArea',boundaryArea,'relativeAreaErrorPerCell',areaRelative, ...
    'relativeTensorErrorPerFluidVolume',momentRelative,'maxRelativeFluidTensorError',max(momentRelative,[],'omitnan'), ...
    'normalizedAreaError',areaError,'normalizedTensorError',tensorError));
cut.closureAreaVector=av; cut.closureTensor=tensorResidual;
cut.normalization=struct('method','independent no-interior-boundary topological full/solid classification; never volume clipping', ...
    'faceCount',faceNormalizationCount,'maxFaceAreaChange',maxFaceAreaChange,'maxFaceTangentialMomentChange',maxFaceMomentChange, ...
    'cellCount',cellNormalizationCount,'maxVolumeChange',maxVolumeChange,'comparisonBoundNormalized',2e-12);
cut.wallOwnerCentroidSignedDistance=distance;
cut.nonpositiveWallOwnerDistancePatchIDs=find(distance<=0);
cut.algorithm='original oriented surface cap integrals and local exact polyhedral divergence moments';
cut.inputGeometry=struct('V',V,'F',F); cut.physicalSolution=false; cut.eligibleForSurfaceCoverage=false;
cut.fluidVolumesImplemented=true; cut.openFaceAperturesImplemented=true; cut.perCellClosureImplemented=true;
end

function id=linearID(g,n)
id=g(1)+(g(2)-1)*n(1)+(g(3)-1)*n(1)*n(2);
end

function poly=recoverPatch(tri,original,sub,n)
% 中文模块说明：重放原模块的共面锁定与半空间次序，保留网格局部有效位。
poly=tri;
for a=1:3
    nearest=round(tri(1,a));
    locked=all(original(:,a)==original(1,a)) && abs(tri(1,a)-nearest)<=64*eps(max([1 n]));
    if locked, continue; end
    poly=clip(poly,a,sub(a)-1,true); poly=clip(poly,a,sub(a),false);
end
end

function m=projectedMoments(poly,a,uv)
% 中文模块说明：凸裁切多边形的有向投影面积及切向一阶矩。
m=zeros(1,3); if size(poly,1)<3, return; end
anchor=poly(1,:);
for k=2:size(poly,1)-1
    b=poly(k,:); c=poly(k+1,:); av=cross(b-anchor,c-anchor)/2;
    area=av(a); m=m+area*[1 (anchor(uv)+b(uv)+c(uv))/3];
end
end

function out=clip(poly,a,s,lower)
% 中文模块说明：原三角面逐半空间精确线性裁切；裁轴显式置于网格面。
if isempty(poly), out=zeros(0,3); return; end
out=zeros(size(poly,1)+2,3); used=0; prev=poly(end,:); d0=prev(a)-s;
if lower, in0=d0>=0; else, in0=d0<=0; end
for k=1:size(poly,1)
    curr=poly(k,:); d1=curr(a)-s; if lower, in1=d1>=0; else, in1=d1<=0; end
    if in0~=in1
        p=prev+d0/(d0-d1)*(curr-prev); p(a)=s;
        if used==0 || any(out(used,:)~=p), used=used+1; out(used,:)=p; end
    end
    if in1 && (used==0 || any(out(used,:)~=curr)), used=used+1; out(used,:)=curr; end
    prev=curr; d0=d1; in0=in1;
end
out=out(1:used,:); if size(out,1)>1 && all(out(1,:)==out(end,:)), out(end,:)=[]; end
end

function p=planeLine(poly,a,s)
% 中文模块说明：截面的正长度线段只用于独立整面拓扑识别，不补造面积。
p=zeros(0,3); prev=poly(end,:); d0=prev(a)-s;
for k=1:size(poly,1)
    curr=poly(k,:); d1=curr(a)-s;
    if d1==0, p(end+1,:)=curr; end %#ok<AGROW>
    if d0*d1<0, v=prev+d0/(d0-d1)*(curr-prev); v(a)=s; p(end+1,:)=v; end %#ok<AGROW>
    prev=curr; d0=d1;
end
if ~isempty(p), p=unique(p,'rows','stable'); end
end

function inside=pointInside(p,V,F)
% 中文模块说明：独立三角面立体角绕数；不依赖体积正负或阈值裁体积。
a=V(F(:,1),:)-p; b=V(F(:,2),:)-p; c=V(F(:,3),:)-p;
la=sqrt(sum(a.^2,2)); lb=sqrt(sum(b.^2,2)); lc=sqrt(sum(c.^2,2));
num=sum(a.*cross(b,c,2),2);
den=la.*lb.*lc+sum(a.*b,2).*lc+sum(b.*c,2).*la+sum(c.*a,2).*lb;
w=sum(2*atan2(num,den))/(4*pi);
if abs(w)<1e-8, inside=false;
elseif abs(w-1)<1e-8, inside=true;
else, error('CutCell:AmbiguousWinding','独立绕数不在0/1附近 %.17g point=%s.',w,mat2str(p,17)); end
end

function validateClosedOrientation(V,F)
% 中文模块说明：闭合且共享边方向相反是有向帽积分的前提。
[~,~,map]=unique(V,'rows'); f=reshape(map(F),size(F));
e=[f(:,[1 2]); f(:,[2 3]); f(:,[3 1])]; sorted=sort(e,2); [~,~,j]=unique(sorted,'rows');
cnt=accumarray(j,1); sign=2*(e(:,1)<e(:,2))-1; balance=accumarray(j,sign);
if any(cnt~=2) || any(balance~=0), error('CutCell:NotClosedOriented','曲面须为无开放边、共享边反向的闭合二维流形。'); end
r=mean(V,1); a=V(F(:,1),:)-r; b=V(F(:,2),:)-r; c=V(F(:,3),:)-r;
sv=sum(sum(a.*cross(b,c,2),2))/6;
if ~(isfinite(sv) && sv>0), error('CutCell:InwardOrientation','实体总体积必须为正，不能自动翻转或裁体积。'); end
end
