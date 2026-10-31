function volume=buildSU2RocketVolume(surface,truth,options)
%BUILDSU2ROCKETVOLUME 将已核验完整箭体母线外域旋转为SU2三维体网格。
% 模块一：复用原生完整箭体及已有星形外域，拒绝突起/任意STL而不改形。
% 模块二：轴节点唯一共享，轴邻近用三棱柱，其余用六面体；没有轴向裂缝。
% 模块三：核对每个分解四面体正体积、内部面共享及完整边界账。
% 模块四：墙面保留原几何与父面候选；体网格通过不赋予壁热或三维物理资格。
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'regionalFlow','axisymmetric'));
if nargin<3,options=struct;end
[fv,domain]=buildRocketAxisymmetricDomain(surface,truth,options);
ring=surface.param.sectionNodeIDs{1};theta=mod(atan2(surface.V(ring,3),surface.V(ring,2)),2*pi);
theta=sort(theta(:));nt=numel(theta);assert(nt>=8 && all(diff([theta;theta(1)+2*pi])>0));
assert(max(abs(diff([theta;theta(1)+2*pi])-2*pi/nt))<1e-12, ...
    'RocketAeroHeat:SU2:Angles','当前旋转域要求原生等角周向分布。');
[ni,nj]=size(fv.X);ids=zeros(ni,nj,nt);V=zeros(((ni-2)*nt+2)*nj,3);cursor=0;
for j=1:nj
    for i=1:ni
        if i==1 || i==ni
            cursor=cursor+1;V(cursor,:)=[fv.X(i,j),0,0];ids(i,j,:)=cursor;
        else
            new=cursor+(1:nt);cursor=cursor+nt;ids(i,j,:)=reshape(new,1,1,nt);
            V(new,:)=[repmat(fv.X(i,j),nt,1),fv.R(i,j)*cos(theta),fv.R(i,j)*sin(theta)];
        end
    end
end
assert(cursor==size(V,1));cells=cell(fv.cellCount*nt,1);types=zeros(numel(cells),1);k=0;
volumeM3=zeros(numel(cells),1);minimumTet=inf;
for angle=1:nt
    next=mod(angle,nt)+1;a=ids(:,:,angle);b=ids(:,:,next);
    for c=1:fv.cellCount
        k=k+1;lower=a(fv.cells(c,:));upper=b(fv.cells(c,:));h=[lower(:)' upper(:)'];
        if h(1)==h(5) && h(2)==h(6)
            h=h([1 4 8 2 3 7]);types(k)=13;
        elseif h(3)==h(7) && h(4)==h(8)
            h=h([4 5 1 3 6 2]);types(k)=13;
        else,types(k)=12;end
        assert(numel(unique(h))==numel(h));cells{k}=h;
        if types(k)==12,tets=[1 2 3 7;1 3 4 7;1 4 8 7;1 8 5 7;1 5 6 7;1 6 2 7];
        else,tets=[1 2 3 4;2 3 4 5;3 4 5 6];end
        q=V(h,:);p=q(tets(:,1),:);v=dot(q(tets(:,2),:)-p,cross(q(tets(:,3),:)-p,q(tets(:,4),:)-p,2),2)/6;
        assert(all(isfinite(v)&v>0),'RocketAeroHeat:SU2:Volume','三维体单元有非正分解体积。');
        volumeM3(k)=sum(v);minimumTet=min(minimumTet,min(v));
    end
end
wall=cell((ni-1)*nt,1);outer=wall;parentBand=zeros(numel(wall),1);regions=cell(numel(wall),1);k=0;
for angle=1:nt
    next=mod(angle,nt)+1;
    for i=1:ni-1
        k=k+1;wall{k}=unique([ids(i,1,angle),ids(i+1,1,angle),ids(i+1,1,next),ids(i,1,next)],'stable');
        outer{k}=unique([ids(i,nj,angle),ids(i+1,nj,angle),ids(i+1,nj,next),ids(i,nj,next)],'stable');
        parentBand(k)=domain.wallSourceFaceBandID(i);regions(k)=domain.wallRegion(i);
    end
end
% bodyWall朝向实体外侧；流体域wall边界须反向，farfield保持外向。
bodyWall=wall;wall=cellfun(@(f)f([1 end:-1:2]),wall,'UniformOutput',false);
[bodyArea,bodyVolume]=surfaceIntegral(V,bodyWall);[~,outerVolume]=surfaceIntegral(V,outer);
assert(abs(bodyArea-sum(surface.area))/sum(surface.area)<2e-12, ...
    'RocketAeroHeat:SU2:Surface','体网格壁面改变了原表面面积。');
assert(abs(bodyVolume-domain.sourceTriangulatedBodyVolumeM3)/bodyVolume<2e-12);
assert(abs(sum(volumeM3)-(outerVolume-bodyVolume))/sum(volumeM3)<2e-12);
topology=checkFaces(cells,types,[wall;outer]);
volume=struct('kind','complete-rocket-three-dimensional-su2-volume','vertices',V,'cells',{cells}, ...
    'cellTypes',types,'cellVolume',volumeM3,'wall',{wall},'bodyWall',{bodyWall},'farfield',{outer}, ...
    'wallRegion',{regions},'wallSourceBand',parentBand,'sourceVertices',surface.V,'sourceFaces',surface.F, ...
    'domainOptions',domain.options,'sourceKind',truth.kind,'azimuthRad',theta, ...
    'bodyArea',bodyArea,'bodyVolume',bodyVolume,'fluidVolume',sum(volumeM3), ...
    'minimumSubTetVolume',minimumTet,'topology',topology, ...
    'sourceMappingStatus','original-band-candidates; field-projection-not-yet-implemented', ...
    'nearWallResolutionValidated',false,'farfieldAdequacyValidated',false,'physicalSolution',false, ...
    'productionQualified',false,'newQualifiedRocketArea',0);
end

function [area,volume]=surfaceIntegral(V,faces)
% 有向面扇形三角化用于独立面积与散度体积账，不给流场补值。
area=0;volume=0;
for k=1:numel(faces)
    f=faces{k};a=V(f(1),:);
    for j=2:numel(f)-1
        b=V(f(j),:);c=V(f(j+1),:);area=area+norm(cross(b-a,c-a))/2;
        volume=volume+dot(a,cross(b,c))/6;
    end
end
end

function info=checkFaces(cells,types,boundary)
% 对全部体单元面计数，核对内面恰好共享两次且外面与声明边界完全一致。
H=[1 4 3 2;5 6 7 8;1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8];
W={ [1 3 2],[4 5 6],[1 2 5 4],[2 3 6 5],[3 1 4 6] };
keys=zeros(6*numel(cells),4);n=0;
for k=1:numel(cells)
    if types(k)==12,local=mat2cell(H,ones(6,1),4);else,local=W;end
    for j=1:numel(local)
        f=sort(cells{k}(local{j}));n=n+1;keys(n,1:numel(f))=f;
    end
end
keys=sort(keys(1:n,:),2);[uniqueFaces,~,groups]=unique(keys,'rows');count=accumarray(groups,1);
assert(all(count<=2),'RocketAeroHeat:SU2:Topology','体网格非流形面。');
declared=zeros(numel(boundary),4);
for k=1:numel(boundary),declared(k,1:numel(boundary{k}))=boundary{k};end
declared=sortrows(sort(declared,2));actual=sortrows(uniqueFaces(count==1,:));
assert(isequal(declared,actual),'RocketAeroHeat:SU2:Topology','声明边界与实际单元边界不一致。');
info=struct('internalFaces',sum(count==2),'boundaryFaces',sum(count==1), ...
    'nonManifoldFaces',sum(count>2),'declaredBoundaryExactlyMatched',true);
end
