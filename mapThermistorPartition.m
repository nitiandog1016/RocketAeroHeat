function map=mapThermistorPartition(mesh,body,tol)
% 模块一：先识别保留的整面，余下的子面只允许落在一个原三角形上。
v=mesh.V;f=mesh.F;r=body.V;rf=body.F;parents=size(rf,1);
offset=(max(v)+min(v)-max(r)-min(r))/2;v=v-offset;
[ok,idx]=ismembertol(v,r,tol,'ByRows',true,'DataScale',1,'OutputAllIndices',true);
ids=zeros(size(ok));knownVertex=ok & cellfun(@numel,idx)==1;ids(knownVertex)=cell2mat(idx(knownVertex));
[whole,order]=ismember(sort(ids(f),2),sort(rf,2),'rows');
used=accumarray(order(whole),1,[parents 1]);duplicate=find(used>1,1);
if ~isempty(duplicate),error('RAH:Bead:Duplicate','标准面 %d 被多个整面覆盖。',duplicate);end
remaining=find(~used);children=find(~whole);bary=zeros(numel(children),9);
a=r(rf(:,1),:);b=r(rf(:,2),:);c=r(rf(:,3),:);
lo=min(min(a,b),c)-tol;hi=max(max(a,b),c)+tol;
refNormal=cross(b-a,c-a,2);parentArea=vecnorm(refNormal,2,2)/2;
partArea=zeros(size(order));partArea(whole)=parentArea(order(whole));
actual=cross(v(f(:,2),:)-v(f(:,1),:),v(f(:,3),:)-v(f(:,1),:),2);
actualArea=vecnorm(actual,2,2)/2;
distance=vecnorm(v(knownVertex,:)-r(ids(knownVertex),:),2,2);maxError=max([0;distance]);
assert(maxError<=tol,'RAH:Bead:Vertex','已识别角点超出几何容差。');
for k=1:numel(children)
    i=children(k);points=v(f(i,:),:);
    candidates=remaining(all(lo(remaining,:)<=min(points),2) & all(hi(remaining,:)>=max(points),2));
    for j=candidates'
        tri=r(rf(j,:),:);basis=(tri(2:3,:)-tri(1,:))';uv=(basis\(points-tri(1,:))')';
        w=[1-sum(uv,2),uv];edge=vecnorm(tri([2 3 1],:)-tri([3 1 2],:),2,2);
        altitude=2*parentArea(j)./edge';
        % 在角点容差内贴合父边。保留内部点，不将跨边三角形投给最近面。
        if any(w.*altitude < -tol,'all'),continue;end
        w(abs(w.*altitude)<=tol)=0;w=w./sum(w,2);
        err=max(vecnorm(w*tri-points,2,2));
        if ~isfinite(err) || err>tol,continue;end
        assert(order(i)==0,'RAH:Bead:Ambiguous','导入面 %d 同时对应多个标准面。',i);
        order(i)=j;bary(k,:)=reshape(w',1,9);maxError=max(maxError,err);
        partArea(i)=det(w)*parentArea(j);
    end
    assert(order(i)>0,'RAH:Bead:Triangle','导入面 %d 没有唯一的共面标准面，或与完整面重叠。',i);
end
cosine=sum((actual./(2*actualArea)).*(refNormal(order,:)./(2*parentArea(order))),2);
bad=find(~(cosine>.999999 & partArea>0 & abs(actualArea./partArea-1)<=1e-5),1);
assert(isempty(bad),'RAH:Bead:NormalArea','导入面 %d 的绕序、平面或面积不符。',bad);
% 模块二：在父面重心坐标平面核对并集面积以及子面面积和，分别识别缺口和重叠。
for j=remaining'
    use=find(order(children)==j);cover=polyshape;
    for k=use'
        w=reshape(bary(k,:),3,3)';cover=union(cover,polyshape(w(:,2),w(:,3)));
    end
    fraction=sum(partArea(children(use)))/parentArea(j);
    assert(abs(fraction-1)<1e-10 && abs(2*area(cover)-1)<1e-10, ...
        'RAH:Bead:Coverage','标准面 %d 的重剖三角形存在缺口或重叠。',j);
end
map=struct('referenceFaceID',order,'translationM',offset,'areaM2',actualArea,'referenceAreaM2',partArea, ...
    'maxVertexErrorM',maxError,'minNormalCosine',min(cosine),'maxRelativeAreaError',max(abs(actualArea./partArea-1)), ...
    'childIDs',children,'childBarycentric',bary,'parentAreaM2',parentArea(order));
end
