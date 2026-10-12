function testCanonicalBody()
% 模块用途：核对闭合性、边方向、球/锥解析面积体积与网格收敛。
addpath(fileparts(fileparts(mfilename('fullpath'))));
kinds={'sphere','cone','coneCylinder','hemisphereCylinder','bluntedCone'};
for j=1:numel(kinds)
    errors=zeros(2,2);
    for k=1:2
        [m,meta]=buildCanonicalBody(kinds{j},struct('resolution',24*2^(k-1)));
        F=m.F;V=m.V;edges=[F(:,[1 2]);F(:,[2 3]);F(:,[3 1])];
        assert(all(diff(m.param.sectionX)>0) && all(m.param.sectionEquivalentRadius>0));
        assert(numel(m.param.faceBandID)==size(F,1));
        [~,~,id]=unique(sort(edges,2),'rows');
        assert(all(accumarray(id,1)==2),'存在开口/非流形边');
        assert(all(accumarray(id,sign(edges(:,2)-edges(:,1)))==0),'相邻面朝向不一致');
        assert(size(V,1)-max(id)+size(F,1)==2,'闭合体 Euler 示性数错误');
        a=cross(V(F(:,2),:)-V(F(:,1),:),V(F(:,3),:)-V(F(:,1),:),2);
        area=vecnorm(a,2,2)/2; assert(all(area>0));
        volume=sum(dot(V(F(:,1),:),cross(V(F(:,2),:),V(F(:,3),:),2),2))/6;
        assert(volume>0 && numel(meta.faceRegion)==size(F,1));
        errors(k,:)=[abs(sum(area)/meta.exactArea-1),abs(volume/meta.exactVolume-1)];
        if k==2,assert(all(errors(k,:)<0.015));end
    end
    assert(all(errors(2,:)<errors(1,:)),'细化未改善解析几何误差');
    fprintf('%s: fine area/volume relative error %.4g / %.4g\n',kinds{j},errors(2,:));
end
fprintf('Canonical body topology and geometric convergence passed.\n');
end
