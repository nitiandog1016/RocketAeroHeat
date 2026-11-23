function bands=thermistorWireBands(vertices,faces,wall,profile,airK,rootK)
% 模块一：沿两根名义直引线，从线根到截点按原面真实轴向边界分段。
% 以面元所属的截段为积分范围，保留入壁热流的原符号，报表中的功率以向空气放热为正。
isSide=ismember(wall.regionID,[2 6]);x=reshape(abs(vertices(faces(:),1)),size(faces))-profile.rootXM;
cuts=uniquetol(x(isSide,:),1e-13,'DataScale',1);cuts=sort(cuts(:));n=numel(cuts)-1;
assert(n>0 && abs(cuts(1))<1e-13 && abs(cuts(end)-profile.lengthM)<1e-13);
cuts(1)=0;cuts(end)=profile.lengthM;mid=mean(x,2);mu=profile.muPerM;
bands=struct('fromM',{},'toM',{},'rootwardK',{},'tipwardK',{},'theorySideW',{},'backSideW',{},'frontSideW',{});
total=zeros(1,2);
for j=1:n
    lower=cuts(j);upper=cuts(j+1);use=isSide & mid>=lower & mid<upper;
    assert(all(all(x(use,:)>=lower-1e-13 & x(use,:)<=upper+1e-13)));
    % 模块二：从一维模型两端轴向热流之差及三维原面壁反力分别取得该段放热。
    power=zeros(1,2);regions=[2 6];
    for k=1:2
        mask=use & wall.regionID==regions(k);assert(any(mask));
        power(k)=-sum(wall.areaM2(mask).*wall.heatFluxWm2(mask));total(k)=total(k)+nnz(mask);
    end
    bands(j)=struct('fromM',lower,'toM',upper,'rootwardK',airK+(rootK-airK)*exp(-mu*lower), ...
        'tipwardK',airK+(rootK-airK)*exp(-mu*upper),'theorySideW',profile.oneLeadRootW*(exp(-mu*lower)-exp(-mu*upper)), ...
        'backSideW',power(1),'frontSideW',power(2)); %#ok<AGROW>
end
assert(isequal(total,[nnz(wall.regionID==2) nnz(wall.regionID==6)]) && abs(sum([bands.theorySideW])-profile.oneLeadSideW)<1e-12);
end
