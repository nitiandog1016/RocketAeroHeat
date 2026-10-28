function bounds=sphericalTrianglePolarBounds(V,F,origin)
%SPHERICALTRIANGLEPOLARBOUNDS 求真实平面三角形相对+X的完整极角界。
% 模块用途：防止只筛顶点/质心而让三角形的边内部跨过条件入口。
% 正X半空间中，将顶点投影到(y/x,z/x)平面；原三角形的角方向恰为
% 投影三角形。最小极角来自其到原点的最短距离，最大值在顶点。
% 不完全位于正X半空间的面返回保守[0,180]，供调用者拒绝，不能外推。
validateattributes(V,{'numeric'},{'real','finite','ncols',3,'nonempty'});
validateattributes(F,{'numeric'},{'real','finite','integer','positive','ncols',3,'<=',size(V,1)});
validateattributes(origin,{'numeric'},{'real','finite','size',[1 3]});
r=V-origin;bounds=repmat([0 180],size(F,1),1);x=r(:,1);
front=all(reshape(x(F(:)),size(F))>0,2);
for j=find(front)'
    p=r(F(j,:),2:3)./r(F(j,:),1);next=[2 3 1];d=p(next,:)-p;
    denominator=sum(d.^2,2);t=zeros(3,1);active=denominator>0;
    t(active)=max(0,min(1,-sum(p(active,:).*d(active,:),2)./denominator(active)));
    closest=p+t.*d;minimumSquared=min(sum(closest.^2,2));
    signed=cross2(d,-p);area=cross2(p(2,:)-p(1,:),p(3,:)-p(1,:));
    if (area>0&&all(signed>=0))||(area<0&&all(signed<=0)),minimumSquared=0;end
    bounds(j,:)=[atand(sqrt(max(0,minimumSquared))),atand(max(vecnorm(p,2,2)))];
end
end

function value=cross2(a,b)
value=a(:,1).*b(:,2)-a(:,2).*b(:,1);
end
