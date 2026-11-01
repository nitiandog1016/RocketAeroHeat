function [volume,info]=su2PyramidGeometry(V,Q)
%SU2PYRAMIDGEOMETRY 五节点金字塔全参考方向、体积及中位面闭合检查。
% 模块一：除去顶点自然退化因子(1-t)^2，双线性Jacobian的四角确定全域正性。
% 模块二：独立用真实面中心三角扇积分，核对体积并拒绝负中位子体积。
if isempty(Q)
    volume=zeros(0,1);info=struct('count',0,'minimumScaledJacobian',NaN,'minimumMedianSubVolume',NaN);return;
end
assert(size(Q,2)==5 && all(Q(:)>=1) && all(Q(:)<=size(V,1)));
p0=V(Q(:,1),:);a=V(Q(:,2),:)-p0;b=V(Q(:,4),:)-p0;
c=V(Q(:,3),:)-V(Q(:,2),:)-V(Q(:,4),:)+p0;apex=V(Q(:,5),:);
corners=[0 0;1 0;1 1;0 1];jac=zeros(size(Q,1),4);
for k=1:4
    r=corners(k,1);s=corners(k,2);point=p0+r*a+s*b+r*s*c;
    jac(:,k)=dot(cross(a+s*c,b+r*c,2),apex-point,2);
end
assert(all(jac>0,'all'),'金字塔全参考归一化Jacobian非正。');volume=sum(jac,2)/12;
center=zeros(size(p0));for k=1:5,center=center+V(Q(:,k),:)/5;end
faces={[1 4 3 2],[1 2 5],[2 3 5],[3 4 5],[4 1 5]};integral=zeros(size(volume));closure=zeros(size(center));area=zeros(size(volume));smallest=inf;
for j=1:5
    ids=Q(:,faces{j});fc=zeros(size(center));for k=1:size(ids,2),fc=fc+V(ids(:,k),:)/size(ids,2);end
    fc=fc-center;
    for k=1:size(ids,2)
        x=V(ids(:,k),:)-center;y=V(ids(:,mod(k,size(ids,2))+1),:)-center;
        av=cross(y-x,fc-x,2)/2;sub=dot(x,av,2)/3;
        assert(all(sub>0),'金字塔中位子体积非正。');smallest=min(smallest,min(sub));
        integral=integral+sub;closure=closure+av;area=area+vecnorm(av,2,2);
    end
end
assert(all(abs(integral-volume)<=2e-12*volume) && all(vecnorm(closure,2,2)<=2e-12*area));
info=struct('count',size(Q,1),'minimumScaledJacobian',min(jac,[],'all'), ...
    'minimumMedianSubVolume',smallest,'maximumVolumeFormulaRelativeError',max(abs(integral-volume)./volume));
end
