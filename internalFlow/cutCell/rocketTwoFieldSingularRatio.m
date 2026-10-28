function ratio=rocketTwoFieldSingularRatio(pages)
%ROCKETTWOFIELDSINGULARRATIO 稳定计算实2×2块的最小/最大奇异值比，仅用于诊断。
% 模块一：先逐块缩放，避免宽尺度系数平方溢出或下溢。
% 模块二：用两个hypot构造最大奇异值，避免近等奇异值时判别式相消。
% 模块三：零块返回比值0；该局部诊断不判断整个稀疏系统的秩。
assert(isreal(pages)&&size(pages,1)==2&&size(pages,2)==2&&all(isfinite(pages),'all'),'FieldRatio:Input','需要有限实2量块。');
scale=max(max(abs(pages),[],1),[],2);scale(scale==0)=1;scaled=pages./scale;
a=reshape(scaled(1,1,:),[],1);b=reshape(scaled(1,2,:),[],1);c=reshape(scaled(2,1,:),[],1);d=reshape(scaled(2,2,:),[],1);
maximum=(hypot(a+d,c-b)+hypot(a-d,b+c))/2;
ratio=abs(a.*d-b.*c)./max(realmin,maximum.^2);
assert(all(isfinite(ratio)),'FieldRatio:Nonfinite','局部奇异值比非有限。');
end
